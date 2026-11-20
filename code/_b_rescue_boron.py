#!/usr/bin/env python
"""Rescue boron-containing actives that fail Vina atom typing (Lead request 2026-11-18).

Vina 1.2.5 rejects AutoDock atom type 'B' (boron) -> all boronic-acid actives fail
(TEM1_Ec 23/91, LepB_Ec 2/27, both 100% boron). Rescue: rewrite ligand PDBQT atom type
B -> C (isosteric surrogate; documented), redock with the SAME settings, update scores CSV.

Outputs per target:
  <target>/rescue_boron_surrogate.csv  (ligand_id,label,original_status,surrogate,score,note)
  <target>_scores.csv updated (status ok, score filled) for rescued rows only
Run: tools\\envs\\boltz\\Scripts\\python.exe code\\trackA\\_b_rescue_boron.py
"""
import csv, json, os, re, subprocess, sys

ROOT = r"D:\harness\projects\d1_gramneg_vs"
OUTROOT = os.path.join(ROOT, "results", "trackA", "full")
VINA = os.path.join(ROOT, "tools", "vina", "vina.exe")
REC = {"TEM1_Ec": "TEM1_Ec_holo.pdbqt", "LepB_Ec": "LepB_Ec_holo.pdbqt"}
BOX = {"TEM1_Ec": [3.088, 13.906, 18.215], "LepB_Ec": [64.921, 22.645, 10.289]}
SIZE = {"TEM1_Ec": 22, "LepB_Ec": 26}
EXH = {"TEM1_Ec": 2, "LepB_Ec": 2}
SEED = 42

def retype_b_to_c(src, dst):
    n = 0
    out = []
    for ln in open(src, encoding="utf-8", errors="replace"):
        if ln.startswith(("ATOM", "HETATM")):
            parts = ln.rstrip("\n").split()
            if parts and parts[-1] == "B":
                # rebuild line preserving columns except final type token replaced
                idx = ln.rfind(" B")
                if idx > 0:
                    ln = ln[:idx] + " C" + ln[idx+2:]
                    n += 1
            out.append(ln.rstrip("\n"))
        else:
            out.append(ln.rstrip("\n"))
    open(dst, "w", encoding="utf-8").write("\n".join(out) + "\n")
    return n

def dock(lig, rec, box, size, exh, logf, outp):
    cmd = [VINA, "--receptor", rec, "--ligand", lig,
           "--center_x", str(box[0]), "--center_y", str(box[1]), "--center_z", str(box[2]),
           "--size_x", str(size), "--size_y", str(size), "--size_z", str(size),
           "--exhaustiveness", str(exh), "--seed", str(SEED), "--num_modes", "3", "--cpu", "8", "--out", outp]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=3600)
        txt = (r.stdout or "") + (r.stderr or "")
    except subprocess.TimeoutExpired:
        txt = "TIMEOUT"
    open(logf, "w", encoding="utf-8", errors="replace").write(txt)
    for ln in txt.splitlines():
        m = re.match(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s*$", ln)
        if m:
            return float(m.group(2))
    return None

def main():
    summary = {}
    for target in ["TEM1_Ec", "LepB_Ec"]:
        d = os.path.join(OUTROOT, target)
        sp = os.path.join(d, f"{target}_scores.csv")
        rows = list(csv.DictReader(open(sp, encoding="utf-8-sig")))
        man = {r["ligand_id"]: r for r in csv.DictReader(open(os.path.join(d, "manifest.csv"), encoding="utf-8-sig"))}
        todo = [r for r in rows if r["status"] == "fail" and str(r["label"]) == "1"]
        print(f"[{target}] failed actives to rescue: {len(todo)}", flush=True)
        rec = os.path.join(ROOT, "data", "receptors", REC[target])
        res = []
        for r in todo:
            lid = r["ligand_id"]
            m = man.get(lid, {})
            src = m.get("pdbqt", "")
            if not src or not os.path.exists(src):
                res.append({"ligand_id": lid, "label": r["label"], "original_status": r["status"],
                            "surrogate": "", "score": "", "note": "source_pdbqt_missing"})
                continue
            sur = os.path.join(d, "ligands", lid + "_B2C.pdbqt")
            nB = retype_b_to_c(src, sur)
            outp = os.path.join(d, "poses", lid + "_B2C_out.pdbqt")
            logf = os.path.join(d, "logs", "rescue_" + lid + ".log")
            score = dock(sur, rec, BOX[target], SIZE[target], EXH[target], logf, outp)
            note = f"B atoms={nB}; atom type B->C isosteric surrogate"
            res.append({"ligand_id": lid, "label": r["label"], "original_status": r["status"],
                        "surrogate": "B->C", "score": "" if score is None else score, "note": note})
            if score is not None:
                r["status"] = "ok"
                r["score"] = score
            else:
                r["status"] = "fail"
            print(f"  {lid}: B={nB} score={score}", flush=True)
        cols = list(rows[0].keys())
        with open(sp, "w", newline="", encoding="utf-8-sig") as f:
            w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)
        rp = os.path.join(d, "rescue_boron_surrogate.csv")
        with open(rp, "w", newline="", encoding="utf-8-sig") as f:
            w = csv.DictWriter(f, fieldnames=["ligand_id", "label", "original_status", "surrogate", "score", "note"])
            w.writeheader(); w.writerows(res)
        n_ok = sum(1 for x in res if x["score"] != "")
        summary[target] = {"n_failed_actives": len(todo), "n_rescued": n_ok, "n_still_fail": len(todo) - n_ok,
                           "rescue_csv": rp}
        print(json.dumps({target: summary[target]}, indent=1), flush=True)
    json.dump(summary, open(os.path.join(OUTROOT, "rescue_boron_summary.json"), "w", encoding="utf-8"),
              indent=1, ensure_ascii=False)

if __name__ == "__main__":
    main()
