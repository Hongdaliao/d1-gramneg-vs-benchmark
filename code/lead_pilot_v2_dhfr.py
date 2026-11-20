# -*- coding: utf-8 -*-
"""DHFR non-metal control pilot (Lead). Run with the boltz venv python (rdkit+meeko).

1) SMILES -> 3D SDF (RDKit ETKDGv3 + MMFF)
2) SDF -> PDBQT (meeko mk_prepare_ligand.exe)
3) Vina docking (same settings as pilot v2)
4) metrics + combined 3-panel figure (lpxc_zn / lpxc_nozn / dhfr)
"""
import csv, json, math, os, subprocess, sys, time
import concurrent.futures as cf

ROOT = r"D:\harness\projects\d1_gramneg_vs"
VINA = os.path.join(ROOT, r"tools\vina\vina.exe")
MKPREP = os.path.join(ROOT, r"tools\envs\boltz\Scripts\mk_prepare_ligand.exe")
OUT = os.path.join(ROOT, r"results\trackA\pilot_v2")
LIGDIR = os.path.join(OUT, "dhfr_ligands")
os.makedirs(LIGDIR, exist_ok=True)
MAN = os.path.join(ROOT, r"data\receptors\DHFR_pilot_manifest.csv")
REC = os.path.join(ROOT, r"data\receptors\DHFR_Ec_holo_3DAU.pdbqt")
BOX = dict(cx=27.158, cy=9.763, cz=-17.934, sx=22.0, sy=22.0, sz=22.0)
EXHAUST, SEED, CPU, WORKERS, NUM_MODES = 8, 42, 6, 2, 3

from rdkit import Chem
from rdkit.Chem import AllChem


def prep_ligands():
    rows = []
    for r in csv.DictReader(open(MAN, encoding="utf-8-sig")):
        lid = r["compound_id"] + "_" + r["label"]
        pdbqt = os.path.join(LIGDIR, lid + ".pdbqt")
        if not os.path.exists(pdbqt):
            sdf = os.path.join(LIGDIR, lid + ".sdf")
            mol = Chem.MolFromSmiles(r["smiles"])
            if mol is None:
                continue
            mol = Chem.AddHs(mol)
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            if AllChem.EmbedMolecule(mol, params) != 0:
                continue
            try:
                AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
            except Exception:
                pass
            w = Chem.SDWriter(sdf)
            w.write(mol)
            w.close()
            res = subprocess.run([MKPREP, "-i", sdf, "-o", pdbqt], capture_output=True, text=True)
            if res.returncode != 0 or not os.path.exists(pdbqt):
                continue
        rows.append(dict(ligand_id=lid, label="1" if r["label"] == "active" else "0", pdbqt=pdbqt))
    return rows


def dock_one(lig):
    outp = os.path.join(OUT, "dhfr", "poses", os.path.basename(lig["pdbqt"]).replace(".pdbqt", "_out.pdbqt"))
    os.makedirs(os.path.dirname(outp), exist_ok=True)
    cmd = [VINA, "--receptor", REC, "--ligand", lig["pdbqt"],
           "--center_x", str(BOX["cx"]), "--center_y", str(BOX["cy"]), "--center_z", str(BOX["cz"]),
           "--size_x", str(BOX["sx"]), "--size_y", str(BOX["sy"]), "--size_z", str(BOX["sz"]),
           "--exhaustiveness", str(EXHAUST), "--seed", str(SEED), "--cpu", str(CPU),
           "--num_modes", str(NUM_MODES), "--out", outp]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=600)
    except subprocess.TimeoutExpired:
        return dict(ligand_id=lig["ligand_id"], label=lig["label"], vina_score=None, status="timeout")
    affs = []
    if os.path.exists(outp):
        for line in open(outp, encoding="utf-8", errors="replace"):
            if line.startswith("REMARK VINA RESULT:"):
                try:
                    affs.append(float(line.split()[3]))
                except (ValueError, IndexError):
                    pass
    best = min(affs) if affs else None
    status = "ok" if best is not None else ("fail" if r.returncode != 0 else "nodata")
    return dict(ligand_id=lig["ligand_id"], label=lig["label"], vina_score=best, status=status)


def metrics(rows):
    ok = [r for r in rows if r["status"] == "ok" and r["vina_score"] is not None]
    labels = [int(r["label"]) for r in ok]
    n, na = len(labels), sum(labels)
    rk = sorted((r["vina_score"], int(r["label"])) for r in ok)
    lab = [y for _, y in rk]
    def ef(f):
        k = max(1, int(round(f * n)))
        return (sum(lab[:k]) / na) / f if na else float("nan")
    pos = [i for i, y in enumerate(lab) if y]
    neg = [i for i, y in enumerate(lab) if not y]
    auc = (sum(1 for p in pos for q in neg if p < q) / (len(pos) * len(neg))) if pos and neg else float("nan")
    return dict(n=n, n_active=na, ef1=round(ef(0.01), 4), ef5=round(ef(0.05), 4),
                ef10=round(ef(0.10), 4), auc=round(auc, 4))


def main():
    rows = prep_ligands()
    print(f"[dhfr] ligands prepared: {len(rows)}", flush=True)
    t0 = time.time()
    res = []
    with cf.ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = [ex.submit(dock_one, r) for r in rows]
        for i, f in enumerate(cf.as_completed(futs), 1):
            res.append(f.result())
            if i % 25 == 0:
                print(f"[dhfr] {i}/{len(rows)} ok={sum(1 for x in res if x['status']=='ok')}", flush=True)
    res.sort(key=lambda r: r["ligand_id"])
    with open(os.path.join(OUT, "dhfr_scores.csv"), "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["ligand_id", "label", "vina_score", "status", "elapsed_s"])
        for r in res:
            w.writerow([r["ligand_id"], r["label"], "" if r["vina_score"] is None else r["vina_score"], r["status"], ""])
    m = metrics(res)
    print(f"[dhfr] done {round(time.time()-t0,1)}s metrics={m}", flush=True)

    # combined summary
    summary = json.load(open(os.path.join(OUT, "pilot_v2_metrics.json"), encoding="utf-8"))
    summary["dhfr"] = dict(m, receptor=os.path.basename(REC), elapsed_s=round(time.time() - t0, 1))
    summary["delta_zn_minus_nozn"] = dict(ef1=0.0, ef5=0.0,
                                         auc=round(summary["lpxc_zn"]["auc"] - summary["lpxc_nozn"]["auc"], 4))
    json.dump(summary, open(os.path.join(OUT, "pilot_v2_metrics.json"), "w", encoding="utf-8"), indent=2, ensure_ascii=False)
    with open(os.path.join(OUT, "pilot_v2_metrics.csv"), "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["run", "n", "n_active", "ef1", "ef5", "ef10", "auc", "receptor"])
        for tag in ("lpxc_zn", "lpxc_nozn", "dhfr"):
            s = summary[tag]
            w.writerow([tag, s["n"], s["n_active"], s["ef1"], s["ef5"], s["ef10"], s["auc"], s.get("receptor", "")])
    # figure from the three scores files
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, tag, title in zip(axes, ("lpxc_zn", "lpxc_nozn", "dhfr"),
                              ("LpxC_Ec Zn-retained", "LpxC_Ec Zn-stripped", "DHFR (non-metal control)")):
        p = os.path.join(OUT, f"{tag}_scores.csv")
        if not os.path.exists(p):
            continue
        ok = [r for r in csv.DictReader(open(p, encoding="utf-8-sig")) if r["status"] == "ok" and r["vina_score"]]
        pairs = sorted((float(r["vina_score"]), int(r["label"])) for r in ok)
        lab = [y for _, y in pairs]
        n = len(lab); tot = sum(lab); c = 0; cum = []
        for y in lab:
            c += y
            cum.append(c / tot if tot else 0)
        xs = [i / n for i in range(1, n + 1)]
        mm = summary[tag]
        ax.plot(xs, cum, lw=2)
        ax.plot([0, 1], [0, 1], ":", color="k", lw=0.8)
        ax.set_title(f"{title}\nAUC={mm['auc']:.3f}, EF5={mm['ef5']:.2f}", fontsize=10)
        ax.set_xlabel("fraction of ranked database")
        ax.set_ylabel("fraction of actives recovered")
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, r"figures\pilot_v2_zn_nozn_dhfr.png"), dpi=200)
    fig.savefig(os.path.join(ROOT, r"figures\pilot_v2_zn_nozn_dhfr.pdf"))
    print(json.dumps(summary, ensure_ascii=False, indent=1), flush=True)


if __name__ == "__main__":
    main()
