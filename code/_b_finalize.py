#!/usr/bin/env python
"""Track A finalizer: add temporal_split/document_year/box_clipped to score tables + metadata.

Per target:
  <target>_scores.csv  -> adds columns: exhaustiveness, temporal_split, document_year, box_clipped
  ligand_split.csv     -> sidecar (ligand_id, temporal_split, document_year) for Amendment #9
  results/trackA/full/target_metadata.csv (upsert row: target, box_size_A, decoy_ratio,
      sampling_seed, decoy_subsample_sha256, exhaustiveness, panel_role, n_act, n_dec, n_rows, status)

Run: tools\\envs\\boltz\\Scripts\\python.exe code\\trackA\\_b_finalize.py --target LpxC_Ec
     ... --all   # regenerate target_metadata.csv for all 9 targets (placeholders where pending)
"""
import argparse, csv, json, os, sys

ROOT = r"D:\harness\projects\d1_gramneg_vs"
OUTROOT = os.path.join(ROOT, "results", "trackA", "full")
SPLIT_DIR = os.path.join(ROOT, "data", "trackA", "splits")
CLIPPED = os.path.join(ROOT, "data", "trackA", "box_clipped_ligands.csv")
TARGETS = ["LpxC_Ec", "LpxC_Pa", "DHFR_Ec", "Gyrase_Ec", "GyrA_Ec", "GyrB_Ec", "Fabi_Ec", "FabI_Ec", "TEM1_Ec", "LepB_Ec"]
BOX_SIZE_A = {"LpxC_Ec": 26, "GyrB_Ec": 26, "Gyrase_Ec": 26, "LepB_Ec": 26}
ROLE = {"LpxC_Ec": "main", "LpxC_Pa": "main", "DHFR_Ec": "main", "Gyrase_Ec": "main",
        "GyrA_Ec": "low_confidence_auxiliary", "GyrB_Ec": "main", "FabI_Ec": "main",
        "TEM1_Ec": "main", "LepB_Ec": "auxiliary"}

def read_splits(target):
    """inchikey -> (split, document_year)."""
    for fn in (f"{target}_temporal_latest20.csv",):
        p = os.path.join(SPLIT_DIR, fn)
        if not os.path.exists(p):
            continue
        m = {}
        for r in csv.DictReader(open(p, encoding="utf-8-sig")):
            ik = (r.get("inchikey") or "").strip()
            if not ik:
                continue
            m[ik] = (r.get("split_latest20") or r.get("split") or "", r.get("document_year") or "")
        return m
    return {}

def clipped_set(target):
    s = set()
    if os.path.exists(CLIPPED):
        for r in csv.DictReader(open(CLIPPED, encoding="utf-8-sig")):
            if r.get("target", "").strip() == target and r.get("inchikey"):
                s.add(r["inchikey"].strip())
    return s

def read_manifest(target):
    p = os.path.join(OUTROOT, target, "manifest.csv")
    if not os.path.exists(p):
        return {}
    return {r["ligand_id"]: r for r in csv.DictReader(open(p, encoding="utf-8-sig"))}

def read_splits_map(target):
    """Two maps: chembl_id -> (split, year) and inchikey -> (split, year)."""
    by_chembl, by_ik = {}, {}
    import glob
    cands = [f"{target}_temporal_latest20.csv"]
    for fn in cands:
        p = os.path.join(SPLIT_DIR, fn)
        if not os.path.exists(p):
            continue
        for r in csv.DictReader(open(p, encoding="utf-8-sig")):
            sp = r.get("split_latest20") or r.get("split") or ""
            yr = r.get("document_year") or ""
            cid = (r.get("compound_chembl_id") or "").strip()
            ik = (r.get("inchikey") or "").strip()
            if cid:
                by_chembl[cid] = (sp, yr)
            if ik:
                by_ik[ik] = (sp, yr)
    return by_chembl, by_ik

def finalize_target(target):
    d = os.path.join(OUTROOT, target)
    score_p = os.path.join(d, f"{target}_scores.csv")
    if not os.path.isdir(d) or not os.path.exists(score_p):
        print(f"[{target}] scores not aggregated yet, skip"); return None
    rows = list(csv.DictReader(open(score_p, encoding="utf-8-sig")))
    man = read_manifest(target)
    by_chembl, by_ik = read_splits_map(target)
    clipped = clipped_set(target)
    out, split_rows = [], []
    n_reused = 0
    for r in rows:
        lid = r["ligand_id"]
        m = man.get(lid, {})
        ik = (m.get("inchikey") or "").strip()
        raw_status = r.get("status", "")
        reused = ("reuse" in raw_status.lower()) or (str(r.get("reused", "")).lower() == "true")
        status = "ok" if raw_status.startswith("ok") else raw_status
        if reused:
            n_reused += 1
        is_active = str(r["label"]) == "1"
        tsplit, yr = "", ""
        if is_active:
            bare = lid[4:] if lid.startswith("act_") else lid
            tsplit, yr = by_chembl.get(bare, ("", ""))
            if not tsplit and ik:
                tsplit, yr = by_ik.get(ik, ("", ""))
        else:
            tsplit = "decoy"
        bc = "True" if (ik and ik in clipped) else m.get("box_clipped", "False")
        if str(bc).lower() != "true":
            bc = "False"
        out.append({"ligand_id": lid, "label": r["label"], "score": r.get("score", ""),
                    "status": status, "reused": "True" if reused else "False",
                    "seed": r.get("seed", 42), "exhaustiveness": r.get("exhaustiveness", ""),
                    "temporal_split": tsplit, "document_year": yr, "box_clipped": bc})
        if is_active:
            split_rows.append({"ligand_id": lid, "temporal_split": tsplit, "document_year": yr})
    cols = ["ligand_id", "label", "score", "status", "reused", "seed", "exhaustiveness",
            "temporal_split", "document_year", "box_clipped"]
    with open(score_p, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(out)
    with open(os.path.join(d, "ligand_split.csv"), "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=["ligand_id", "temporal_split", "document_year"]); w.writeheader(); w.writerows(split_rows)
    from collections import Counter
    ok_rows = [r for r in out if r["status"] == "ok"]
    n_act_ok = sum(1 for r in ok_rows if r["label"] == "1")
    n_dec_ok = sum(1 for r in ok_rows if r["label"] != "1")
    n_prep_missing = sum(1 for r in out if r["status"] == "prep_missing")
    n_box_clipped = sum(1 for r in out if r["box_clipped"] == "True")
    dist = Counter(r["temporal_split"] for r in split_rows)
    stats = {"n_rows": len(out), "n_ok": len(ok_rows), "n_fail": len(out) - len(ok_rows),
             "n_reused": n_reused, "n_prep_missing": n_prep_missing,
             "n_box_clipped": n_box_clipped, "n_act_ok": n_act_ok, "n_dec_ok": n_dec_ok,
             "decoy_ratio_observed": round(n_dec_ok / n_act_ok, 2) if n_act_ok else "",
             "temporal_split_actives": {"train": dist.get("train", 0), "test": dist.get("test", 0),
                                        "missing": dist.get("", 0)},
             "n_actives_with_split": sum(1 for r in split_rows if r["temporal_split"])}
    sp = os.path.join(d, "summary.json")
    if os.path.exists(sp):
        s = json.load(open(sp, encoding="utf-8"))
    else:
        s = {}
    s.update({k: stats[k] for k in ("n_prep_missing", "n_box_clipped", "decoy_ratio_observed",
                                    "temporal_split_actives", "n_reused")})
    s.update({"n_rows": stats["n_rows"], "n_ok": stats["n_ok"], "n_fail": stats["n_fail"],
              "n_act_ok": stats["n_act_ok"], "n_dec_ok": stats["n_dec_ok"]})
    s["finalize_date"] = "2026-11-17"
    json.dump(s, open(sp, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    print(json.dumps({target: stats}, indent=1), flush=True)
    return stats

def upsert_metadata(target, stats):
    p = os.path.join(OUTROOT, "target_metadata.csv")
    cols = ["target", "box_size_A", "decoy_ratio", "decoy_ratio_observed", "sampling_seed",
            "decoy_subsample_sha256", "exhaustiveness", "panel_role", "n_act", "n_dec", "n_rows",
            "n_ok", "n_prep_missing", "n_box_clipped", "temporal_train", "temporal_test", "status"]
    rows = []
    if os.path.exists(p):
        rows = list(csv.DictReader(open(p, encoding="utf-8-sig")))
    by = {r["target"]: r for r in rows}
    d = os.path.join(OUTROOT, target)
    meta = json.load(open(os.path.join(d, "sampling_meta.json"), encoding="utf-8")) if os.path.exists(os.path.join(d, "sampling_meta.json")) else {}
    summ = json.load(open(os.path.join(d, "summary.json"), encoding="utf-8")) if os.path.exists(os.path.join(d, "summary.json")) else {}
    ts = summ.get("temporal_split_actives", {}) if isinstance(summ.get("temporal_split_actives"), dict) else {}
    by[target] = {"target": target, "box_size_A": BOX_SIZE_A.get(target, 22),
                  "decoy_ratio": meta.get("decoy_ratio", ""),
                  "decoy_ratio_observed": summ.get("decoy_ratio_observed", ""),
                  "sampling_seed": meta.get("sampling_seed", ""),
                  "decoy_subsample_sha256": meta.get("decoy_subsample_sha256", ""),
                  "exhaustiveness": summ.get("exhaustiveness", ""),
                  "panel_role": ROLE.get(target, "main"),
                  "n_act": summ.get("n_act", meta.get("n_act", "")),
                  "n_dec": summ.get("n_dec", meta.get("n_dec_sampled", "")),
                  "n_rows": summ.get("n_rows", meta.get("n_rows", "")),
                  "n_ok": summ.get("n_ok", ""),
                  "n_prep_missing": summ.get("n_prep_missing", ""),
                  "n_box_clipped": summ.get("n_box_clipped", ""),
                  "temporal_train": ts.get("train", sys.__dict__.get("__empty__", "")) if ts else "",
                  "temporal_test": ts.get("test", "") if ts else "",
                  "status": "done" if summ else ("prepared" if meta else "pending")}
    order = ["LpxC_Ec","LpxC_Pa","DHFR_Ec","Gyrase_Ec","GyrA_Ec","GyrB_Ec","FabI_Ec","TEM1_Ec","LepB_Ec"]
    out = [by[t] for t in order if t in by]
    with open(p, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(out)
    print("wrote", p, len(out), "rows", flush=True)



ALL_TARGETS = ["LpxC_Ec", "GyrB_Ec", "Gyrase_Ec", "LpxC_Pa", "DHFR_Ec", "GyrA_Ec", "FabI_Ec", "TEM1_Ec", "LepB_Ec"]

def write_panel_ready():
    """When all 9 main-panel targets have summary.json, write PANEL_READY.json + notify message draft."""
    import datetime
    infos = {}
    all_done = True
    for t in ALL_TARGETS:
        sp = os.path.join(OUTROOT, t, "summary.json")
        if os.path.exists(sp):
            s = json.load(open(sp, encoding="utf-8"))
            infos[t] = {"status": "done", "n_rows": s.get("n_rows"), "n_act": s.get("n_act"),
                        "n_dec": s.get("n_dec"), "n_ok": s.get("n_ok"), "n_fail": s.get("n_fail"),
                        "n_prep_missing": s.get("n_prep_missing", s.get("n_fail")),
                        "n_box_clipped": s.get("n_box_clipped", 0),
                        "decoy_ratio_observed": s.get("decoy_ratio_observed", ""),
                        "exhaustiveness": s.get("exhaustiveness"), "box_size": s.get("box_size"),
                        "elapsed_s": s.get("elapsed_s"), "ligands_per_min": s.get("ligands_per_min")}
        else:
            infos[t] = {"status": "missing"}
            all_done = False
    if not all_done:
        missing = [k for k, v in infos.items() if v["status"] != "done"]
        print("[panel_ready] not yet; missing:", missing, flush=True)
        return False
    payload = {"ready": True, "tag": "full", "generated_at": datetime.datetime.now().isoformat(timespec="seconds"),
               "targets": infos, "n_targets": len(infos),
               "n_total_rows": sum(v["n_rows"] or 0 for v in infos.values()),
               "n_total_ok": sum(v["n_ok"] or 0 for v in infos.values()),
               "instructions": "D-line: run `python -m maes.tracka --tag full` (or maes.tracka --tag full); inputs at results/trackA/full/<target>/<target>_scores.csv + target_metadata.csv + ligand_split.csv",
               "source": "B-line _b_finalize.py write_panel_ready()"}
    out = os.path.join(OUTROOT, "PANEL_READY.json")
    json.dump(payload, open(out, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    msg = f"""PANEL_READY: Track A full panel 9/9 done ({payload['generated_at']})
total rows={payload['n_total_rows']} ok={payload['n_total_ok']}
inputs: results/trackA/full/<target>/<target>_scores.csv + target_metadata.csv + ligand_split.csv
next: D-line run `maes.tracka --tag full`"""
    open(os.path.join(OUTROOT, "PANEL_READY_MESSAGE.md"), "w", encoding="utf-8").write(msg + "\n")
    print("[panel_ready] WRITTEN", out, flush=True)
    return True

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--target")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    if a.all:
        for t in ["LpxC_Ec","LpxC_Pa","DHFR_Ec","Gyrase_Ec","GyrA_Ec","GyrB_Ec","FabI_Ec","TEM1_Ec","LepB_Ec"]:
            done = finalize_target(t)
            upsert_metadata(t, done)
        write_panel_ready()
    elif a.target:
        done = finalize_target(a.target)
        upsert_metadata(a.target, done)
        write_panel_ready()
    else:
        sys.exit("need --target or --all")
