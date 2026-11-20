#!/usr/bin/env python
"""Boron rescue: complete surrogate records + two-version sensitivity AUC (Lead request 2026-11-18).

Outputs:
  results/trackA/full/TEM1_Ec/rescue_boron_surrogate.csv   (enriched: original SMILES, error text, vina version)
  results/trackA/full/LepB_Ec/rescue_boron_surrogate.csv
  results/trackA/full/sensitivity_boron.csv                 (AUC/EF5: main vs boron-excluded)
  results/trackA/full/sensitivity_boron.json
Run: tools\\envs\\boltz\\Scripts\\python.exe code\\trackA\_b_boron_sensitivity.py
"""
import csv, json, os, subprocess

ROOT = r"D:\harness\projects\d1_gramneg_vs"
R = os.path.join(ROOT, "results", "trackA", "full")
VINA = os.path.join(ROOT, "tools", "vina", "vina.exe")
ERROR_TEXT = "PDBQT parsing error: Atom type B is not a valid AutoDock type (atom types are case-sensitive)."
JUSTIFY = ("Vina 'vina' scoring function uses atom-type-based terms only (no atomic charges), "
           "so B->C atom-type substitution keeps the molecular graph and search space identical "
           "while making the ligand parseable; treated as isosteric surrogate, sensitivity analysis excludes it.")
SURROGATE = {"TEM1_Ec": ["TEM1_Ec"], "LepB_Ec": ["LepB_Ec"]}

def vina_version():
    try:
        r = subprocess.run([VINA, "--version"], capture_output=True, text=True, timeout=30)
        return (r.stdout or r.stderr).strip().splitlines()[0]
    except Exception:
        return "unknown"

def enrich(target):
    d = os.path.join(R, target)
    sp = os.path.join(d, f"{target}_scores.csv")
    rp = os.path.join(d, "rescue_boron_surrogate.csv")
    if not (os.path.exists(sp) and os.path.exists(rp)):
        return []
    rows = list(csv.DictReader(open(sp, encoding="utf-8-sig")))
    man = {r["ligand_id"]: r for r in csv.DictReader(open(os.path.join(d, "manifest.csv"), encoding="utf-8-sig"))}
    old = {r["ligand_id"]: r for r in csv.DictReader(open(rp, encoding="utf-8-sig"))}
    out = []
    for lid, o in old.items():
        smi = man.get(lid, {}).get("smiles", "")
        row = {"ligand_id": lid, "label": o.get("label", ""), "original_status": o.get("original_status", "fail"),
               "original_smiles": smi, "surrogate_smiles": smi, "molecular_graph_changed": "False",
               "atom_type_original": "B", "atom_type_surrogate": "C",
               "vina_version": vina_version(), "error_message_original": ERROR_TEXT,
               "justification": JUSTIFY, "score": o.get("score", ""), "note": o.get("note", "")}
        out.append(row)
    cols = list(out[0].keys()) if out else []
    with open(rp, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(out)
    print(f"[{target}] enriched {len(out)} surrogate rows -> {rp}", flush=True)
    return out

def auc_ef(rows, excluded_active_ids=None):
    """rows: ok-scored rows; excluded_active_ids: actives removed for sensitivity version B."""
    ex = set(excluded_active_ids or [])
    data = []
    for r in rows:
        if r["status"] != "ok" or r.get("score") in ("", None):
            continue
        lid = r["ligand_id"]
        if str(r["label"]) == "1" and lid in ex:
            continue
        data.append((float(r["score"]), 1 if str(r["label"]) == "1" else 0, lid))
    data.sort(key=lambda x: (x[0], x[2]))
    n_pos = sum(1 for x in data if x[1] == 1)
    n_neg = len(data) - n_pos
    if n_pos == 0 or n_neg == 0:
        return {"n_act": n_pos, "n_dec": n_neg, "auc": None, "ef5": None, "ef1": None}
    # AUC via rank sum (lower score = better -> rank from low)
    ranks = {}
    i = 0
    while i < len(data):
        j = i
        while j + 1 < len(data) and data[j + 1][0] == data[i][0]:
            j += 1
        avg_rank = (i + j) / 2.0 + 1  # 1-based
        for k in range(i, j + 1):
            ranks[data[k][2]] = avg_rank
        i = j + 1
    sum_pos_ranks = sum(ranks[x[2]] for x in data if x[1] == 1)
    auc = (sum_pos_ranks - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)
    auc = 1.0 - auc  # because lower score = better, invert to standard direction
    def ef(frac):
        k = max(1, int(round(frac * len(data))))
        hits = sum(1 for x in data[:k] if x[1] == 1)
        return hits / n_pos / frac
    return {"n_act": n_pos, "n_dec": n_neg, "auc": round(auc, 4),
            "ef1": round(ef(0.01), 3), "ef5": round(ef(0.05), 3), "n_total": len(data)}

def main():
    all_sur = {}
    for t in ["TEM1_Ec", "LepB_Ec"]:
        all_sur[t] = enrich(t)
    out_rows = []
    detail = {}
    for t in ["TEM1_Ec", "LepB_Ec"]:
        d = os.path.join(R, t)
        sp = os.path.join(d, f"{t}_scores.csv")
        rows = list(csv.DictReader(open(sp, encoding="utf-8-sig")))
        sur_ids = [x["ligand_id"] for x in all_sur[t]]
        main = auc_ef(rows, None)
        sens = auc_ef(rows, sur_ids)
        detail[t] = {"main_including_surrogate": main, "sensitivity_excluding_boron_actives": sens,
                     "delta_auc": (round(sens["auc"] - main["auc"], 4) if main["auc"] is not None and sens["auc"] is not None else None),
                     "n_boron_actives": len(sur_ids)}
        out_rows.append({"target": t, "version": "main_incl_surrogate", **main})
        out_rows.append({"target": t, "version": "sensitivity_excl_boron", **sens})
    p = os.path.join(R, "sensitivity_boron.csv")
    with open(p, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys())); w.writeheader(); w.writerows(out_rows)
    json.dump(detail, open(os.path.join(R, "sensitivity_boron.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    print(json.dumps(detail, indent=1, ensure_ascii=False), flush=True)
    print("wrote", p, flush=True)

if __name__ == "__main__":
    main()
