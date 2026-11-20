#!/usr/bin/env python
"""Prepare frozen inputs for Amendment #13 (DrugCLIP) and #14 (Uni-Dock).

PREP ONLY (no screening, no docking). Produces:
  results/trackA/method_axis/drugclip/inputs/library_806.csv (+ .lmdb built by prepare_mols.py)
  results/trackA/method_axis/drugclip/inputs/manifest.json
  results/trackA/method_axis/unidock/inputs/{target}.csv
  results/trackA/method_axis/unidock/inputs/unidock_inputs_all.csv
  results/trackA/method_axis/unidock/inputs/unidock_manifest.json
  results/trackA/method_axis/unidock/inputs/run_unidock_806.template.sh
All sha256 recorded for the amendment freeze chain.
"""
import csv
import hashlib
import json
import os
import subprocess

ROOT = r"D:\harness\projects\d1_gramneg_vs"
B = os.path.join(ROOT, "results", "trackA", "method_axis")
SUBSETS = os.path.join(B, "subsets_diffdock.csv")
TARGET_META = os.path.join(ROOT, "results", "trackA", "full", "target_metadata.csv")
REC_MANIFEST = os.path.join(ROOT, "data", "receptors", "receptor_manifest.csv")
DRUGCLIP_ENV = r"F:\Drugclip\Drugclip_run\env\python.exe"
PREPARE_MOLS = r"F:\Drugclip\Drugclip_run\prepare_mols.py"
TARGETS = ["LpxC_Ec", "GyrB_Ec", "Gyrase_Ec", "LpxC_Pa", "DHFR_Ec",
           "GyrA_Ec", "FabI_Ec", "TEM1_Ec", "LepB_Ec"]


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_csv(p):
    return list(csv.DictReader(open(p, encoding="utf-8-sig")))


def wsl_path(p):
    # D:\a\b -> /mnt/d/a/b
    p = os.path.abspath(p)
    drive = p[0].lower()
    return "/mnt/" + drive + p[2:].replace("\\", "/")


def main():
    subs = read_csv(SUBSETS)
    tmeta = {r["target"]: r for r in read_csv(TARGET_META)}
    recs = {(r["target_key"], r["role"]): r for r in read_csv(REC_MANIFEST)}
    rep = {"subsets_diffdock_sha256": sha256(SUBSETS), "targets": {}}

    # ---------------- DrugCLIP library ----------------
    ddir = os.path.join(B, "drugclip", "inputs")
    os.makedirs(ddir, exist_ok=True)
    lib_csv = os.path.join(ddir, "library_806.csv")
    with open(lib_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["mol_id", "target", "ligand_id", "label", "smiles"])
        for r in subs:
            w.writerow([f"{r['target']}__{r['ligand_id']}", r["target"], r["ligand_id"],
                        r["label"], r["smiles"]])
    lib_lmdb = os.path.join(ddir, "library_806.lmdb")
    print("[drugclip] building LMDB:", lib_lmdb)
    subprocess.run([DRUGCLIP_ENV, PREPARE_MOLS, lib_csv, lib_lmdb],
                   check=True, capture_output=True, text=True)
    rep["drugclip"] = {
        "library_csv": os.path.relpath(lib_csv, ROOT), "library_csv_sha256": sha256(lib_csv),
        "library_lmdb": os.path.relpath(lib_lmdb, ROOT),
        "library_lmdb_sha256": sha256(lib_lmdb),
        "n_molecules": sum(1 for _ in subs),
        "pocket_manifest": "results/trackA/method_axis/drugclip/pocket_manifest.json",
        "pocket_manifest_sha256": sha256(os.path.join(B, "drugclip", "pocket_manifest.json")),
        "env_python": DRUGCLIP_ENV,
    }

    # ---------------- Uni-Dock inputs ----------------
    udir = os.path.join(B, "unidock", "inputs")
    os.makedirs(udir, exist_ok=True)
    all_rows = []
    for t in TARGETS:
        rk = "GyrB_Ec" if t == "Gyrase_Ec" else t
        rec = recs.get((rk, "holo"))
        if not rec:
            raise RuntimeError(f"no holo receptor manifest row for {t}")
        box_size = float(tmeta[t]["box_size_A"])
        center = (float(rec["grid_center_x"]), float(rec["grid_center_y"]),
                  float(rec["grid_center_z"]))
        rec_pdbqt = os.path.join(ROOT, rec["pdbqt"])
        rows = []
        for r in subs:
            if r["target"] != t:
                continue
            lig_pdbqt = os.path.join(ROOT, r["ligand_pdbqt"])
            rows.append({
                "subset_id": r["subset_id"], "target": t, "ligand_id": r["ligand_id"],
                "label": r["label"], "ligand_pdbqt": os.path.relpath(lig_pdbqt, ROOT),
                "receptor_pdbqt": os.path.relpath(rec_pdbqt, ROOT),
                "center_x": center[0], "center_y": center[1], "center_z": center[2],
                "size_x": box_size, "size_y": box_size, "size_z": box_size,
            })
            all_rows.append(rows[-1])
        out = os.path.join(udir, f"{t}.csv")
        with open(out, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        n_missing = sum(1 for r in rows if not os.path.exists(os.path.join(ROOT, r["ligand_pdbqt"])))
        rep["targets"][t] = {
            "n_ligands": len(rows), "n_missing_ligand_pdbqt": n_missing,
            "input_csv": os.path.relpath(out, ROOT), "input_csv_sha256": sha256(out),
            "receptor_pdbqt": os.path.relpath(rec_pdbqt, ROOT),
            "receptor_pdbqt_sha256": sha256(rec_pdbqt),
            "grid_center": list(center), "box_size_A": box_size,
            "pdb_id": rec.get("pdb_id", ""), "grid_source": rk,
        }
    all_csv = os.path.join(udir, "unidock_inputs_all.csv")
    with open(all_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
        w.writeheader()
        w.writerows(all_rows)
    rep["unidock"] = {"unidock_inputs_all": os.path.relpath(all_csv, ROOT),
                      "unidock_inputs_all_sha256": sha256(all_csv),
                      "n_ligands": len(all_rows)}
    # WSL command template
    tmpl = ["#!/bin/bash",
            "# Amendment #14 Uni-Dock 806-subset template (Lead runs in WSL; paths = /mnt/d/...)",
            "set -euo pipefail",
            'UNIDOCK="${UNIDOCK_BIN:-unidock}"',
            'PY="python3"',
            'TRACKA="/mnt/d/harness/projects/d1_gramneg_vs/results/trackA"',
            'OUT="/mnt/d/harness/projects/d1_gramneg_vs/results/trackA/method_axis/unidock/scores"',
            'mkdir -p "$OUT"',
            "while IFS=, read -r subset_id target ligand_id label ligand_pdbqt receptor_pdbqt cx cy cz sx sy sz; do",
            '  [ "$ligand_id" = "ligand_id" ] && continue',
            '  "$UNIDOCK" --receptor "$TRACKA/../$receptor_pdbqt" --ligand "$TRACKA/../$ligand_pdbqt" \\',
            '    --center_x "$cx" --center_y "$cy" --center_z "$cz" \\',
            '    --size_x "$sx" --size_y "$sy" --size_z "$sz" \\',
            '    --scoring vina --exhaustiveness 4 --seed 42 --num_modes 9 \\',
            '    --out "$OUT/${target}__${ligand_id}_out.pdbqt" --log "$OUT/${target}__${ligand_id}_log.txt"',
            'done < /mnt/d/harness/projects/d1_gramneg_vs/results/trackA/method_axis/unidock/inputs/unidock_inputs_all.csv',
            ""]
    tmpl_path = os.path.join(udir, "run_unidock_806.template.sh")
    open(tmpl_path, "w", newline="\n", encoding="utf-8").write("\n".join(tmpl))
    rep["unidock"]["command_template"] = os.path.relpath(tmpl_path, ROOT)
    rep["unidock"]["command_template_sha256"] = sha256(tmpl_path)
    rep["unidock"]["note"] = ("paths in template are /mnt/d/... WSL paths; CSV rows use repo-relative "
                              "paths so prefix TRACKA root accordingly; adjust template before running")
    man = os.path.join(udir, "unidock_manifest.json")
    json.dump(rep, open(man, "w"), indent=1, ensure_ascii=False)
    print("wrote", man, "sha256:", sha256(man))
    print("drugclip library:", rep["drugclip"]["n_molecules"], "mols |",
          rep["drugclip"]["library_lmdb_sha256"][:16])
    print("unidock ligands:", rep["unidock"]["n_ligands"], "| inputs sha:",
          rep["unidock"]["unidock_inputs_all_sha256"][:16])


if __name__ == "__main__":
    main()
