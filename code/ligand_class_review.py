import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "data/trackB/pose_manifest_curated.csv"
OUT = ROOT / "data/trackB/ligand_class_review.csv"
OUT_SUMMARY = ROOT / "data/trackB/ligand_class_review_summary.json"

manifest = pd.read_csv(MANIFEST)
main = manifest[manifest["inhibitor_main"].astype(str).str.lower().isin(["true", "1"])].copy()

true_inhibitors = {"5QE", "5QF", "5QG", "MBX"}
antibiotic_substrates = {"MIY", "DM2", "FUA", "RBT", "ERY", "LFX", "DXT", "CLM", "NOV"}
covalent_beta_lactams = {"AZR", "0W0", "PFV"}
fragments = {"1K8", "1KE", "A1AN8", "A1AOE", "A1AOF", "A1AON", "XE9"}


def classify(ligand_id: str):
    if ligand_id in true_inhibitors:
        return "inhibitor", "high", "RCSB chemical component name plus established AcrB pyranopyridine inhibitor series (5QE/5QF/5QG/MBX)"
    if ligand_id in antibiotic_substrates:
        return "antibiotic-substrate", "high", "RCSB chemical component name; known AcrB/MexB substrate or antibiotic (independent class review)"
    if ligand_id in covalent_beta_lactams:
        return "covalent", "high", "RCSB chemical component name; beta-lactam PBP3 covalent inhibitor (ligand_class_scheme §7)"
    if ligand_id in fragments:
        return "fragment", "medium", "RCSB chemical component name; MW below the declared inhibitor window / fragment-like; no inhibitor evidence in manifest"
    if ligand_id == "95X":
        return "substrate", "low", "RCSB chemical component name; functional class unresolved in available data; not counted as inhibitor evidence"
    return "substrate", "low", "RCSB chemical component name; no explicit substrate/inhibitor annotation in manifest"


rows = []
for _, row in main.iterrows():
    ligand_id = str(row["ligand_id"])
    suggested, confidence, evidence = classify(ligand_id)
    known = (
        "known inhibitor (pyranopyridine AcrB series)" if suggested == "inhibitor"
        else "known antibiotic/substrate" if suggested == "antibiotic-substrate"
        else "known covalent beta-lactam" if suggested == "covalent"
        else "fragment/low-MW; no inhibitor evidence" if suggested == "fragment"
        else "unresolved; not counted as inhibitor"
    )
    rows.append({
        "pose_id": row["pose_id"],
        "target": row["target"],
        "pdb_id": row["pdb_id"],
        "ligand_id": ligand_id,
        "ligand_name": row["ligand_name"],
        "mw": row["mw"],
        "current_class": row["ligand_class"],
        "current_class_confidence": row["class_confidence"],
        "evidence_source": evidence,
        "known_substrate_or_inhibitor": known,
        "suggested_class": suggested,
        "confidence": confidence,
        "note": row.get("notes", ""),
    })

review = pd.DataFrame(rows)
review.to_csv(OUT, index=False)
summary = {
    "n_poses": int(len(review)),
    "n_current_inhibitor": int((review["current_class"].astype(str).str.lower() == "inhibitor").sum()),
    "suggested_class_counts": review["suggested_class"].value_counts().to_dict(),
    "confidence_counts": review["confidence"].value_counts().to_dict(),
    "by_target": review.groupby(["target", "suggested_class"]).size().reset_index(name="n").to_dict("records"),
    "by_ligand_id": review.groupby(["target", "ligand_id", "suggested_class", "confidence"]).size().reset_index(name="n").to_dict("records"),
    "notes": [
        "Current manifest labels all 40 poses as inhibitor; this review separates true pump inhibitors (4 pyranopyridines) from known antibiotic substrates, covalent beta-lactams, and fragment-like unresolved ligands.",
        "No docking was rerun; evidence is RCSB chemical-component names from the manifest and the known ligand classes.",
        "The 55 non-polymer AcrB poses are not all inhibitor-class; 34 inhibitor-class poses were selected for the main set.",
    ],
}
OUT_SUMMARY.write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))
print("wrote", OUT, OUT_SUMMARY)
