from pathlib import Path
p=Path(r"code/maes/build_jcim_style_draft.py"); t=p.read_text(encoding="utf-8")
repls=[
("The panel also includes BamA and LptD entries contributed no inhibitor-class pose that passed curation", "The panel also includes BamA and LptD entries contributed no inhibitor-class pose that passed curation"),
("BamA and LptD entries contributed no inhibitor-class pose that passed curation, and lipid or substrate analogues were used only for sensitivity analyses (68 clean poses in total).",
 "BamA and LptD entries contributed no inhibitor-class pose that passed curation, and lipid or substrate analogues were used only for sensitivity analyses (68 clean poses in total). The co-crystal panel is therefore intentionally heterogeneous rather than a matched inhibitor set."),
("The seven primary targets run from AUC 0.300 (LpxC_Ec) to 0.852 (LpxC_Pa);",
 "Table 1 reports all nine targets, including the two auxiliary targets. The seven primary targets run from AUC 0.300 (LpxC_Ec) to 0.852 (LpxC_Pa);"),
("This work addresses a gap between how SBVS is used in Gram-negative antibacterial discovery and how it is benchmarked.",
 "This work addresses a gap between how SBVS is used in Gram-negative antibacterial discovery and how it is benchmarked. It separates benchmark-construction artefacts from target biology because the same panel supports both a full enrichment table and a balanced paired subset."),
("The design was organised around four questions. First,", "The questions were framed before the method-class results were inspected, and the answers are reported as descriptive contrasts rather than a ranking. First,"),
("The nine-target panel is broad for Gram-negative SBVS but still small,", "The unbalanced target sizes also mean that pooled means weight each target equally, not each active. The nine-target panel is broad for Gram-negative SBVS but still small,"),
("Both choices affect absolute enrichment.", "Matched decoys reduce property mismatch but cannot remove chemotype bias, and both choices affect absolute enrichment."),
("Practical utility is bounded by cost and data access.", "Practical utility is bounded by cost and data access. The practical cost is asymmetric: CPU-only Vina remains inexpensive, while GPU arms require scheduling, model-weight provenance and deterministic environment disclosure."),
("These three modes are reported together because they change conclusions in different ways.", "These three modes are reported together because they change conclusions in different ways."),
("Three benchmark/tool failure modes recur.", "These three modes are reported together because they change conclusions in different ways. Three benchmark/tool failure modes recur."),
("Nine targets were retained (seven primary and two auxiliary);", "Targets with fewer than 50 eligible actives were retained but flagged descriptive in the method axis. Nine targets were retained (seven primary and two auxiliary);"),
]
for old,new in repls:
 if old not in t: print('MISSING',old[:80])
 else: t=t.replace(old,new,1)
p.write_text(t,encoding='utf-8'); print('expanded')
