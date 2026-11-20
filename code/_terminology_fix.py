from pathlib import Path

p = Path(r"manuscript/v3/manuscript_v3.md")
t = p.read_text(encoding="utf-8")
repls = [
    # Abstract (keep <=250 words)
    ("shows ΔRMSD −3.493 Å and ΔPLIF −0.289",
     "shows improved ligand RMSD (−3.493 Å) but reduced interface-contact recovery (−0.289)"),
    ("and the two data regimes they occupy", "and the two data regimes"),
    ("We analyse two regimes with regime-appropriate endpoints.", "We analyse two regimes with appropriate endpoints."),
    ("All data and code will be released upon publication.", "Data and code will be released upon publication."),
    # L47
    ("Uni-Dock [19] was run in WSL (v1.2.0) and DrugCLIP [20]",
     "Uni-Dock [19] v1.2.0 was run under Windows Subsystem for Linux (WSL), and DrugCLIP [20]"),
    ("ROC-AUC and BEDROC20 are secondary metrics", "ROC-AUC and BEDROC (α = 20) are secondary metrics"),
    # L55
    ("predicted AFDB receptors", "predicted AlphaFold DB (AFDB) receptors"),
    ("H2 is negative/underpowered: the original acceptable-pose-rate endpoint (MAES >= 0.5) difference was 0; the endpoint was replaced by paired RMSD/PLIF after the Amendment #1 trigger, which is disclosed as a protocol deviation. H1 (co-folding subset) and H3 (membrane MD) were not evaluated; definitions and reporting status are provided in SI S24.",
     "The predicted-receptor (structure-source) hypothesis was negative/underpowered: the original acceptable-pose-rate endpoint (MAES >= 0.5) difference was 0; the endpoint was replaced by paired RMSD/PLIF after the pre-specified QC trigger, which is disclosed as a protocol deviation. The co-folding subset and membrane-MD analyses were not evaluated; their definitions and reporting status are given in S24."),
    # L59
    ("with symmetry-aware ligand RMSD, PLIF/contact recovery and PoseBusters",
     "with symmetry-aware ligand RMSD, protein–ligand interaction fingerprint (PLIF) recovery and PoseBusters"),
    # W12 headings/methods
    ("### 2.5 MAES and the pre-registered W12 calibration gate", "### 2.5 The pre-registered MAES calibration gate"),
    ("### 3.6 W12 calibration gate: NO_GO", "### 3.6 MAES calibration gate: NO_GO"),
    ("The W12 gate used 146 pose observations across four method/arm labels: vina 40, vinardo 40, vina_afdb 40 and boltz2 26",
     "The pre-registered MAES calibration gate used 146 pose observations across four method/arm labels: Vina (holo receptors) 40, Vinardo 40, Vina with predicted (AFDB) receptors 40 and Boltz-2 26"),
    # L73
    ("Amendment A2 used a 20:1", "The pre-specified protocol used a 20:1"),
    ("AC2 was evaluated on the nine sets with n_val_act >= 20;",
     "The anchor's acceptance criteria AC1–AC7 are listed in Table 4. AC2 was evaluated on the nine sets with n_val_act >= 20;"),
    # L77
    ("This is a post-hoc exploratory analysis (Amendment #12, #13 and #14; sha256 `9ce4342e45dee330243e2af64b331232b0ba381b7f5406ab41cfb30c1ca222a6`, `d93d8ee16bb9ff8eb7737e7e96707da1d377b3b89a40bd599dabf96a9607478b` and `eebedbf484bb751d651622fc2ef96b8564d785abd861032d586cd437ea3a4dff`; Addendum 12a and 14a).",
     "This is a post-hoc exploratory analysis (pre-specified amendments to the analysis plan; sha256s and addenda are listed in S24)."),
    # L85
    ("W12 wild-cluster AUROC uses 4,999 Rademacher draws and W12 rho uses 2,000 cluster bootstrap resamples",
     "the MAES calibration-gate wild-cluster AUROC uses 4,999 Rademacher draws and the MAES calibration-gate rho uses 2,000 cluster bootstrap resamples"),
    ("We report MDES/descriptive power", "We report minimum detectable effect sizes (MDES) and descriptive power"),
    ("the W12 manifest is frozen", "the MAES calibration-gate manifest is frozen"),
    # L123
    ("H2 is negative/underpowered: the original acceptable-pose-rate endpoint difference was 0, and the endpoint was replaced after Amendment #1; predicted receptors do not support a confirmatory transfer claim.",
     "The predicted-receptor hypothesis was negative/underpowered: the original acceptable-pose-rate endpoint difference was 0, and the endpoint was replaced after the pre-specified QC trigger; predicted receptors do not support a confirmatory transfer claim."),
    # L127
    ("(vina 40, vinardo 40, vina_afdb 40, boltz2 26; 19 pose-level evaluation positives; 40 systems)",
     "(Vina 40, Vinardo 40, Vina-AFDB 40, Boltz-2 26; 19 pose-level evaluation positives; 40 systems)"),
    ("(BSS −0.281 full-sample", "(Brier skill score (BSS; negative values indicate worse than the no-skill baseline) −0.281 full-sample"),
    ("(vina 0.572, vina_afdb 0.454; holo 0.572, predicted_afdb 0.454; PBP3",
     "(Vina 0.572 versus Vina-AFDB 0.454; holo 0.572 versus predicted-receptor 0.454; PBP3"),
    # L129
    ("**Table 2. W12 calibration gate", "**Table 2. MAES calibration gate"),
    # L143
    ("H2 is negative/underpowered; the original acceptable-rate endpoint difference was 0.",
     "The predicted-receptor hypothesis was negative/underpowered; the original acceptable-rate endpoint difference was 0."),
    # L153
    ("exploratory; H2 negative", "exploratory; predicted-receptor hypothesis negative"),
    # L177
    ("13 NA samples", "13 not-assessable (NA) samples"),
    # L193
    ("H2 is negative/underpowered and its endpoint was replaced after Amendment #1.",
     "The predicted-receptor hypothesis was negative/underpowered, and its endpoint was replaced after the pre-specified QC trigger."),
    # L197
    ("The W12 gate is NO_GO", "The MAES calibration gate is NO_GO"),
    # L205
    ("H2 is negative/underpowered and its endpoint was replaced after Amendment #1; H1 and H3 were not reported.",
     "The predicted-receptor hypothesis is negative/underpowered, and its endpoint was replaced after the pre-specified QC trigger; the co-folding subset and membrane-MD analyses were not reported."),
    ("The W12 MAES gate is NO-GO", "The MAES calibration gate is NO-GO"),
    ("4,999 W12 wild; 2,000 W12 rho", "4,999 MAES calibration-gate wild; 2,000 MAES calibration-gate rho"),
    ("per-target MDES/descriptive power", "per-target minimum detectable effect sizes (MDES) and descriptive power"),
    # L217
    ("The W12 calibration manifest and report", "The MAES calibration-gate manifest and report"),
    ("W12 metrics", "MAES calibration-gate metrics"),
    # L316
    ("**Figure 6. W12 calibration gate", "**Figure 6. MAES calibration gate"),
]
for old, new in repls:
    if old not in t:
        print("MISSING:", old[:160])
    else:
        t = t.replace(old, new, 1)
p.write_text(t, encoding="utf-8")
print("written", p.stat().st_size)
