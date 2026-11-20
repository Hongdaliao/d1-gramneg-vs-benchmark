from pathlib import Path

p = Path("manuscript/v3/manuscript_v3.md")
t = p.read_text(encoding="utf-8")

# M1 exact automated-label sentence in Results 3.5 and Discussion 4
auto = 'The `ligand_class` field in the manifest is an automated annotation (all 40 rows labelled "inhibitor"); expert re-annotation against RCSB names (Table S18) gives 4 true pump inhibitors, 24 antibiotic substrates, 3 covalent β-lactams, 8 fragment-like and 1 unresolved/low-confidence ligand (95X; annotated as substrate in Table S18, not counted as an inhibitor).'
result_old = 'The ligand-class composition is 4 true inhibitors, 24 antibiotic substrates, 3 covalent β-lactams, 8 fragment-like ligands and 1 substrate-like unresolved (95X); AcrB contributes 55 non-polymer poses in total; 34 carry the automated label "inhibitor", but expert re-annotation identifies only four true pump inhibitors. Class-stratified descriptive medians'
result_new = auto + ' AcrB contributes 55 non-polymer poses in total; 34 carry the automated label "inhibitor". Class-stratified descriptive medians'
if result_old in t:
    t = t.replace(result_old, result_new)
else:
    print('result sentence pattern missing')

disc_old = 'Track B is not an inhibitor benchmark: among 40 poses, only four are true pump inhibitors, 24 are antibiotic substrates, three are covalent β-lactams, eight are fragment-like and one is unresolved; AcrB contributes 55 non-polymer poses in total; 34 carry the automated label inhibitor, but expert re-annotation identifies only four true pump inhibitors.'
disc_new = 'Track B is not an inhibitor benchmark. ' + auto + ' AcrB contributes 55 non-polymer poses in total; 34 carry the automated label "inhibitor".'
if disc_old in t:
    t = t.replace(disc_old, disc_new)
else:
    print('discussion sentence pattern missing')

# Normalize residual 95X wording and duplicate
t = t.replace('**1 substrate-like unresolved (95X)** (95X)', '**1 unresolved/low-confidence ligand (95X; annotated as substrate in Table S18, not counted as an inhibitor)**')
t = t.replace('1 substrate-like unresolved (95X) ligand', '1 unresolved/low-confidence ligand (95X; annotated as substrate in Table S18, not counted as an inhibitor)')
t = t.replace('1 substrate-like unresolved (95X)', '1 unresolved/low-confidence ligand (95X; annotated as substrate in Table S18, not counted as an inhibitor)')
t = t.replace('1 is unresolved', '1 is an unresolved/low-confidence ligand (95X; annotated as substrate in Table S18, not counted as an inhibitor)')
t = t.replace('one is unresolved', 'one is an unresolved/low-confidence ligand (95X; annotated as substrate in Table S18, not counted as an inhibitor)')

# m1 mode-1 wording
t = t.replace('4/8 receptors with mode-1 redock >3 Å (GyrA 13.54 Å, TEM1 11.26 Å, LepB 6.59 Å, GyrB 3.11 Å); a best-of-nine criterion gives 2/8.',
              '4/8 receptors with top-ranked (mode-1) redock RMSD >3 Å (GyrA 13.54 Å, TEM1 11.26 Å, LepB 6.59 Å, GyrB 3.11 Å); a best-of-nine criterion gives 2/8.')
t = t.replace('receptor-preparation quality (4/8 mode-1 redock >3 Å; best-of-nine 2/8)',
              'receptor-preparation quality (4/8 top-ranked (mode-1) redock RMSD >3 Å; best-of-nine 2/8)')

# m4 BSS annotation
t = t.replace('BSS −0.281 full / −0.290 complete-case',
              'BSS −0.281 full-sample (146 rows, missing RMSD as negative; complete-case −0.290)')

# m5 version
t = t.replace('Manuscript version: submission-v1 (2026-11-19)', 'Manuscript version: v3-submission (2026-11-19)')

# m6 H1/H3 exact guidance
t = t.replace('Definitions and reporting status of H1 and H3 are provided in SI S24.',
              'H1 (co-folding subset) and H3 (membrane MD) were not evaluated; definitions and reporting status are provided in SI S24.')

p.write_text(t, encoding="utf-8")
print('v3 E fixes step 2 applied')
