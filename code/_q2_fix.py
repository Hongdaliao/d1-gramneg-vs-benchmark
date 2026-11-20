from pathlib import Path

p = Path(r"manuscript/discussion_draft.md")
t = p.read_text(encoding="utf-8")
t = t.replace("Discussion 中 MAES 统一写 descriptive only，不写 fitted/general-purpose。",
              "Discussion 中 MAES 统一写 descriptive only，不写概率/排名承诺。")
p.write_text(t, encoding="utf-8")

p = Path(r"manuscript/main_skeleton.md")
t = p.read_text(encoding="utf-8")
t = t.replace("do **not** describe the anchor as validation claim or use it as a statistical-power argument.",
              "do **not** describe the anchor as a validation claim or use it as a statistical-power argument.")
p.write_text(t, encoding="utf-8")

p = Path(r"manuscript/methods_anchor_draft.md")
t = p.read_text(encoding="utf-8")
t = t.replace(
    'It is **not validation claim**, is **not part of the main statistical-power argument**, and does **not** support statements such as "validated on LIT-PCBA" or "validation claim".',
    'It is **not a validation claim**, is **not part of the main statistical-power argument**, and does **not** support a validation claim ("validated on LIT-PCBA").',
)
p.write_text(t, encoding="utf-8")

p = Path(r"manuscript/results_trackA_draft.md")
t = p.read_text(encoding="utf-8")
t = t.replace("**not validation claim** and **not part of the main statistical-power argument**.",
              "**not a validation claim** and **not part of the main statistical-power argument**.")
t = t.replace("and is not counted as validation claim.", "and is not counted as a validation claim.")
p.write_text(t, encoding="utf-8")

print("grammar fixes applied")
