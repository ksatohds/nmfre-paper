# Table 5 of the main text (topic model) with era labels matching the coding of the era indicators.
# NEW file, 2026-11-08 (work memo A183).  No recomputation and no change to the numbers: reads
# tab_topic_theta_v2.tex (written by emit_tab_examples_boundary_r2_20261102.R, left untouched) and
# replaces only the era labels.  The coding is in make_topic_data_20261031.R: Early = Year <= 1865,
# Industrial = 1865 < Year <= 1945, Postwar = Year > 1945; the addresses in the corpus fall in
# 1789-1865 (20), 1869-1945 (20) and 1949-2021 (19).
# Usage: python -X utf8 relabel_tab_topic_theta_20261108.py   Output: tab_topic_theta_v3_20261108.tex
import io, os
here = os.path.dirname(os.path.abspath(__file__))
src = os.path.join(here, "tab_topic_theta_v2.tex")
dst = os.path.join(here, "tab_topic_theta_v3_20261108.tex")
t = io.open(src, encoding="utf-8", newline="").read()
reps = [("Early (1775--1865)     ", "Early (1789--1865)     "),
        ("Industrial (1865--1945)", "Industrial (1869--1945)"),
        ("Postwar (1945--)       ", "Postwar (1949--2021)   ")]
for a, b in reps:
    assert t.count(a) == 3, (a, t.count(a))
    t = t.replace(a, b)
io.open(dst, "w", encoding="utf-8", newline="").write(t)
# check: apart from the labels, every cell is unchanged
old = io.open(src, encoding="utf-8").read().splitlines()
new = t.splitlines()
same = all(o.split("&", 1)[1] == n.split("&", 1)[1] for o, n in zip(old, new)) and len(old) == len(new)
print("rows", len(new), "| numbers unchanged:", same)
print(t)
