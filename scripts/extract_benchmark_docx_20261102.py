# Plain text of the four benchmark reports (Rmd rendered to .docx; the Ubuntu re-run of 2026-03-01) so that the numbers
# quoted in Appendix G can be checked and tagged with check_src_tags.py (work memo A95).  NEW file, 2026-11-02.
# Reads the .docx (zip) and writes one .txt per report: paragraphs and table rows, cells separated by " | ".
# No computation.  Usage: python -X utf8 -I extract_benchmark_docx_20261102.py   (run in scripts_results/)
# Output: benchmark_docx_text_20261102/<name>.txt and extract_benchmark_docx_20261102.log
import hashlib, os, sys, zipfile
import xml.etree.ElementTree as ET
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
NAMES = ["LMM線形混合モデル20260112", "MIMIC潜在変数回帰20260111", "PPR射影追跡回帰20260111", "STMトピックモデル20260112"]
W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
OUT = os.path.join(HERE, "benchmark_docx_text_20261102")
os.makedirs(OUT, exist_ok=True)
def para_text(p):
    return "".join(t.text or "" for t in p.iter(W + "t"))
log = []
for name in NAMES:
    path = os.path.join(SRC, name + ".docx")
    raw = open(path, "rb").read()
    md5 = hashlib.md5(raw).hexdigest()
    with zipfile.ZipFile(path) as z:
        root = ET.fromstring(z.read("word/document.xml"))
    body = root.find(W + "body")
    lines = []
    for el in body:
        if el.tag == W + "p":
            s = para_text(el)
            if s.strip():
                lines.append(s)
        elif el.tag == W + "tbl":
            for tr in el.iter(W + "tr"):
                cells = [" ".join(para_text(p) for p in tc.iter(W + "p")).strip() for tc in tr.iter(W + "tc")]
                lines.append(" | ".join(cells))
    out = os.path.join(OUT, name + ".txt")
    open(out, "w", encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")
    log.append("%s.docx md5 %s -> %s (%d lines)" % (name, md5, os.path.relpath(out, HERE), len(lines)))
open(os.path.join(HERE, "extract_benchmark_docx_20261102.log"), "w", encoding="utf-8").write("\n".join(log) + "\nDONE\n")
print("\n".join(log))
