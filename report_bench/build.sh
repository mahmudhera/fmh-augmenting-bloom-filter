#!/usr/bin/env bash
# Two-pass groff build of the report (pass 1 collects headings for the TOC).
set -euo pipefail
cd "$(dirname "$0")"
GROFF="groff -s -t -e -k -Tps -dpaper=letter -P-pletter -M . -m fmhdoc"
make_toc() {
  python3 - "$1" > gen/toc.ms <<'PY'
import sys
lines = [l.rstrip("\n").split("|") for l in open(sys.argv[1]) if l.startswith("TOC|")]
print(".nf\n.ta 0.45i 6.5iR\n.lc .\n.ps 10\n.vs 11.9p")
for _, lvl, num, title, page in lines:
    if lvl == "1":
        print(".sp 0.25v")
        print(f"\\f[B]{num}\t{title}\\f[]\x01{page}")
    else:
        print(f"\\h'0.45i'{num}\\h'|1.05i'{title}\x01{page}")
print(".lc\n.ps \\n(PS\n.vs \\n(VSp\n.fi")
PY
}
: > gen/toc.ms
$GROFF main.ms > main.ps 2> build1.log
make_toc build1.log
$GROFF main.ms > main.ps 2> build.log
grep -v '^TOC|' build.log || true
ps2pdf -dPDFSETTINGS=/prepress main.ps ../aif_benchmark_report.pdf
rm -f main.ps build1.log
echo "built ../aif_benchmark_report.pdf"
