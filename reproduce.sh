#!/usr/bin/env bash
# Regenerate every result, figure, generated table and the manuscript PDF.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT/simulation/src/separation"
R=../../results/separation
for c in finite kantorovich costs audit phase designs boundary coverage curve firstorder rdcompare rdlimit; do
  echo "== run.py $c"; python run.py "$c"
done
python run.py designs --rate --ns 10000 100000 1000000 --reps 200
python run.py --out "$R/coverage_n1e6_r5000" coverage --ns 1000000 --reps-list 5000 --seed 20260920
python plot_figures.py
cd "$ROOT/paper/tex"
python make_separation_tables.py
if command -v latexmk >/dev/null; then latexmk -pdf main_jci.tex
elif command -v tectonic >/dev/null; then tectonic main_jci.tex
else echo "No LaTeX engine found; skipped PDF build."; fi
