#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

latexmk -C -outdir=build/latex vampire-manual.tex
latexmk -xelatex -interaction=nonstopmode -halt-on-error -file-line-error -outdir=build/latex vampire-manual.tex
cp build/latex/vampire-manual.pdf vampire-manual.pdf
