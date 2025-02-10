#!/bin/bash

set -euxo pipefail

latexmk main
latexmk -pdf -cd thesis/thesis.tex
latexmk -cd slides/slides.tex

gcloud storage cp chen-evanchen-phd-math-2025-thesis.pdf gs://web.evanchen.cc/textbooks/chen-evanchen-phd-math-2025-thesis.pdf --cache-control="private,max-age=0"
gcloud storage cp slides/slides.pdf gs://web.evanchen.cc/textbooks/thesis-slides.pdf --cache-control="private,max-age=0"
gcloud storage cp main.pdf gs://web.evanchen.cc/textbooks/phd-notebook.pdf --cache-control="private,max-age=0"
