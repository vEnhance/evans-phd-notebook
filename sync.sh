#!/bin/bash

set -euxo pipefail

# latexmk main
latexmk -pdf -cd thesis/thesis.tex

gcloud storage cp thesis/thesis.pdf gs://web.evanchen.cc/upload/thesis-draft.pdf --cache-control="private,max-age=0"
