#!/bin/bash

set -euxo pipefail

latexmk main
latexmk -cd thesis/thesis.tex
PYTHONWARNINGS="ignore" gsutil cp main.pdf gs://web.evanchen.cc/upload/evans-phd-notebook.pdf
PYTHONWARNINGS="ignore" gsutil cp thesis/thesis.pdf gs://web.evanchen.cc/upload/thesis-draft.pdf
PYTHONWARNINGS="ignore" gsutil setmeta -h 'Cache-Control:private, max-age=0, no-transform' gs://web.evanchen.cc/upload/evans-phd-notebook.pdf
PYTHONWARNINGS="ignore" gsutil setmeta -h 'Cache-Control:private, max-age=0, no-transform' gs://web.evanchen.cc/upload/thesis-draft.pdf
