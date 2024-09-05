#!/bin/bash

set -euxo pipefail

latexmk main
latexmk -cd thesis/thesis.tex
gsutil cp main.pdf gs://web.evanchen.cc/upload/evans-phd-notebook.pdf
gsutil cp thesis/thesis.pdf gs://web.evanchen.cc/upload/thesis-draft.pdf
gsutil setmeta -h 'Cache-Control:private, max-age=0, no-transform' gs://web.evanchen.cc/upload/evans-phd-notebook.pdf
gsutil setmeta -h 'Cache-Control:private, max-age=0, no-transform' gs://web.evanchen.cc/upload/thesis-draft.pdf
