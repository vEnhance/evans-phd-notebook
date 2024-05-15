#!/bin/bash

set -euxo pipefail

latexmk main
gsutil cp main.pdf gs://web.evanchen.cc/upload/evans-phd-notebook.pdf
gsutil setmeta -h 'Cache-Control:private, max-age=0, no-transform' gs://web.evanchen.cc/upload/evans-phd-notebook.pdf
