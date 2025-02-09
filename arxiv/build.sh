#!/usr/bin/env bash
set -euxo pipefail

./clean.py
latexmk main

tar cfhvz upload.tar.gz -- main-*.pdf *.tex main.bbl derivative.sty
