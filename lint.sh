#!/bin/bash

set -euxo pipefail
codespell ./notes/*.tex
codespell ./thesis/*.tex
codespell ./slides/*.tex
