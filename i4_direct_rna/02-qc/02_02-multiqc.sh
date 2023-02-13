#!/bin/bash
set -euo pipefail

# input:
DIR_WITH_ALL_PYCOQCS="$1"

multiqc --pdf "$DIR_WITH_ALL_PYCOQCS"
