#!/bin/bash
set -euo pipefail

# input:
COMPARISON="$1"
# output:
METHYLKIT_RTL_CSV="$2"

TEMP_JSON=$(mktemp)
echo '{
    "postVpre1": "2021-06,2021-08,2021-09-1pre,2021-09-2post",
    "postVpre2": "2021-06,2021-08,2021-09-1pre,2021-09-2post,2021-11,2021-12",
    "postVpre3": "2021-06,2021-08,2021-09-1pre,2021-09-2post,2021-11,2021-12,2022-03",
    "recVpre1": "2021-06,2021-08,2021-09-1pre,2021-11",
    "recVpre2": "2021-06,2021-08,2021-09-1pre,2021-12",
    "recVpre3": "2021-06,2021-08,2021-09-1pre,2022-03",
    "recVpre4": "2021-06,2021-08,2021-09-1pre,2021-12,2022-03",
    "recVpre5": "2021-06,2021-08,2021-09-1pre,2021-11,2021-12,2022-03",
    "postVrec1": "2021-09-2post,2021-11",
    "postVrec2": "2021-09-2post,2021-12",
    "postVrec3": "2021-09-2post,2022-03",
    "postVrec4": "2021-09-2post,2021-12,2022-03",
    "postVrec5": "2021-09-2post,2021-11,2021-12,2022-03"
}' > "$TEMP_JSON"
TARGETS=$(python -c "import json; h=open('$TEMP_JSON'); print(json.load(h)['$COMPARISON']); h.close()")
rm "$TEMP_JSON"

Rscript --vanilla ./r/diff-m6a.r ./r/samples.tsv ALL \
    "$TARGETS" "$METHYLKIT_RTL_CSV"
