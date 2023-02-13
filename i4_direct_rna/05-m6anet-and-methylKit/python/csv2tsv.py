#!/usr/bin/env python
import sys, pandas

pandas.read_csv(sys.argv[1]).to_csv(
    sys.argv[2], sep="\t", index=False, compression="gzip",
)
