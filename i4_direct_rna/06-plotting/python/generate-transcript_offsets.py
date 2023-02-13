#!/usr/bin/env python3
import pandas as pd
import re


def gtf_attr_ENSG(a):
    m = re.search(r'ENSG[^"]+', a)
    return m.group() if m else float("nan")


def gtf_attr_ENST(a):
    m = re.search(r'ENST[^"]+', a)
    return m.group() if m else float("nan")


gtf = pd.read_csv(
    "gencode.v41.annotation.gtf.gz", sep="\t", comment="#", names=[
        "name", "source", "feature", "start", "end",
        "score", "strand", "frame", "attribute",
    ],
)
gtf["ENSG"] = gtf["attribute"].apply(gtf_attr_ENSG)
gtf["ENST"] = gtf["attribute"].apply(gtf_attr_ENST)


G = gtf.loc[
    gtf["feature"].isin({"gene", "transcript"}),
    ["ENSG", "ENST", "feature", "name", "start", "end", "strand"],
].copy()
assert set(G[G["feature"]=="transcript"][["ENST", "start"]].groupby("ENST", as_index=False).aggregate(
    lambda vv: len(list(vv))
)["start"].value_counts().index) == {1}
assert len(G[["ENST", "name", "start"]].dropna()) == len(G[["ENST", "name", "start"]].dropna().drop_duplicates())


GT = pd.merge(
    G.loc[G["ENST"].isnull(), ["ENSG", "name", "start"]].rename(columns={"start": "ENSG_start"}),
    G[["ENSG", "ENST", "name", "start"]].dropna().rename(columns={"start": "ENST_start"}),
)
assert len(G[["ENSG", "ENST", "name", "start"]].dropna()) == len(GT)


transcript_offsets = GT.copy()
transcript_offsets["offset"] = transcript_offsets["ENST_start"] - transcript_offsets["ENSG_start"]
transcript_offsets = transcript_offsets[["ENST", "ENSG", "offset"]].copy()
transcript_offsets.to_csv("transcript_offsets.tsv.xz", sep="\t", index=False)
