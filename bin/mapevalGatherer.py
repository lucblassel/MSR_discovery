#!/usr/bin/env python3

import json
import re
import gzip

import pandas as pd
from tqdm import tqdm

def getMetadata(filename, gzipped=True):

    opener = gzip.open if gzipped else open
    with opener(filename,  "rb") as file:
        for line in file:
            reduction, simulator, mapper = line.decode("utf-8").strip().split(", ")
            break
    return dict(reduction=reduction, simulator=simulator, mapper=mapper)

def readFiles(filenames, nReads=None, gzipped=True):
    cols = ["mapType", "threshold", "nMapped", "nErrors", "cumErrorRate", "cumNum"]
    short = {"raw": "raw", "homopolymerCompression": "HC"}

    regex = re.compile(r"^uniqueReduction_(\d)out_(.*)_partition(\d+)_configuration(\d+)$")

    if nReads is not None:
        read_num = json.load(open(nReads, "r"))

    dfs = []
    for filename in tqdm(filenames):

        metadata = getMetadata(filename, gzipped)
        reduction = metadata["reduction"]
        match = regex.match(reduction)
        d = pd.read_csv(filename, sep="\t", header=None, skiprows=1)

        d.columns = cols
        for name, val in metadata.items():
            d[name] = val

        s = short.get(reduction, "RF")
        d["type"] = s
        if nReads is not None:
            d["fracReads"] = d["cumNum"] / read_num.get(metadata["simulator"])
        if match is not None:
            gs = match.groups()
            d["renamed"] = f"R({gs[0]},{gs[1]},{gs[2]},{gs[3]})"
        else:
            d["renamed"] = s
        dfs.append(d)

    return pd.concat(dfs)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--files", nargs="+", required=True)
    parser.add_argument("--nreads", required=True, type=str)
    parser.add_argument("--output", required=True)
    parser.add_argument("--unzipped", required=False, default=False, action="store_true")
    args = parser.parse_args()


    evals = readFiles(args.files, args.nreads, not args.unzipped)
    evals.to_csv(args.output)

