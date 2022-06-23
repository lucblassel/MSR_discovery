#!/usr/bin/env python3

import re


def formatFile(file, count, msr, mapper, simulator, organism):
    regex = re.compile(
        r"^uniqueReduction_(\d)out_(.*)_partition(\d+)_configuration(\d+)$"
    )
    short = {"raw": "raw", "homopolymerCompression": "HPC"}
    type_ = short.get(msr, "MSR")
    renamed = type_

    if (match := regex.match(msr)) is not None:
        gs = match.groups()
        renamed = f"R({gs[0]},{gs[1]},{gs[2]},{gs[3]})"

    with open(file, "r") as inFile:
        for line in inFile:
            fields = line.strip().split()
            frac = int(fields[-1]) / count
            print(*fields, frac, type_, renamed, organism, simulator, mapper, sep="\t")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--file", "-f", required=True)
    parser.add_argument("--count", "-c", type=int, required=True)
    parser.add_argument("--msr", "-m", required=True)
    parser.add_argument("--mapper", "-p", required=True)
    parser.add_argument("--simulator", "-s", required=True)
    parser.add_argument("--organism", "-o", required=True)
    args = parser.parse_args()

    formatFile(**args.__dict__)
