#!/usr/bin/env python3

import pandas as pd

import argparse
import json
import os


def select_fraction(evals):
    return (
        evals.groupby("reduction")
        .apply(lambda df: df.sort_values(by="threshold").iloc[-1]["fracReads"])
        .sort_values()
        .iloc[-20:]
        .index.tolist()
    )


def select_error(evals):
    return (
        evals.groupby("reduction")
        .apply(lambda df: df.sort_values(by="threshold").iloc[-1]["cumErrorRate"])
        .sort_values()
        .iloc[:20]
        .index.tolist()[::-1]
    )


def select_percentage(evals):
    return (
        evals.groupby("reduction")["reduction"]
        .count()
        .sort_values()
        .iloc[-20:]
        .index.tolist()
    )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--csv", "-c", help="csv file with gathered mapeval results", required=True
    )
    parser.add_argument(
        "--top",
        "-t",
        help="how many MSRs to keep in each category",
        required=False,
        default=20,
    )
    parser.add_argument(
        "--dir",
        "-d",
        help="directory with SSR json definitions.\nIf specified, the program will copy selected MSRs in a subfolder",
        required=False,
    )

    args = parser.parse_args()

    evals = pd.read_csv(args.csv)

    hpc_err, hpc_frac = evals[
        (evals["simulator"] == "nanosim")
        & (evals["renamed"] == "HC")
        & (evals["threshold"] == 60)
    ][["cumErrorRate", "fracReads"]].iloc[0]

    subset = evals[
        (evals["simulator"] == "nanosim")
        & (evals["cumErrorRate"] <= hpc_err)
        & (evals["fracReads"] >= hpc_frac)
        & (evals["type"] == "RF")
    ]

    MSRs_e = select_error(subset)
    MSRs_f = select_fraction(subset)
    MSRs_p = select_percentage(subset)

    MSRs = set(MSRs_e + MSRs_p + MSRs_f)

    MSR_f = (
        evals[
            (evals["simulator"] == "nanosim")
            & (evals["reduction"].isin(MSRs_f))
            & (evals["threshold"] == 0)
        ]
        .sort_values(by="fracReads")
        .iloc[-1]["reduction"]
    )

    if args.dir is not None:

        dest_dir = os.path.abspath(os.path.join(args.dir, "MSRs"))
        os.makedirs(dest_dir, exist_ok=True)

        for msr in MSRs.union({"raw", "homopolymerCompression"}):
            src = os.path.abspath(os.path.join(args.dir, f"{msr}.json"))
            dest = os.path.join(dest_dir, f"{msr}.json")
            os.symlink(src, dest)

    else:
        with open("msrs.json", "w") as file:
            json.dump(
                dict(msrs=list(MSRs), msr_e=MSRs_e[-1], msr_f=MSR_f, msr_p=MSRs_p[-1]),
                file,
            )


if __name__ == "__main__":
    main()
