#!/usr/bin/env python3
import os
import json

def partitionK(collection, minimum, k):
    """
    Generates all partitions in collection between sizes minimum and k
    """
    if len(collection) == 1:
        yield [collection]
        return

    first = collection[0]
    for smaller in partitionK(collection[1:], minimum - 1, k):
        if len(smaller) > k:
            continue
        if len(smaller) >= minimum:
            for n, subset in enumerate(smaller):
                yield smaller[:n] + [[first] + subset] + smaller[n + 1 :]
        if len(smaller) < k:
            yield [[first]] + smaller


def kPartitions(collection, k):
    """
    Generate all k-partitions of a given collection
    """
    yield from partitionK(collection, k, k)


def getPerms(k):
    """
    Generate output configuration for k nucleotide outputs
    """
    return {
        1: [["A"]],
        2: [["A", "T"], ["A", "C"]],
        3: [["A", "T", "C"], ["A", "C", "T"], ["C", "A", "T"]],
        4: [["A", "T", "C", "G"], ["A", "C", "T", "G"], ["A", "C", "G", "T"]],
    }.get(k, [])


def getOutputs(k, deletions=False):
    """
    Generate output configuration for k outputs (with or without deletion)
    """
    if k > 4 and not deletions:
       return []
    if k > 5:
        raise ValueError("k must be <= 5")

    perms = getPerms(k - 1) if deletions else getPerms(k)

    if deletions:
        newPerms = []
        for perm in perms:
            for i in range(len(perm) + 1):
                m = [x for x in perm]
                m.insert(i, ".")
                newPerms.append(m)
        return newPerms

    return perms


def reverseComplement(seq):
    """
    Returns reverse complement in nucleotides of input string sequence
    """
    complements = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(complements.get(x, x) for x in seq[::-1])


def completeFunction(partialFunction):
    """
    Return reverse complement robust reduction function from unique input mapping
    """

    if len(partialFunction) != 6:
        raise ValueError(f"Partial mapping must be of length 6 (got {len(partialFunction)})")

    finalFunction = dict()
    for k, v in partialFunction.items():
        finalFunction[k] = v
        finalFunction[reverseComplement(k)] = reverseComplement(v)
    for diNuc in ["AT", "CG", "GC", "TA"]:
        finalFunction[diNuc] = "."
    return finalFunction


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--output", "-o", help="output directory")
    args = parser.parse_args()

    if not os.path.isdir(args.output):
        os.makedirs(args.output)

    inputs = ["AA", "AC", "AG", "CA", "CC", "GA"]

    funcCounter = 0

    for nOut in range(2, 6):
        for i, mapping in enumerate(kPartitions(inputs, nOut)):
            for deletion in [False, True]:
                if nOut == 5 and not deletion:
                    continue
                configurations = getOutputs(nOut, deletion)
                for j, config in enumerate(configurations):
                    func = completeFunction({k: v for ks, v in zip(mapping, config) for k in ks})
                    json.dump(
                        func,
                        open(
                            os.path.join(
                                args.output,
                                f"uniqueReduction_{nOut}out_{deletion}_partition{i}_configuration{j}.json",
                            ),
                            "w",
                        ),
                    )
                    funcCounter += 1

    print(f"Saved {funcCounter} functions to disk")
