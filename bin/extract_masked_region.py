#!/usr/bin/env python3

import sys

import pysam
import pyBigWig


Non_TEs = set(['Beta', 'Satellite', 'Simple_repeat', 'Low_complexity', "Unknown"]) 

def main():
    fa_fname, bb_fname, region = sys.argv[1:4]
    fa = pysam.FastaFile(fa_fname)
    bb = pyBigWig.open(bb_fname)
    chrom, coords = region.split(':')
    start, end = coords.split('-')
    start, end = int(start), int(end)
    seq = fa.fetch(region=region)
    elements = bb.entries(chrom, start, end)
    if elements is None:
        elements = []
    for i in range(len(elements))[::-1]:
        if elements[i][2].split('\t')[7] in Non_TEs:
            elements.pop(i)
        else:
            elements[i] = (max(start, elements[i][0]), min(end, elements[i][1]), elements[i][2])
    if len(elements) > 0:
        masked = []
        if elements[0][0] > start:
            masked.append(seq[:(elements[0][0] - start)])
        for i in range(len(elements) - 1):
            masked.append("N" * (elements[i][1] - elements[i][0]))
            masked.append(seq[(elements[i][1] - start):(elements[i + 1][0] - start)])
        masked.append("N" * (elements[-1][1] - elements[-1][0]))
        if elements[-1][1] < end:
            masked.append(seq[(elements[-1][1] - start):end])
        seq = "".join(masked)
    print(">{}".format(region), file=sys.stdout)
    i = 0
    while i < len(seq):
        print(seq[i:min(i + 80, len(seq))], )
        i += 80


if __name__ == "__main__":
    main()
