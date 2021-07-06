#!/usr/bin/env python3

import sys

import pyBigWig
import numpy


def main():
    bw_fname, out_fname, kmer = sys.argv[1:4]
    kmer = int(kmer)
    bw = pyBigWig.open(bw_fname)
    chrom_sizes = bw.chroms()
    chroms = list(chrom_sizes.keys())
    chroms.sort()
    output = open(out_fname, 'w')
    for chrom in chroms:
        get_unique_seqs(chrom, 0, chrom_sizes[chrom], kmer, bw, output)
    output.close()


def get_unique_seqs(chrom, start, end, kmer, bw, output):
    print("{}:{}-{}".format(chrom, start, end))
    if end <= start:
        return
    counts = numpy.array(bw.values(chrom, start, end), numpy.float32)
    counts[numpy.where(numpy.isnan)] = 0
    starts = numpy.where((counts[:-1] != 1) & (counts[1:] == 1))[0] + 1 + start
    ends = numpy.minimum(numpy.where(
        (counts[:-1] == 1) & (counts[1:] != 1))[0] + kmer + start, end)
    if counts[0] == 1:
        starts = numpy.r_[start, starts]
    if counts[-1] == 1:
        ends = numpy.r_[ends, end]
    new_starts = numpy.zeros(starts.shape[0], numpy.int32)
    new_ends = numpy.zeros(starts.shape[0], numpy.int32)
    new_starts[0] = starts[0]
    new_ends[0] = ends[0]
    pos = 0
    for i in range(1, starts.shape[0]):
        if starts[i] <= new_ends[pos]:
            new_ends[pos] = ends[i]
        else:
            pos += 1
            new_starts[pos] = starts[i]
            new_ends[pos] = ends[i]
    starts = new_starts[:(pos + 1)]
    ends = new_ends[:(pos + 1)]
    for j in range(starts.shape[0]):
        print("{}\t{}\t{}".format(chrom, starts[j], ends[j]), file=output)


if __name__ == "__main__":
    main()
