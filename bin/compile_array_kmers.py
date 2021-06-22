#!/usr/bin/env python3

import sys

def main():
    prefix, regions, groups = sys.argv[1:4]
    regions = regions.split(',')
    groups = groups.split(',')
    kmers = {}
    for r in regions:
        for g in groups:
            for line in open("{}{}/{}.txt".format(prefix, r, g)):
                kmer, count = line.rstrip().split()[:2]
                kmers.setdefault(kmer, (g, []))
                kmers[kmer][1].append((r, count))
    print("Kmer\tType\tRegions\tCounts")
    keys = list(kmers.keys())
    keys.sort()
    for k in keys:
        r = []
        c = []
        for region, count in kmers[k][1]:
            r.append(region)
            c.append(count)
        print("{}\t{}\t{}\t{}".format(k, kmers[k][0], ",".join(r), ",".join(c)))

if __name__ == "__main__":
    main()
