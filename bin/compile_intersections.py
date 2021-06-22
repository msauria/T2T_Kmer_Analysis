#!/usr/bin/env python3

import sys


def main():
    kmer, regions, arrays, out_fname = sys.argv[1:5]
    regions = regions.split(',')
    arrays = arrays.split(',')
    N = len(regions)
    data = []
    print(N)
    for i in range(N):
        data.append(open("intersections_{}/chm13v1_{}/compiled.txt".format(
            kmer, regions[i])).readline().rstrip('\n').split(","))
        print(len(data[-1]))
    for i in range(N - 1):
        for j in range(i + 1, N):
            data[j][i + 1] = data[i][j + 1]
    output = open(out_fname, 'w')
    temp = [""] + regions + arrays + ['centro', 'kmers']
    print(",".join(temp), file=output)
    for i, r1 in enumerate(regions):
        print(",".join(data[i]), file=output)
    output.close()


if __name__ == "__main__":
    main()
