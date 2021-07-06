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
            kmer, regions[i])).readline().rstrip('\n').split(",")[1:])
        for j in range(len(data[-1])):
            if len(data[i][j]) > 0:
                data[i][j] = int(data[i][j])
    for i in range(N):
        for j in range(i, N):
            data[j][i] = data[i][j]
    for i in range(N):
        u = data[i][-2]
        for j in range(N, len(data[i]) - 2):
            data[i][j] -= u
        index = arrays.index(regions[i].split("_")[0]) + N
        for j in range(N, index):
            data[i][j] -= data[i][index]
        for j in range(index + 1, len(data[i]) - 3):
            data[i][j] -= data[i][index]
        for j in range(N, len(data[i]) - 3):
            data[i][-3] -= data[i][j]
    output = open(out_fname, 'w')
    temp = [""] + regions + arrays + ['centromere', 'unique', 'total']
    print(",".join(temp), file=output)
    for i, r1 in enumerate(regions):
        print("{},{}".format(r1, ",".join([str(x) for x in data[i]])), file=output)
    output.close()


if __name__ == "__main__":
    main()
