#!/usr/bin/env python3

import sys

import numpy

def main():
    anno_fname, kmer, aname, results_fname = sys.argv[1:5]
    arrays, anames, chroms, colors = load_arrays(anno_fname)
    data, anames, regions = load_data(results_fname)
    triu = numpy.triu_indices(len(regions), 1)
    minval = numpy.inf
    maxval = -numpy.inf
    nonzero = numpy.where(data['counts'][triu] > 0)[0]
    minval = min(minval, numpy.amin(data['counts'][triu][nonzero]))
    maxval = max(maxval, numpy.amax(data['counts'][triu][nonzero]))

    # Write colors
    write_colors(colors, aname, kmer)

    # Write ideogram
    write_ideogram(arrays, chroms, aname, kmer)

    # Write links
    hist_set = write_links(data, arrays, regions, aname, kmer, minval, maxval)

    # Write histograms
    write_histograms(hist_set, data, regions, arrays, anames, aname, kmer)

def load_arrays(fname):
    data = []
    chroms = []
    colors = {}
    for line in open(fname):
        chrom, start, end, region, array, color = line.rstrip().split()[:6]
        if chrom.encode('utf8') not in chroms:
            chroms.append(chrom.encode('utf8'))
        colors[array] = color
        data.append((chrom, int(start), int(end), region,
                     region.split("(")[0], array, color))
    data = numpy.array(data, dtype=numpy.dtype([('chrom', 'S6'),
                                                ('start', numpy.int32),
                                                ('end', numpy.int32),
                                                ('fullregion', 'S20'),
                                                ('region', 'S20'),
                                                ('array', 'S20'),
                                                ('color', 'S20')]))
    arrays = numpy.unique(data['array'])
    return data, arrays, chroms, colors

def load_data(fname):
    data = []
    fs = open(fname)
    temp = fs.readlines()
    header = temp.pop(0).rstrip('\n').split(',')
    n = len(temp)
    m = len(header) - n - 3
    for i, line in enumerate(temp):
        line = line.rstrip('\n').split(',')
        data.append(tuple([line[0], line[0].split('_')[0], tuple(line[1:(n + 1)]),
                    tuple(line[(n + 1):-2]), line[-2], line[-1]]))
    regions = header[1:(n + 1)]
    anames = header[(n + 1):-2]
    data = numpy.array(data, dtype=numpy.dtype([('region', 'S20'),
                                                ('array', 'S20'),
                                                ('counts', numpy.float32, (n,)),
                                                ('acounts', numpy.float32, (len(anames),)),
                                                ('centro', numpy.float32),
                                                ('total', numpy.int32)]))
    data['acounts'] -= data['counts'][numpy.arange(data.shape[0]),
                                      numpy.arange(data.shape[0])].reshape(-1, 1)
    for i in range(data.shape[0]):
        index = anames.index(data['array'][i].decode('utf8'))
        data['acounts'][i, :index] -= data['acounts'][i, index]
        data['acounts'][i, (index + 1):] -= data['acounts'][i, index]
        data['centro'][i] -= numpy.sum(data['acounts'][i, :]) + data['counts'][i, i]
    return data, anames, regions

def write_colors(colors, aname, kmer):
    output = open("circos_etc/colors_{}_{}.conf".format(aname, kmer), 'w')
    print("<<include circos_etc/cenSatColors.conf>>", file=output)
    print("<colors>", file=output)
    for k, v in colors.items():
        print("{} = {}".format(k, v), file=output)
    print("</colors>", file=output)
    output.close()

def write_ideogram(arrays, chroms, aname, kmer):
    output = open("circos_data/karyotype_{}_{}.txt".format(
        aname, kmer), 'w')
    for chrom in chroms:
        dchrom = chrom.decode('utf8')
        where = numpy.where(arrays['chrom'] == chrom)[0]
        if where.shape[0] == 0:
            continue
        start = arrays['start'][where[0]]
        end = arrays['end'][where[-1]]
        print("chr - hs{} {} {} {} {}".format(
            dchrom.lstrip('chr'), dchrom.lstrip('chr'),
            start, end, dchrom), file=output)
    for chrom in chroms:
        dchrom = chrom.decode('utf8')
        where = numpy.where(arrays['chrom'] == chrom)[0]
        if where.shape[0] == 0:
            continue
        for i in range(where.shape[0]):
            region = arrays['region'][where[i]].decode('utf8')
            array = arrays['array'][where[i]].decode('utf8')
            start = arrays['start'][where[i]]
            end = arrays['end'][where[i]]
            if i > 0:
                end2 = arrays['end'][where[i - 1]]
                if end2 < start:
                    print("band hs{} {} {} {} {} other".format(
                        dchrom.lstrip('chr'), "other", "other",
                        end2, start), file=output)
            print("band hs{} {} {} {} {} {}".format(
                dchrom.lstrip('chr'), region, region, start,
                end, array), file=output)
    output.close()

def write_links(data, arrays, regions, aname, kmer, minval, maxval):
    w = 9.5
    w2 = w / 19
    span = maxval - minval
    if aname == 'ct':
        cfunc = get_ct_color
    elif aname == 'mismatch':
        cfunc = get_mismatch_color
    elif aname == 'hor':
        cfunc = get_hor_color
    elif aname == 'mon':
        cfunc = get_mon_color
    elif aname == 'hsat':
        cfunc = get_hsat_color
    elif aname == 'bsat':
        cfunc = get_bsat_color
    elif aname == 'censat':
        cfunc = get_censat_color
    else:
        cfunc = get_all_color
    triu = numpy.triu_indices(len(regions), 1)
    nonzero = numpy.where(data['counts'][triu] > 0)[0]
    links = numpy.zeros(nonzero.shape[0], dtype=numpy.dtype([
        ('c1', 'S20'), ('s1', numpy.int32),('e1', numpy.int32),
        ('c2', 'S20'), ('s2', numpy.int32),('e2', numpy.int32),
        ('a1', 'S20'), ('a2', 'S20'), ('color', 'S20'),
        ('count', numpy.int32)]))
    for i, j in enumerate(nonzero):
        reg1 = numpy.where(arrays['region'] ==
                           regions[triu[0][j]].encode('utf8'))[0][0]
        c1, s1, e1, fr1, r1, a1, _ = arrays[reg1]
        reg2 = numpy.where(arrays['region'] ==
                           regions[triu[1][j]].encode('utf8'))[0][0]
        c2, s2, e2, fr2, r2, a2, _ = arrays[reg2]
        m1 = (s1 + e1) // 2
        m2 = (s2 + e2) // 2
        t1 = triu[0][j]
        t2 = triu[1][j]
        count = data['counts'][t1, t2]
        dc1 = c1.decode('utf8').replace('chr','hs')
        dc2 = c2.decode('utf8').replace('chr','hs')
        dr1, dr2 = fr1.decode('utf8'), fr2.decode('utf8')
        da1, da2 = a1.decode('utf8'), a2.decode('utf8')
        s1 = int(round(m1 - (count - minval) * w - w2 * span))
        e1 = int(round(m1 + (count - minval) * w + w2 * span))
        s2 = int(round(m2 - (count - minval) * w - w2 * span))
        e2 = int(round(m2 + (count - minval) * w + w2 * span))
        col = cfunc(da1, da2, dr1, dr2)
        links[i] = tuple([c1, s1, e1, c2, s2, e2, a1, a2, col, count])
    if aname == 'ct':
        where = numpy.where((links['a1'] == b'ct') | (links['a2'] == b'ct'))[0]
        hist_set = set(['ct'])
    elif aname == 'mismatch':
        where = numpy.where(((links['a1'] != links['a2'])) &
                            ((links['a1'] != b'ct') & (links['a2'] != b'ct')))[0]
        hist_set = set(['hor', 'mon', 'dhor', 'hor', 'mon', 'hsat1', 'hsat2',
                        'hsat3', 'bsat', 'gsat', 'censat'])
    elif aname == 'hor':
        where = numpy.where((((links['a1'] == b'hor') & (links['a2'] == b'hor')) |
                            ((links['a1'] == b'hor') & (links['a2'] == b'dhor')) |
                            ((links['a1'] == b'dhor') & (links['a2'] == b'dhor')) |
                            ((links['a1'] == b'hor') & (links['a2'] == b'dhor'))) &
                            ((links['a1'] != b'ct') & (links['a2'] != b'ct')))[0]
        hist_set = set(['hor', 'dhor'])
    elif aname == 'mon':
        where = numpy.where(((links['a1'] == b'mon') | (links['a2'] == b'mon')) &
                            ((links['a1'] != b'ct') & (links['a2'] != b'ct')))[0]
        hist_set = set(['mon'])
    elif aname == 'hsat':
        where = numpy.where(((links['a1'] == b'hsat1') | (links['a1'] == b'hsat2') |
                             (links['a1'] == b'hsat3')) & ((links['a2'] == b'hsat1') |
                             (links['a2'] == b'hsat2') | (links['a2'] == b'hsat3')) &
                            ((links['a1'] != b'ct') & (links['a2'] != b'ct')))[0]
        hist_set = set(['hsat1', 'hsat2', 'hsat3'])
    elif aname == 'bsat':
        where = numpy.where(((links['a1'] == b'bsat') | (links['a2'] == b'bsat')) &
                            ((links['a1'] != b'ct') & (links['a2'] != b'ct')))[0]
        hist_set = set(['bsat'])
    elif aname == 'censat':
        where = numpy.where(((links['a1'] == b'gsat') | (links['a2'] == b'gsat') |
                             (links['a1'] == b'censat') | (links['a2'] == b'censat')) &
                            ((links['a1'] != b'ct') & (links['a2'] != b'ct')))[0]
        hist_set = set(['gsat', 'censat'])
    else:
        where = numpy.where((links['a1'] == links['a2']) & (links['a1'] != b'ct') &
                            (links['a2'] != b'ct'))[0]
        hist_set = set(['hor', 'mon', 'dhor', 'hor', 'mon', 'hsat1', 'hsat2',
                        'hsat3', 'bsat', 'gsat', 'censat'])
    links = links[where]
    links = links[numpy.argsort(links['count'])]
    output = open("circos_data/links_{}_{}.txt".format(aname, kmer), 'w')
    for i in range(links.shape[0]):
        c1, s1, e1, c2, s2, e2, a1, a2, col, _ = links[i]
        c1 = c1.decode('utf8').replace('chr', 'hs')
        c2 = c2.decode('utf8').replace('chr', 'hs')
        col = col.decode('utf8')
        print("{} {} {} {} {} {} color={},z={}".format(
            c1, s1, e1, c2, e2, s2, col, i), file=output)
    output.close()
    return hist_set

def get_ct_color(a1, a2, r1, r2):
    return 'ct'

def get_mismatch_color(a1, a2, r1, r2):
    return 'mismatch'

def get_hor_color(a1, a2, r1, r2):
    supra1 = set(['hor_10_1', 'hor_10_2', 'hor_10_3', 'hor_12_2', 'hor_16_2',
                  'hor_19_8', 'hor_1_5', 'hor_3_1', 'hor_3_2', 'hor_3_4',
                  'hor_5_2', 'hor_5_5', 'hor_6_1', 'hor_7_2', 'dhor_3_2',
                  'dhor_6_1'])
    supra2 = set(['hor_13_3', 'hor_14_3', 'hor_15_3', 'hor_18_1',
                  'hor_18_3', 'hor_18_5', 'hor_20_1', 'hor_20_2', 'hor_20_3',
                  'hor_21_3', 'hor_22_9', 'hor_2_2', 'hor_4_1', 'hor_4_2',
                  'hor_4_3', 'hor_8_2', 'hor_9_1'])
    supra3 = set(['hor_11_2', 'hor_17_1', 'hor_17_2', 'hor_1_3', 'hor_1_6', 'hor_X_1'])
    supra45 = set(['hor_15_1', 'hor_15_2', 'hor_19_2', 'hor_19_4', 'hor_20_4', 'hor_22_2',
                   'hor_22_6', 'hor_4_4', 'hor_5_6', 'hor_5_8', 'hor_7_1',
                   'dhor_12_8'])
    r1 = r1.split('(')[0]
    r2 = r2.split('(')[0]
    if r1 in supra1 and r2 in supra1:
        return 'supra1'
    if r1 in supra2 and r2 in supra2:
        return 'supra2'
    if r1 in supra3 and r2 in supra3:
        return 'supra3'
    if r1 in supra45 and r2 in supra45:
        return 'supra45'
    return 'mismatch'

def get_mon_color(a1, a2, r1, r2):
    if a1 == 'mon':
        return a2
    else:
        return a1

def get_hsat_color(a1, a2, r1, r2):
    hsats = set(['hsat1', 'hast2', 'hsat3'])
    if a1 in hsats and a2 in hsats:
        if a1 == a2:
            return a1
        else:
            return 'mismatch'
    elif a1 in hsats:
        return a2
    else:
        return a1

def get_bsat_color(a1, a2, r1, r2):
    if a1 == a2:
        if r1.count('LSAU') > 0:
            if r2.count('LSAU') > 0:
                return 'bsat2'
            else:
                return 'mismatch'
        if r2.count('LSAU') > 0:
            return 'mismatch'
        else:
            return 'bsat'
    elif a1 == 'bsat':
        return a2
    else:
        return a1

def get_censat_color(a1, a2, r1, r2):
    censats = set(['censat', 'gsat'])
    if a1 in censats and a2 in censats:
        if a1 == a2:
            return a1
        else:
            return 'mismatch'
    elif a1 in censats:
        return a2
    else:
        return a1

def get_all_color(a1, a2, r1, r2):
    return a1

def write_histograms(hist_set, data, regions, arrays, anames, aname, kmer):
    output = open("circos_data/histogram_{}_{}.txt".format(aname, kmer), 'w')
    w = 250000
    w2 = w * 1.07
    keep = []
    for i, a in enumerate(data['array']):
        if a.decode('utf8') in hist_set:
            keep.append(i)
    keep = numpy.array(keep)
    mids = numpy.zeros(keep.shape[0], dtype=numpy.int32)
    chroms = numpy.zeros(keep.shape[0], dtype='S20')
    for i, j in enumerate(keep):
        r = numpy.where(arrays['region'] == regions[j].encode('utf8'))[0][0]
        mids[i] = (arrays['start'][r] + arrays['end'][r]) // 2
        chroms[i] = arrays['chrom'][r]
    for chrom in numpy.unique(chroms):
        where = numpy.where(chroms == chrom)[0]
        dist = mids[where[1:]] - mids[where[:-1]]
        i = 0
        s = None
        while i < dist.shape[0]:
            while i < dist.shape[0]:
                if dist[i] < w2:
                    break
                i += 1
            if i < dist.shape[0]:
                s = i
                while i < dist.shape[0]:
                    if dist[i] >= w2:
                        break
                    i += 1
                e = i
                m = numpy.mean(mids[where[s]:(where[e] + 1)])
                mids[where[s]:(where[e] + 1)] = numpy.round(numpy.linspace(
                    m - (e - s) * w2 / 2, m + (e - s) * w2 / 2,
                    e - s + 1)).astype(numpy.int32)
    prev_chrom = ""
    for i, j in enumerate(keep):
        mid = mids[i]
        s = mid - w // 2
        e = mid + w // 2
        r = numpy.where(data['region'][j] == arrays['region'])[0][0]
        array = arrays['array'][r].decode('utf8')
        total = data['total'][j]
        unique = data['counts'][j, j]
        chrom = arrays['chrom'][r].decode('utf8').replace('chr', 'hs')
        if chrom != prev_chrom:
            count = numpy.zeros(data['acounts'].shape[1] + 3)
            prev_chrom = chrom
        index = anames.index(array)
        acounts = numpy.cumsum(numpy.r_[unique, data['acounts'][j, :], data['centro'][j]] / total)
        print("{} {} {} {} z={},fill_color={}".format(
              chrom, s, e, acounts[0], acounts.shape[0], array), file=output)
        for k in range(1, acounts.shape[0]):
            if (acounts[k - 1] == acounts[k]) and (count[acounts.shape[0] - k] >= 2):
                continue
            if k - 1 == index:
                col = 'white'
            elif k == acounts.shape[0] - 1:
                col = 'centro'
            else:
                col = anames[k - 1]
            print("{} {} {} {} z={},fill_color={}".format(
                  chrom, s, e, acounts[k], acounts.shape[0] - k, col), file=output)
            count[acounts.shape[0] - k] += 1
        print("{} {} {} {} z={},fill_color={}".format(
              chrom, s, e, 1, 0, 'black'), file=output)
    output.close()

if __name__ == "__main__":
    main()
