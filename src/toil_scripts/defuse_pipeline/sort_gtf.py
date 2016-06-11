#!/usr/bin/env python2.7

import csv
import sys
import os

from collections import defaultdict
from operator import itemgetter

def main():

    filename = sys.argv[1]

    reader = csv.reader(open(filename, 'rb'), delimiter='\t')
    line = reader.next()
    header = []
    while line[0].startswith('##'):
        header.append(line)
        line = reader.next()

    chrom = 0
    start = 3

    data = defaultdict(list)
    while True:
        if len(line) == 9:
            data[line[chrom]].append((line, int(line[start])))
        else:
            raise ValueError("GTF file must have 9 features")
        try:
            line = reader.next()
        except StopIteration:
            break

    chrom_format = 'chr'

    chroms = ['{}{}'.format(chrom_format, num) for num in range(1, 23) + ['X', 'Y', 'M']]

    base, ext = os.path.splitext(filename)
    output = '{}.sorted{}'.format(base, ext)

    out = open(output, 'w')

    for head in header:
        out.write('\t'.join(head) + '\n')

    for chrom in chroms:
        gtfs = sorted(data[chrom], key=itemgetter(1))
        lines = ['\t'.join(gtf) for gtf, _ in gtfs]
        for line in lines:
            out.write(line + '\n')

if __name__ == '__main__':
    main()