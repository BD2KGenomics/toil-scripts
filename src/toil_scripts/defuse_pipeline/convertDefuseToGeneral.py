#!/usr/bin/env python2.7

import sys
import csv

def main():
    tsvFile = sys.argv[1]
    outpath = sys.argv[2]
    with open(tsvFile, 'rb') as f, open(outpath, 'wb') as g:
        reader = csv.reader(f, delimiter='\t')
        writer = csv.writer(g, delimiter='\t')
        header = reader.next()
        header = {key: index for index, key in enumerate(header)}

        # Also needs tot_reads but we need to sum splitr_count and
        defuse_features = ['gene_name1', 'gene1', 'genomic_strand1', 'gene_chromosome1',
                           'gene_name2', 'gene2', 'genomic_strand2', 'gene_chromosome2',
                           'gene_start1', 'genomic_break_pos1', 'genomic_break_pos2', 'gene_end2',
                           'splitr_count']
        pagasus_features = ['#5p_symbol', '5p_ensembl','5p_strand', '5p_chr',
                            '3p_symbol', '3p_ensembl', '3p_strand', '3p_chr',
                            '5p_start', '5p_end', '3p_start', '3p_end',
                            'split_reads']

        writer.writerow(pagasus_features)
        for line in reader:
            general = [line[header[feature]] for feature in defuse_features]
            tot_reads = str(int(line[header['splitr_count']]) + int(line[header['span_count']]))
            general.append(tot_reads)
            writer.writerow(general)



if __name__ == '__main__':
    main()
