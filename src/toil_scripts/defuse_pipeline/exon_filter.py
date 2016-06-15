#!/usr/bin/env python2.7

import re
import os
import sys
import csv
import cPickle

import matplotlib.pyplot as plt

from collections import defaultdict, namedtuple


class GTF(object):

    attributes = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

    def __init__(self, line):
        if isinstance(line, str):
            line = line.strip().split('\t')
        if len(line) != 9:
            msg = '\n'.join(self.attributes)
            raise ValueError('GTF file format has changed. Must have these attributes: {}'.format(msg))
        for attr, value in zip(self.attributes, line):
            setattr(self, attr, value)
        for feature in self.attribute.split(';'):
            try:
                attr, value = feature.split()
                setattr(self, attr, value.replace('"', ''))
            except ValueError:
                pass


def get_gtfs(gtf_fname):
    with open(gtf_fname, 'rb') as f:
        reader = csv.reader(f, delimiter='\t')
        line = reader.next()
        while line[0].startswith('##'):
            line = reader.next()
        while True:
            yield GTF(line)
            try:
                line = reader.next()
            except StopIteration:
                break


def get_pickled_data(fname):
    try:
        with open(fname, 'rb') as f:
            return cPickle.load(f)
    except IOError:
        raise IOError("Pickle file does not exist")


def pickle_data(fname, obj):
    with open(fname, 'wb') as f:
        cPickle.dump(obj, f)

def get_transcript_data(gtf_fname):
    gtf_stats = namedtuple('GTF_STATS', 'seqname source feature start end score strand frame')
    transcript_data = defaultdict(list)
    reader = csv.reader(open(gtf_fname, 'rb'), delimiter='\t')
    while True:
        line = reader.next()
        if line[0].startswith('##'):
            line = reader.next()
        else:
            break
    while True:
        if line[2] != 'exon':
            pass
        else:
            gtf_line = gtf_stats(*line[:8])
            transcript_match = re.search('transcript_id "(\w+)\.*(\d*)"', line[8])
            exon_match = re.search('exon_id "(\w+)\.*(\d*)"', line[8])
            exon_number_match = re.search('exon_number (\d+)', line[8])
            if transcript_match and exon_match:
                transcript_name = transcript_match.group(1)
                transcript_version = transcript_match.group(2)
                exon_id = exon_match.group(1)
                exon_number = exon_number_match.group(1)
                transcript_data[transcript_name].append((exon_id, gtf_line, exon_number))
                reverse = False if gtf_line.strand == '+' else True
                sorted(transcript_data[transcript_name], key=lambda x: int(x[2]), reverse=reverse)
            else:
                pass
        try:
            line = reader.next()
        except StopIteration:
            break
    return transcript_data


def get_exon_coverage(exon_coverage):
    reader = csv.reader(open(exon_coverage, 'rb'), delimiter='\t')
    exon_data = {}
    for line in reader:
        exon_match = re.search('exon_id "(\w+)\.*(\d*)"', line[8])
        if exon_match:
            exon_data[exon_match.group(1)] = float(line[9]) / float(line[11])
    return exon_data



def read_rsem_data(rsem_path):
    rsem_data = namedtuple('RSEM_DATA', 'gene_name gene_version transcript_name transcript_version fpkm')
    rsem_reader = csv.reader(open(rsem_path, 'rb'), delimiter='\t')
    header = rsem_reader.next()
    header = {key: index for index, key in enumerate(header)}
    data = defaultdict(list)
    for line in rsem_reader:
        gene_id = line[header['gene_id']]
        gene_name, gene_version = gene_id.split('.')
        transcript_id = line[header['transcript_id']]
        transcript_name, transcript_version = transcript_id.split('.')
        fpkm =  line[header['FPKM']]
        data[gene_name].append(rsem_data(gene_name, gene_version, transcript_name, transcript_version, fpkm))
    return data



def main():
    defuse_output = '/home/jacob/munge/defuse/data/results.tsv'
    gtf_path = '/home/jacob/munge/defuse/data/gencode.v19.annotation.gtf'
    rsem_path = '/home/jacob/munge/defuse/data/rsem.isoforms.results'
    cov_path = '/home/jacob/munge/defuse/data/SRR1657557rnaAligned.sortedByCoord.out.bam.cov'

    rsem_data = read_rsem_data(rsem_path)

    exon_data = get_exon_coverage(cov_path)

    transcripts = get_transcript_data(gtf_path)
    reader = csv.reader(open(defuse_output, 'rb'), delimiter='\t')
    header = reader.next()
    line = reader.next()
    defuse_header = {key: index for index, key in enumerate(header)}
    defuse_features = ['gene1', 'gene2', 'genomic_break_pos1', 'genomic_break_pos2']
    for feature in defuse_features:
        print feature, line[defuse_header[feature]]

    rsem_datum1 = rsem_data[line[defuse_header['gene1']]][0]
    rsem_datum2 = rsem_data[line[defuse_header['gene2']]][0]

    for exon, gtf, number in transcripts[rsem_datum1.transcript_name]:
        print exon_data[exon], gtf.start, number

    print '##############################', line[defuse_header['genomic_break_pos1']], \
        line[defuse_header['genomic_break_pos2']]

    for exon, gtf, number in transcripts[rsem_datum2.transcript_name]:
        print exon_data[exon], gtf.start, number



if __name__ == '__main__':
     main()
