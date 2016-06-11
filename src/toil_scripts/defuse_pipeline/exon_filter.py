#!/usr/bin/env python2.7

import re
import csv
import cPickle

from collections import defaultdict, Counter

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
    try:
        transcripts = get_pickled_data('transcripts_{}.pkl'.format(gtf_fname))
    except IOError:
        # exons = get_exons(gtf_fname)
        transcripts = defaultdict(list)
        for gtf in get_gtfs(gtf_fname):
            if gtf.feature == 'exon':
                transcripts[gtf.transcript_id].append(gtf)
        pickle_data('transcripts_{}.pkl'.format(gtf_fname), transcripts)
    return transcripts


def main():
    gtf_fname = 'gencode.v19.annotation.gtf'
    defuse_output = 'results.tsv'

    transcripts = get_transcript_data(gtf_fname)
    reader = csv.reader(open(defuse_output, 'rb'), delimiter='\t')
    header = reader.next()
    line = reader.next()
    defuse_header = {key: index for index, key in enumerate(header)}
    defuse_features = ['gene1', 'gene2', 'genomic_break_pos1', 'genomic_break_pos2']
    for feature in defuse_features:
        print feature, line[defuse_header[feature]]
    print transcripts[line[defuse_header['gene1']]]



if __name__ == '__main__':
     main()
