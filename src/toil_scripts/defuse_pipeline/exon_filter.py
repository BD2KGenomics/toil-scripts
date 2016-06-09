#!/usr/bin/env python2.7

import csv
import cPickle

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
                setattr(self, attr, value)
            except ValueError:
                pass


def get_exons(gtf_fname):
    pickle_name = 'exons_{}'.format(gtf_fname)
    exons = {}
    try:
        with open(pickle_name, 'rb') as f:
            exons = cPickle.load(f)
    except IOError:
        with open(gtf_fname, 'rb') as f, open(pickle_name, 'wb') as g:
            reader = csv.reader(f, delimiter='\t')
            line = reader.next()
            while line[0].startswith('##'):
                line = reader.next()
            while True:
                gtf = GTF(line)
                if gtf.feature == 'exon':
                    exons[gtf.gene_id] = gtf
                try:
                    line = reader.next()
                except StopIteration:
                    break
            cPickle.dump(exons, g)
    return exons




def main():
    gtf_fname = 'gencode.v19.annotation.gtf'

    exons = get_exons(gtf_fname)
    print exons

if __name__ == '__main__':
     main()
