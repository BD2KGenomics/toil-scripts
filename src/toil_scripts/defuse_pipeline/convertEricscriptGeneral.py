#!/usr/bin/env python2.7

import sys
import csv
import re
import cPickle

class GTF(object):
    def __init__(self, line):
        """
        Converts a GTF record into an python object

        :param line str|list: GTF record
        """
        if isinstance(line, str):
            line = line.strip().split('\t')
        if len(line) != 9:
            msg = '\n'.join(self.attributes)
            raise ValueError('GTF file format has changed. Must have these attributes: {}'.format(msg))
        for (attr, func), value in zip(self.attributes, line):
            setattr(self, attr, func(value))
        for feature in self.attribute.split(';'):
            try:
                attr, value = feature.split()
                setattr(self, attr, value.replace('"', ''))
            except ValueError:
                pass

    attributes = [('seqname', str), ('source', str), ('feature', str), ('start', int), ('end', int),
                  ('score', str), ('strand', str), ('frame', str), ('attribute', str)]

    def __repr__(self):
        if self.feature == 'gene':
            return 'GTF({}, gene:{}, start:{}, length:{})'.format(self.seqname, self.gene_name,
                                                                  self.start, self.end-self.start+1)
        elif self.feature == 'start_codon':
            return 'GTF({}, start codon:{}, start:{}, length:{})'.format(self.seqname,
                                                                         self.gene_name,
                                                                         self.start,
                                                                         self.end - self.start + 1)
        elif self.feature == 'transcript':
            return 'GTF({}, transcript:{}, start:{}, length:{})'.format(self.seqname,
                                                                        self.transcript_name,
                                                                        self.start,
                                                                        self.end-self.start+1)
        elif self.feature == 'exon':
            return 'GTF({}, exon:{}, start:{}, length:{}, id:{})'.format(self.seqname,
                                                                         self.exon_number,
                                                                         self.start,
                                                                         self.end-self.start+1,
                                                                         self.exon_id)
        elif self.feature == 'CDS':
            return 'GTF({}, CDS:{}, start:{}, length:{}, id:{})'.format(self.seqname,
                                                                         self.exon_number,
                                                                         self.start,
                                                                         self.end-self.start+1,
                                                                         self.exon_id)
        elif self.feature == 'start_codon':
            return 'GTF({}, start codon:{}, start:{}, length:{})'.format(self.seqname,
                                                                         self.gene_name,
                                                                         self.start,
                                                                         self.end - self.start + 1)

    def __hash__(self):
        if self.feature == 'gene':
            return hash(self.gene_id)
        elif self.feature == 'transcript':
            return hash(self.transcript_id)
        elif self.feature == 'CDS':
            return hash(self.exon_id)
        elif self.feature == 'exon':
            return hash(self.exon_id)

    def __eq__(self, other):
        try:
            if self.feature == 'gene':
                return self.gene_id == other.gene_id
            elif self.feature == 'transcript':
                return self.transcript_id == other.transcript_id
            elif self.feature == 'CDS':
                return self.exon_id == other.exon_id
            elif self.feature == 'exon':
                return self.exon_id == other.exon_id
        except AttributeError:
            return False


def get_gene_gtf(infile):
    """
    Creates a dictionary for identifying coding sequences

    :param infile:
    :return:
    """
    gene_model = {}
    for line in infile:
        if line.startswith('##'):
            continue
        gtf = GTF(line)
        if gtf.feature == 'gene':
            gene_model[gtf.gene_name] = gtf
            gene_model[gtf.gene_id] = gtf
    return gene_model



def get_pickled_data(fname):
    try:
        with open(fname, 'rb') as f:
            return cPickle.load(f)
    except IOError:
        raise IOError("Pickle file does not exist")


def pickle_data(fname, obj):
    with open(fname, 'wb') as f:
        cPickle.dump(obj, f)


def main():
    tsvFile = sys.argv[1]
    outpath = sys.argv[2]
    with open(tsvFile, 'rb') as f, open(outpath, 'wb') as g:
        reader = csv.reader(f, delimiter='\t')
        writer = csv.writer(g, delimiter='\t')
        header = reader.next()
        header = {key: index for index, key in enumerate(header)}

        gtf_path = '/home/jacob/munge/reference/gencode.v19.annotation.gtf'
        gene_model_pickle = '/home/jacob/munge/reference/gencode.v19.annotation.gtf.gene.pkl'

        try:
            gene_model = get_pickled_data(gene_model_pickle)
        except IOError:
            gene_model = get_gene_gtf(open(gtf_path, 'r'))
            pickle_data(gene_model_pickle, gene_model)

        # Also needs tot_reads but we need to sum splitr_count and
        mapped_features = ['GeneName1', 'EnsemblGene1', 'strand1', 'chr1',
                           'GeneName2', 'EnsemblGene2', 'strand2', 'chr2',
                           'gene_start1', 'Breakpoint1', 'Breakpoint2', 'gene_end2',
                           'spanningreads', 'tot_reads']

        pagasus_features = ['#5p_symbol', '5p_ensembl','5p_strand', '5p_chr',
                            '3p_symbol', '3p_ensembl', '3p_strand', '3p_chr',
                            '5p_start', '5p_end', '3p_start', '3p_end',
                            'split_reads', 'tot_reads']

        general_format = dict(zip(pagasus_features, mapped_features))

        print(general_format)

        eric_features = ['GeneName1', 'EnsemblGene1', 'strand1', 'chr1',
                         'GeneName2', 'EnsemblGene2', 'strand2', 'chr2',
                         'Breakpoint1', 'Breakpoint2', 'spanningreads',
                         'crossingreads']

        writer.writerow(pagasus_features)
        general = []
        for line in reader:
            if 'Unable to predict breakpoint position' in line:
                continue
            data = dict((key, line[header[key]]) for key in eric_features)
            data['tot_reads'] = str(int(data['spanningreads']) + int(data['crossingreads']))
            data['gene_start1'] = gene_model[data['GeneName1']].start
            data['gene_end2'] = gene_model[data['GeneName2']].end
            general = [data[general_format[pagasus_feature]] for pagasus_feature in pagasus_features]
            writer.writerow(general)

if __name__ == '__main__':
    main()
