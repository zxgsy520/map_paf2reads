#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import json
import gzip
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.1.1"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_fasta(file):
    '''Read fasta file'''

    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fasta") or file.endswith(".fa"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = ''

    for line in fp:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if not seq:
            seq += "%s\n" % line.strip(">").split()[0]
            continue
        if line.startswith(">"):
            line = line.strip(">").split()[0]
            seq = seq.split('\n')

            yield [seq[0], seq[1]]
            seq = ''
            seq += "%s\n" % line
        else:
            seq += line

    seq = seq.split('\n')
    if len(seq)==2:
        yield [seq[0], seq[1]]
    fp.close()


def read_fastq(file):
    '''Read fastq file'''

    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fastq") or file.endswith(".fq"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = []

    for line in fp:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if line.startswith('@') and (len(seq)==0 or len(seq)>=5):
            seq.append(line.strip("@").split()[0])
            continue
        if line.startswith('@') and len(seq)==4:
            yield seq
            seq = []
            seq.append(line.strip("@").split()[0])
        else:
            seq.append(line)

    if len(seq)==4:
        yield seq
    fp.close()


def read_tsv(file, sep=None):
    """
    read tsv joined with sep
    :param file: file name
    :param sep: separator
    :return: list
    """
    LOG.info("reading message from %r" % file)

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def read_stdin(sep=None):

    for line in sys.stdin:
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def read_paf(file, identity, coverage):

    if file == "":
        fp = read_stdin('\t')
    elif file.endswith(".paf"):
        fp = read_tsv(file, '\t')
    else:
        raise Exception("%r file format error" % file)

    for line in fp:
        rlen = int(line[6])
        rsite = [int(line[7]), int(line[8])]
        mlen = max(rsite)-min(rsite)
        match = [int(line[9]), int(line[10])]

        similar =min(match) *100.0/max(match)
        lenratio = mlen*100.0/rlen

        if similar<identity or lenratio<coverage:
            continue

        yield line[0], line[2], line[3], similar, lenratio


def sort_start_end(start, end, seqlen):

    nstart, nend = start, end
    if start>=end:
        nstart, nend = send, start
    nstart = nstart - 5
    nend = nend + 5
    if nstart<0:
        nstart = 0
    if nend>seqlen:
        nend = seqlen

    return nstart, nend


def map_paf2reads(pafile, file, stat, identity, coverage):


    data = {}
    n = 0
    fs = open(stat, 'w')

    for seqid, start, end, ident, cover in read_paf(pafile, identity, coverage):
        if seqid not in data:
            data[seqid] = []
            data[seqid].append([seqid, int(start), int(end)])
            nseqid = seqid
            n = 0
        else:
            n += 1
            nseqid = "%s_%s" % (seqid, n)
            data[seqid].append([nseqid, int(start), int(end)])
        fs.write("{0}\t{1:.2f}\t{2:.3f}\n".format(nseqid, ident, cover))
    fs.close()

    if file.endswith(".fastq") or file.endswith(".fq") or file.endswith(".fastq.gz") or file.endswith(".fq.gz"):
        fh = read_fastq(file)
    elif file.endswith(".fasta") or file.endswith(".fa") or file.endswith(".fasta.gz") or file.endswith(".fa.gz"):
        fh = read_fasta(file)
    else:
        raise Exception("%r file format error" % file)

    for line in fh:
        seqid = line[0]
        if seqid not in data:
            continue
        seq = line[1]
        slen = len(seq)

        for line in data[seqid]:
            start, end = sort_start_end(line[1], line[2], slen)
            print(">%s\n%s" % (line[0], seq[start:end]))


def add_hlep_args(parser):

    parser.add_argument("-i", "--input", metavar='FILE', type=str, default="",
        help="Input files in paf")
    parser.add_argument("-r", "--read", metavar='FILE', type=str, required=True,
        help="Input files in fastq or fasta format")
    parser.add_argument("-id", "--identity", metavar='FLOAT', type=float, default=75.0,
        help="Comparison of identity values, default=75.0")
    parser.add_argument("-c", "--coverage", metavar='FLOAT', type=float, default=80.0,
        help="Comparison coverage, default=80.0")
    parser.add_argument("-s", "--stat", metavar='STR', type=str, default='out.stv',
        help="The name of the output file")

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''
name:
    map_paf2reads  Choose appropriate reads for subsequent correction

attention:
    map_paf2reads -i input.paf -r reads.fastq -s identity_coverage.tsv >out_reads.fa
    map_paf2reads -i input.sam -g genome.fa -p name
    minimap2 -secondary=no -c -x map-ont genome.fa reads.fastq|map_paf2reads -r reads.fastq -s identity_coverage.tsv >out_reads.fa

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    map_paf2reads(args.input, args.read, args.stat, args.identity, args.coverage)


if __name__ == "__main__":

    main()
