#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def read_tsv(file, sep='\t'):

    if file.endswith(".gz"):
        ft = gzip.open(file)
    else:
        ft = open(file)

    for line in ft:
        if isinstance(line, bytes):
            line = line.decode("utf-8")

        line = line.strip()
        if not line or line.startswith("#"):
            continue

        yield line.split(sep)
    ft.close()


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
            line = line.decode("utf-8")
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
            line = line.decode("utf-8")
        line = line.strip()

        if not line:
            continue
        if line.startswith('@') and (len(seq)==0 or len(seq)>=5):
            seq.append(line.strip("@").split()[0])
            continue
        if line.startswith('@') and len(seq)==4:
            yield seq[0], seq[1]
            seq = []
            seq.append(line.strip("@").split()[0])
        else:
            seq.append(line)

    if len(seq)==4:
        yield seq[0], seq[1]
    fp.close()


def stat_length(file):

    r = {}

    if file.endswith(".fastq") or file.endswith(".fq") or file.endswith(".fastq.gz") or file.endswith(".fq.gz"):
        fh = read_fastq(file)
    elif file.endswith(".fasta") or file.endswith(".fa") or file.endswith(".fasta.gz") or file.endswith(".fa.gz"):
        fh = read_fasta(file)
    else:
        raise Exception("%r file format error" % file)

    for seqid, seq in fh:
        r[seqid] = len(seq)

    return r


def stat_full_len(file, pafile, identity=80, coverage=80):

    data = {}
    r = stat_length(file)

    for line in read_tsv(pafile):
        if line[5] not in r:
            continue

        qlen = int(line[1])
        qstart, qend = int(line[2]), int(line[3])
        if qstart >= qend:
            qstart, qend = qend, qstart

        if (qend-qstart+1)*100.0/qlen < 75:
            continue

        rstart, rend = int(line[7]), int(line[8])
        if rstart >= rend:
            rstart, rend = rend, rstart

        if (rend-rstart+1)*100/r[line[5]] < coverage:
            continue

        unmap = abs((qend-qstart+1) - (rend-rstart+1))
        if (r[line[5]]-unmap)*100.0/r[line[5]] < identity:
            continue
        if line[5] not in data:
            data[line[5]] = 0
        data[line[5]] += 1
        LOG.info('\t'.join(line))
    print("#Reference id\tReads number")
    for i in r:
        k = 0
        if i in data:
            k = data[i]
        print("{0}\t{1:,.2f}".format(i, k))

    return 0


def add_hlep_args(parser):

    parser.add_argument("reference", metavar='FILE', type=str,
        help="Input files in fastq or fasta format")
    parser.add_argument("-p", "--paf", metavar='FILE', type=str, required=True,
        help="Input files in paf")
    parser.add_argument("-id", "--identity", metavar='FLOAT', type=float, default=80.0,
        help="Comparison of identity values, default=80.0")
    parser.add_argument("-c", "--coverage", metavar='FLOAT', type=float, default=80.0,
        help="Comparison coverage, default=80.0")

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
    stat_full_len.py  Count the number of full-length reads with reference pairs

attention:
    stat_full_len.py ref.fa -p map.paf >map.stat
    

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    stat_full_len(args.reference, args.paf, args.identity, args.coverage)


if __name__ == "__main__":

    main()
