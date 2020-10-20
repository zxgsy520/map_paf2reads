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

__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
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
        if type(line) == type(b''):
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
        if type(line) == type(b''):
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


def read_paf(file):

    if file == "":
        fp = read_stdin('\t')
    elif file.endswith(".paf"):
        fp = read_tsv(file, '\t')
    else:
        raise Exception("%r file format error" % file)

    reads = set()

    for line in fp:
        reads.add(line[0])

    return reads


def paf2reads(pdfile, file):

    reads = read_paf(pdfile)

    if file.endswith(".fastq") or file.endswith(".fq") or file.endswith(".fastq.gz") or file.endswith(".fq.gz"):
        fh = read_fastq(file)
    elif file.endswith(".fasta") or file.endswith(".fa") or file.endswith(".fasta.gz") or file.endswith(".fa.gz"):
        fh = read_fasta(file)
    else:
        raise Exception("%r file format error" % file)

    for line in fh:
        if line[0] not in reads:
            continue
        print(">%s\n%s" % (line[0], line[1]))

    return 0


def add_hlep_args(parser):

    parser.add_argument("-i", "--input", metavar='FILE', type=str, default="",
        help="Input files in paf")
    parser.add_argument("-r", "--read", metavar='FILE', type=str, required=True,
        help="Input files in fastq or fasta format")

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
    paf2reads  Choose appropriate reads for subsequent correction

attention:
    paf2reads -i input.paf -r reads.fastq  >out_reads.fa
    paf2reads -i input.sam -g genome.fa >out_reads.fa
    minimap2 -secondary=no -c -x map-ont genome.fa reads.fastq| paf2reads -r reads.fastq >out_reads.fa

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    paf2reads(args.input, args.read)


if __name__ == "__main__":

    main()
