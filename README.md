# map_paf2reads
Processing nanopore amplicon sequencing reads with reference sequence

The alignment files are used to interrupt self-linked and chimeric reads to obtain more full-length amplicon sequences.
## Installation
<pre><code>
wget -c https://github.com/zxgsy520/map_paf2reads/archive/v1.0.0.tar.gz
tar -zxvf v1.0.0.tar.gz
cd map_paf2reads-1.0.0
 chmod 755 map_paf2reads
</code></pre>
or
<pre><code>
git clone https://github.com/zxgsy520/map_paf2reads.git
cd map_paf2reads
chmod 755 map_paf2reads
</code></pre>
## Instructions
<pre><code>
usage: map_paf2reads [-h] [-i FILE] -r FILE [-id FLOAT] [-c FLOAT] [-s STR]

name:
    map_paf2reads  Choose appropriate reads for subsequent correction

attention:
    map_paf2reads -i input.paf -r reads.fastq -s identity_coverage.tsv >out_reads.fa
    map_paf2reads -i input.sam -g genome.fa -p name
    minimap2 -secondary=no -c -x map-ont genome.fa reads.fastq|map_paf2reads -r reads.fastq -s identity_coverage.tsv >out_reads.fa

version: 1.0.0
contact:  Xingguo Zhang <113178210@qq.com>        

optional arguments:
  -h, --help            show this help message and exit
  -i FILE, --input FILE
                        Input files in paf
  -r FILE, --read FILE  Input files in fastq or fasta format
  -id FLOAT, --identity FLOAT
                        Comparison of identity values, default=75.0
  -c FLOAT, --coverage FLOAT
                        Comparison coverage, default=80.0
  -s STR, --stat STR    The name of the output file
</code></pre>
### example
<pre><code>
minimap2 -t 4 -secondary=no -c -x map-ont reference.fasta reads.fastq |map_paf2reads -r reads.fastq -s identity_coverage.tsv >clean_reads.fa
</code></pre>
identity_coverage.tsv: Output the identity and coverage value files of each read
clean_reads.fa:Output trimmed reads
