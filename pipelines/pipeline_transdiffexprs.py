"""========================================================
RNA-Seq pseudocounting pipeline for differential expression
===========================================================

The RNA-Seq differential expression pipeline performs differential
expression analysis using pseudocounting methods. It requires three inputs:
   1. A geneset in :term:`gtf` formatted file
   2. Unaligned reads in :term:`fastq` formatted files
   3. Design files as :term:`tsv`-separated format

This pipeline works on a single genome.

Overview
========

The pipeline performs the following 

   * Gene expression estimates (TPM and counts) at the transcript and
     gene level. The following alignment-free expression estimation
     methods are implemented:
      * kallisto_
      * salmon_
      * sailfish_

   * Perform differential expression analysis using DeSeq2
   
Installation
============

If you are on a mac then you will need to also install the R dependancy wasabi
   
Usage
=====

Configuration
-------------

The pipeline requires a pipeline.yml configuration file. This is located
within the pipeline_transdiffexpress/ directory.

Input
-----

Reads
+++++

Reads are imported by placing :term:`fastq` formatted files in the :term:`working directory`.

The default file format assumes the following convention::
   <samplename>.fastq.gz (fastq.1.gz (and fastq.2.gz for second read of paired data) are also accepted for raw reads)
   
   
Geneset
++++++++
The Geneset is specified by the "geneset" parameter

Design matrices
+++++++++++++++
Design matrices are imported by placing :term:`tsv` formatted files
into the :term:`working directory`. A design matrix describes the
experimental design to test. The design files should be named
design*.tsv. An example can be found in the pipeline_transdiffexprs/ directory.

Each design file has at leasr four columns but may contain any number
of columns after the 'pair' column:

      track   include group   pair
      CW-CD14-R1      0       CD14    1
      CW-CD14-R2      0       CD14    1
      CW-CD14-R3      1       CD14    1
      CW-CD4-R1       1       CD4     1
      FM-CD14-R1      1       CD14    2
      FM-CD4-R2       0       CD4     2
      FM-CD4-R3       0       CD4     2
      FM-CD4-R4       0       CD4     2
      
track
     name of track - should correspond to a sample name.
include
     flag to indicate whether or not to include this data
group
     group indicator - experimental group
pair
     pair that sample belongs to (for paired tests) - set to 0 if the
     design is not paired.

Requirements
------------

The pipeline requires installation using conda and instructions are set out in the main repository README.


Pipeline output
===============

Quantification
--------------

The quantification estimates from each method are outputted to:
[method].dir/[sample]/[level].tsv.gz, where [method] is the quantification method, [sample] is the sample
name and [level] is the feature level (transcript or gene)

Each tool also generates specific log files etc which are outputted, along with the raw quantification outfile in the directory:
[method].dir/[sample]

For each method, the merged counts are outputted to:
[method].dir/[level].tsv.gz

Differential gene expression results
-------------------------------------

Results are stored per method in subdirectories
such as :file:`deseq.dir`

Plots from the differential expression analyses are also contained
within the directories.


Code
====

"""


# Load modules
from ruffus import *

import sys
import os
import re
import glob
import pandas as pd
import sqlite3
import cgat.GTF as GTF
import cgatcore.iotools as iotools
import cgatcore.experiment as E
from cgatcore import pipeline as P

import tasks.mapping as mapping
import tasks.rnaseq as rnaseq


# load options from the config file
P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

PARAMS = P.PARAMS

@mkdir('geneset.dir')
@transform(PARAMS['geneset'],
           regex("(\S+).gtf.gz"),
           r"geneset.dir/\1.fa")
def buildReferenceTranscriptome(infile, outfile):
    '''
    Builds a reference transcriptome from the provided GTF geneset - generates
    a fasta file containing the sequence of each feature labelled as
    "exon" in the GTF.
    --fold-at specifies the line length in the output fasta file
    Parameters
    ----------
    infile: str
        path to the GTF file containing transcript and gene level annotations
    genome_dir: str
        :term: `PARAMS` the directory of the reference genome
    genome: str
        :term: `PARAMS` the filename of the reference genome (without .fa)
    outfile: str
        path to output file
    '''

    genome_file = os.path.abspath(
        os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fa"))

    statement = '''
    zcat < %(infile)s |
    awk '$3=="exon"'|
    cgat gff2fasta
    --is-gtf --genome-file=%(genome_file)s --fold-at=60 -v 0
    --log=%(outfile)s.log > %(outfile)s;
    samtools faidx %(outfile)s
    '''

    P.run(statement)


@transform(buildReferenceTranscriptome,
           suffix(".fa"),
           ".kallisto.index")
def buildKallistoIndex(infile, outfile):
    '''
    Builds a kallisto index for the reference transcriptome
    Parameters
    ----------
    infile: str
       path to reference transcriptome - fasta file containing transcript
       sequences
    kallisto_kmer: int
       :term: `PARAMS` kmer size for Kallisto.  Default is 31.
       Kallisto will ignores transcripts shorter than this.
    outfile: str
       path to output file
    '''

    job_memory = "12G"

    statement = '''
    kallisto index -i %(outfile)s -k %(kallisto_kmer)s %(infile)s
    '''

    P.run(statement)


@transform(buildReferenceTranscriptome,
           suffix(".fa"),
           ".salmon.index")
def buildSalmonIndex(infile, outfile):
    '''
    Builds a salmon index for the reference transriptome
    Parameters
    ----------
    infile: str
       path to reference transcriptome - fasta file containing transcript
       sequences
    salmon_kmer: int
       :term: `PARAMS` kmer size for sailfish.  Default is 31.
       Salmon will ignores transcripts shorter than this.
    salmon_index_options: str
       :term: `PARAMS` string to append to the salmon index command to
       provide specific options e.g. --force --threads N
    outfile: str
       path to output file
    '''

    job_memory = "6G"
    # need to remove the index directory (if it exists) as ruffus uses
    # the directory timestamp which wont change even when re-creating
    # the index files
    statement = '''
    rm -rf %(outfile)s;
    salmon index %(salmon_index_options)s -t %(infile)s -i %(outfile)s
    -k %(salmon_kmer)s
    '''

    P.run(statement)


@originate("transcript2geneMap.tsv")
def getTranscript2GeneMap(outfile):
    ''' Extract a 1:1 map of transcript_id to gene_id from the geneset '''

    iterator = GTF.iterator(iotools.open_file(PARAMS['geneset']))
    transcript2gene_dict = {}

    for entry in iterator:

        # Check the same transcript_id is not mapped to multiple gene_ids!
        if entry.transcript_id in transcript2gene_dict:
            if not entry.gene_id == transcript2gene_dict[entry.transcript_id]:
                raise ValueError('''multipe gene_ids associated with
                the same transcript_id %s %s''' % (
                    entry.gene_id,
                    transcript2gene_dict[entry.transcript_id]))
        else:
            transcript2gene_dict[entry.transcript_id] = entry.gene_id

    with iotools.open_file(outfile, "w") as outf:
        outf.write("transcript_id\tgene_id\n")
        for key, value in sorted(transcript2gene_dict.items()):
            outf.write("%s\t%s\n" % (key, value))

DATADIR = "."

SEQUENCESUFFIXES = ("*.fastq.1.gz",
                    "*.fastq.gz",
                    "*.sra")
SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])

# enable multiple fastqs from the same sample to be analysed together
if "merge_pattern_input" in PARAMS and PARAMS["merge_pattern_input"]:
    SEQUENCEFILES_REGEX = regex(
        r"%s/%s.(fastq.1.gz|fastq.gz|sra)" % (
            DATADIR, PARAMS["merge_pattern_input"].strip()))

    # the last expression counts number of groups in pattern_input
    SEQUENCEFILES_KALLISTO_OUTPUT = [
        r"kallisto.dir/%s/transcripts.tsv.gz" % (
            PARAMS["merge_pattern_output"].strip()),
        r"kallisto.dir/%s/genes.tsv.gz" % (
            PARAMS["merge_pattern_output"].strip())]

    SEQUENCEFILES_SALMON_OUTPUT = [
        r"salmon.dir/%s/transcripts.tsv.gz" % (
            PARAMS["merge_pattern_output"].strip()),
        r"salmon.dir/%s/genes.tsv.gz" % (
            PARAMS["merge_pattern_output"].strip())]

else:
    SEQUENCEFILES_REGEX = regex(
        "(\S+).(fastq.1.gz|fastq.gz|sra)")

    SEQUENCEFILES_KALLISTO_OUTPUT = [
        r"kallisto.dir/\1/transcripts.tsv.gz",
        r"kallisto.dir/\1/genes.tsv.gz"]

    SEQUENCEFILES_SALMON_OUTPUT = [
        r"salmon.dir/\1/transcripts.tsv.gz",
        r"salmon.dir/\1/genes.tsv.gz"]


@follows(mkdir("kallisto.dir"))
@collate(SEQUENCEFILES,
         SEQUENCEFILES_REGEX,
         add_inputs(buildKallistoIndex, getTranscript2GeneMap),
         SEQUENCEFILES_KALLISTO_OUTPUT)
def runKallisto(infiles, outfiles):
    '''
    Computes read counts across transcripts and genes based on a fastq
    file and an indexed transcriptome using Kallisto.
    Runs the kallisto "quant" function across transcripts with the specified
    options.  Read counts across genes are counted as the total in all
    transcripts of that gene (based on the getTranscript2GeneMap table)
    Parameters
    ----------
    infiles: list
        list with three components
        0 - list of strings - paths to fastq files to merge then quantify
        across using Kallisto
        1 - string - path to Kallisto index file
        2 - string - path totable mapping transcripts to genes
    kallisto_threads: int
       :term: `PARAMS` the number of threads for Kallisto
    kallisto_memory: str
       :term: `PARAMS` the job memory for Kallisto
    kallisto_options: str
       :term: `PARAMS` string to append to the Kallisto quant command to
       provide specific
       options, see https://pachterlab.github.io/kallisto/manual
    kallisto_bootstrap: int
       :term: `PARAMS` number of bootstrap samples to run.
       Note, you need to bootstrap for differential expression with sleuth
       if there are no technical replicates. If you only need point estimates,
       set to 1.  Note that bootstrap must be set to at least 1
    kallisto_fragment_length: int
       :term: `PARAMS` Fragment length for Kallisto, required for single end
       reads only
    kallisto_fragment_sd: int
       :term: `PARAMS` Fragment length standard deviation for Kallisto,
       required for single end reads only.
    outfiles: list
       paths to output files for transcripts and genes
    '''

    # TS more elegant way to parse infiles and index?
    fastqfile = [x[0] for x in infiles]
    index = infiles[0][1]
    transcript2geneMap = infiles[0][2]

    transcript_outfile, gene_outfile = outfiles
    Quantifier = rnaseq.KallistoQuantifier(
        infile=fastqfile,
        transcript_outfile=transcript_outfile,
        gene_outfile=gene_outfile,
        annotations=index,
        job_threads=PARAMS["kallisto_threads"],
        job_memory=PARAMS["kallisto_memory"],
        options=PARAMS["kallisto_options"],
        bootstrap=PARAMS["kallisto_bootstrap"],
        fragment_length=PARAMS["kallisto_fragment_length"],
        fragment_sd=PARAMS["kallisto_fragment_sd"],
        transcript2geneMap=transcript2geneMap)

    Quantifier.run_all()


@follows(mkdir("salmon.dir"))
@collate(SEQUENCEFILES,
         SEQUENCEFILES_REGEX,
         add_inputs(buildSalmonIndex, getTranscript2GeneMap),
         SEQUENCEFILES_SALMON_OUTPUT)
def runSalmon(infiles, outfiles):
    '''
    Computes read counts across transcripts and genes based on a fastq
    file and an indexed transcriptome using Salmon.
    Runs the salmon "quant" function across transcripts with the specified
    options.  Read counts across genes are counted as the total in all
    transcripts of that gene (based on the getTranscript2GeneMap table)
    Parameters
    ----------
    infiles: list
        list with three components
        0 - list of strings - paths to fastq files to merge then quantify
        across using sailfish
        1 - string - path to sailfish index file
        2 - string - path to table mapping transcripts to genes
    salmon_threads: int
       :term: `PARAMS` the number of threads for salmon
    salmon_memory: str
       :term: `PARAMS` the job memory for salmon
    salmon_options: str
       :term: `PARAMS` string to append to the salmon quant command to
       provide specific
       options, see http://sailfish.readthedocs.io/en/master/salmon.html
    salmon_bootstrap: int
       :term: `PARAMS` number of bootstrap samples to run.
       Note, you need to bootstrap for differential expression with sleuth
       if there are no technical replicates. If you only need point estimates,
       set to 1.
    salmon_libtype: str
       :term: `PARAMS` salmon library type
       as for sailfish - use
       http://sailfish.readthedocs.io/en/master/library_type.html#fraglibtype
    outfiles: list
       paths to output files for transcripts and genes
    '''
    fastqfile = [x[0] for x in infiles]
    index = infiles[0][1]
    transcript2geneMap = infiles[0][2]

    transcript_outfile, gene_outfile = outfiles
    Quantifier = rnaseq.SalmonQuantifier(
        infile=fastqfile,
        transcript_outfile=transcript_outfile,
        gene_outfile=gene_outfile,
        annotations=index,
        job_threads=PARAMS["salmon_threads"],
        job_memory=PARAMS["salmon_memory"],
        options=PARAMS["salmon_options"],
        bootstrap=PARAMS["salmon_bootstrap"],
        libtype=PARAMS['salmon_libtype'],
        transcript2geneMap=transcript2geneMap)

    Quantifier.run_all()


def full():
    ''' dummy task for full ruffus tasks'''
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))


   