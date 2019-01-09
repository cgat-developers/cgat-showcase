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


   * Perform differential expression analysis using DeSeq2
   
Installation
============

If you are on a mac then you will need to also install the R dependancy wasabi
   
Usage
=====

Configuration
-------------

The pipeline requires a pipeline.yml configuration file. This is located
within the ?? directory.

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
design*.tsv. An example can be found ???.

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




Code
====

"""


# Load modules
from ruffus import *
from ruffus.combinatorics import *

import sys
import os
import pandas as pd
import cgatcore.iotools as iotools
import cgatcore.experiment as E
from cgatcore import pipeline as P
import cgatshowcase.ModuleTransdiffexprs as ModuleTransdiffexprs

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

    #####################################
    # Domonstraion of logging and warning
    #####################################
    # infomation can be output to the standard out using E.info()
    E.info(genome_file)
    # warning messages can also be displayed to stdout using E.warn()
    E.warn("This is an example of a message used to warn people")


    statement = '''
    zcat < %(infile)s |
    cgat gff2fasta
    --feature=exon --is-gtf --genome-file=%(genome_file)s --fold-at=60 -v 0
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


#################################################
# Run alignment free quantification - kallisto
#################################################


@transform(tuple([suffix_name for suffix_name in ("*.fastq.1.gz",
                                                  "*.fastq.gz")]),
           regex("(\S+).(fastq.1.gz|fastq.gz)"),
           add_inputs(buildKallistoIndex),
           r"kallisto.dir/\1/abundance.h5")
def runKallisto(infiles, outfile):
    '''
    Computes read counts across transcripts and genes based on a fastq
    file and an indexed transcriptome using Kallisto.
    Runs the kallisto "quant" function across transcripts with the specified
    options.  Read counts across genes are counted as the total in all
    transcripts of that gene.
    '''

    fastqfile = infiles[0]
    index = infiles[1]

    # check for paired end files and overwrite fastqfile if True
    fastqfile = ModuleTransdiffexprs.check_paired_end(fastqfile)

    outfile = outfile.replace("/abundance.h5","")

    statement = '''mkdir %(outfile)s && kallisto quant
                   -i %(index)s
                   -t %(kallisto_threads)s
                   %(kallisto_options)s
                   -o %(outfile)s
                   %(fastqfile)s
                   &> %(outfile)s/kallisto.log
                   > %(outfile)s/kallisto.sdtout'''

    P.run(statement)


###################################################
# Differential Expression
###################################################


@mkdir("DEresults.dir")
@merge(runKallisto,
       "DEresults.dir/counts.csv")
def run_deseq2(infiles, outfiles):
    ''' run DESeq2 to identify differentially expression'''

    # this may be a bug, but when you supply a list of output files
    # using @merge you need to unpack them before running P.run(). I will
    # look into this AC

    R_ROOT = os.path.join(os.path.dirname(__file__), "R")

    if PARAMS["deseq2_detest"] == "lrt":

        statement = '''Rscript %(R_ROOT)s/DESeq2_lrt.R
                       --design=design.tsv
                       --contrast=%(deseq2_contrast)s
                       --refgroup=%(deseq2_control)s
                       --fdr=%(deseq2_fdr)s
                       --biomart=%(deseq2_biomart)s'''

    elif PARAMS['deseq2_detest'] == "wald":
        
        statement = '''Rscript  %(R_ROOT)s/DESeq2_wald.R
                       --design=design.tsv
                       --contrast=%(deseq2_contrast)s
                       --refgroup=%(deseq2_control)s
                       --fdr=%(deseq2_fdr)s
                       --biomart=%(deseq2_biomart)s'''

    P.run(statement)


###################################################
# Convert counts to tpm values
###################################################

# active_if decorators control the flow of a pipeline
@active_if(PARAMS["tpm_run"] == 1)
@transform(run_deseq2,
           suffix(".csv"),
           ".tsv")
def counts2tpm(infile, outfile):
    '''Converts counts to tpm values using a gist ada-pted from slowkow:
       https://gist.github.com/slowkow/c6ab0348747f86e2748b

       Parameters
       ----------
       infiles: output of deseq2 counts and length
       genome: the genome name in ensembl format
       meanfraglength: the mean fragment length - this can be estimated from bioanalyzer
    '''

    R_ROOT = os.path.join(os.path.dirname(__file__), "R")

    # on our cluster there is latency before the deseq2 counts are saved
    # therefore I have added a sleep to make sure output from previous task
    # has closed the connection fully.
    statement = '''sleep 30 && Rscript %(R_ROOT)s/counts2tpm.R
                           --counts=%(infile)s
                           --genome=%(tpm_genome_version)s
                           --meanfraglength=%(tpm_frag_length)s
                           --effectivelength=DEresults.dir/length.csv'''

    P.run(statement)


##################################################
# Add output of counts to a database
##################################################

@transform(counts2tpm, suffix(".tsv"), ".load")
def load_tpm(infile, outfile):
    '''This function uploads a tsv seperated file to an sqlite database
       Parameters
       ----------
       infile: term tsv file containing a table of tpm outputs
       outfile: .load file
    '''

    P.load(infile,outfile)


###################################################
# Generate a report
###################################################

# both multiqc and Rmarkdown
@follows(mkdir("MultiQC_report.dir"))
@follows(run_deseq2)
def run_multiqc_report():
    '''This will generate a mutiqc report '''
    statement = (
        "export LC_ALL=en_GB.UTF-8 && "
        "export LANG=en_GB.UTF-8 && "
        "multiqc . -f && "
        "mv multiqc_report.html MultiQC_report.dir/")
    P.run(statement)


@follows(mkdir("R_report.dir"))
@follows(run_deseq2)
def run_rmarkdown_report():
    '''This will generate a rmarkdown report '''

    report_path = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                               'pipeline_docs',
                                               'pipeline_transdiffexprs',
                                               'R_report.dir'))

    statement = '''cp %(report_path)s/* R_report.dir ; cd R_report.dir ; R -e "rmarkdown::render_site(encoding = 'UTF-8')"'''
    P.run(statement)

###################################################
# target functions for code execution             #
###################################################
@follows(run_multiqc_report, run_rmarkdown_report, counts2tpm)
def full():
    'dummy task for full ruffus tasks'
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))


   
