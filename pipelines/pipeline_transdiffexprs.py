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

import module_transdiffexprs as tdexp


# load options from the config file
P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])


def full():
    ''' dummy task for full ruffus tasks'''
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))


   