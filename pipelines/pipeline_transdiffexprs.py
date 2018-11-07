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
from ruffus.combinatorics import *

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
import tasks.tracks as tracks


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


DATADIR = "."

SEQUENCESUFFIXES = ("*.fastq.1.gz",
                    "*.fastq.gz")
SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])

SEQUENCEFILES_REGEX = regex(
        "(\S+).(fastq.1.gz|fastq.gz)")


@follows(mkdir("kallisto.dir"))
@collate(SEQUENCEFILES,
         SEQUENCEFILES_REGEX,
         add_inputs(buildKallistoIndex),
         [r"kallisto.dir/\1/abundance.h5",r"kallisto.dir/\1/abundance.tsv"])
def run_kallisto(infiles, outfiles):
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

    fastqfile = infiles[0][0]
    index = infiles[0][1]

	# check for paired end files and overwrite fastqfile if True
    if fastqfile.endswith(".fastq.1.gz"):
    	bn = P.snip(fastqfile, ".fastq.1.gz")
    	infile1 = "%s.fastq.1.gz" % bn
    	infile2 = "%s.fastq.2.gz" % bn
    	if not os.path.exists(infile2):
    		raise ValueError(
    			"can not find paired ended file "
    			"'%s' for '%s'" % (infile2, infile))
    	fastqfile = infile1 + " " + infile2

    outfile = outfiles[0].replace("abundance.h5","")
    statement = """kallisto quant -i %(index)s %(kallisto_options)s -o %(outfile)s %(fastqfile)s"""
    P.run(statement)


###################################################
# Generate transcript2gene infomation
###################################################

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


###################################################
###################################################
# Create quantification targets
###################################################

@collate(run_kallisto,
         regex("(\S+).dir/(\S+)/abundance.h5"),
         [r"\1.dir/counts.tsv.gz"])
def merge_tpm(infiles, outfile):
    ''' merge counts across all samples - this is not relly used
        within the downstream tasks of the pipeline and is there
        for reference only'''

    transcript_infiles = [x[1] for x in infiles]

    final_df = pd.DataFrame()

    for infile in transcript_infiles:
    	path = os.path.normpath(infile)
    	folder_name = path.split("/")[1]

    	tmp_df = pd.read_table(infile, sep="\t", index_col=0)
    	# Only select tpm values and rename with name of folder
    	tmp_df = tmp_df[["tpm"]]
    	tmp_df.columns = [folder_name]
    	final_df = final_df.merge(tmp_df, how="outer",  left_index=True, right_index=True)
    final_df = final_df.round()
    final_df.sort_index(inplace=True)
    final_df.to_csv(outfile[0], sep="\t", compression="gzip")


###################################################
# Differential Expression
###################################################

Sample = tracks.AutoSample
DESIGNS = tracks.Tracks(Sample).loadFromDirectory(
    glob.glob("design*.tsv"), "design(\S+).tsv")

@mkdir("DEresults.dir/deseq2")
@product(merge_tpm,
         formatter(".*/(?P<QUANTIFIER>\S+).dir/transcripts.tsv.gz"),
         ["design%s.tsv" % x.asFile() for x in DESIGNS],
         formatter(".*/design(?P<DESIGN>\S+).tsv$"),
         ["DEresults.dir/deseq2/{QUANTIFIER[0][0]}_{DESIGN[1][0]}_transcripts_results.tsv",
          "DEresults.dir/deseq2/{QUANTIFIER[0][0]}_{DESIGN[1][0]}_genes_results.tsv"],
         "{DESIGN[1][0]}")
def runDESeq2(infiles, outfiles, design_name):
    ''' run DESeq2 to identify differentially expression transcripts/genes'''

    design_name = design_name.lower()
    counts, design = infiles
    transcripts, genes = counts
    transcript_out, gene_out = outfiles

    transcript_prefix = P.snip(transcript_out, ".tsv")
    transcript_log = transcript_prefix + ".log"

    gene_prefix = P.snip(gene_out, ".tsv")
    gene_log = gene_prefix + ".log"

    model = PARAMS.get('deseq2_model%s' % design_name, None)
    contrast = PARAMS.get('deseq2_contrast%s' % design_name, None)
    refgroup = PARAMS.get('deseq2_refgroup%s' % design_name, None)

    if model is None:
        raise ValueError("deseq2_model{} is not specified".format(
            design_name))
    if contrast is None:
        raise ValueError("deseq2_contrast{} is not specified".format(
            design_name))
    if refgroup is None:
        raise ValueError("deseq2_refgroup{} is not specified".format(
            design_name))
# in future when it is a package run as: python -m pipelines.tasks.counts2table
    statement = '''
    python ../cgat-showcase/pipelines/tasks/counts2table.py
    --tag-tsv-file=%(transcripts)s
    --design-tsv-file=%(design)s
    --method=deseq2
    --de-test=%(deseq2_detest)s
    --output-filename-pattern=%(transcript_prefix)s
    --model=%(model)s
    --contrast=%(contrast)s
    --reference-group=%(refgroup)s
    --fdr=%(deseq2_fdr)s
    --log=%(transcript_log)s
    -v 0
    > %(transcript_out)s;
    '''
    P.run(statement)

    statement = '''
    python ../cgat-showcase/pipelines/tasks/counts2table.py
    --tag-tsv-file=%(genes)s
    --design-tsv-file=%(design)s
    --method=deseq2
    --de-test=%(deseq2_detest)s
    --output-filename-pattern=%(gene_prefix)s
    --model=%(model)s
    --contrast=%(contrast)s
    --reference-group=%(refgroup)s
    --fdr=%(deseq2_fdr)s
    --log=%(gene_log)s
    -v 0
    > %(gene_out)s;
    '''

    P.run(statement)


@mkdir("DEresults.dir/sleuth")
@product(merge_tpm,
         formatter(
             ".*/(?P<QUANTIFIER>(kallisto|salmon|sailfish)).dir/transcripts.tsv.gz"),
         ["design%s.tsv" % x.asFile() for x in DESIGNS],
         formatter(".*/design(?P<DESIGN>\S+).tsv$"),
         ["DEresults.dir/sleuth/{QUANTIFIER[0][0]}_{DESIGN[1][0]}_transcripts_results.tsv",
          "DEresults.dir/sleuth/{QUANTIFIER[0][0]}_{DESIGN[1][0]}_genes_results.tsv"],
         "{DESIGN[1][0]}",
         "{QUANTIFIER[0][0]}")
def runSleuth(infiles, outfiles, design_name, quantifier):
    ''' run sleuth to identify differentially expression transcripts/genes'''

    design_name = design_name.lower()
    counts, design = infiles
    transcripts, genes = counts
    transcript_out, gene_out = outfiles

    transcript_prefix = P.snip(transcript_out, ".tsv")
    transcript_log = transcript_prefix + ".log"

    gene_prefix = P.snip(gene_out, ".tsv")
    gene_log = gene_prefix + ".log"

    model = PARAMS['sleuth_model%s' % design_name]
    E.info(model)
    reduced_model = PARAMS['sleuth_reduced_model%s' % design_name]

    contrast = PARAMS['sleuth_contrast%s' % design_name]
    refgroup = PARAMS['sleuth_refgroup%s' % design_name]
    detest = PARAMS['sleuth_detest']
    transcripts = os.path.join("geneset.dir",
                               P.snip(PARAMS['geneset'], ".gtf.gz") + ".fa")

    # to estimate sleuth memory, we need to know the number of
    # samples, transcripts and boostraps
    number_transcripts = 0
    with iotools.open_file(transcripts, "r") as inf:
        for line in inf:
            if line.startswith(">"):
                number_transcripts += 1

    Design = Expression.ExperimentalDesign("design%s.tsv" % design_name)
    number_samples = sum(Design.table['include'])

    job_memory = rnaseq.estimateSleuthMemory(
        PARAMS["%(quantifier)s_bootstrap" % locals()],
        number_samples, number_transcripts)

    statement = '''
    python -m cgatpipelines.tasks.counts2table
    --design-tsv-file=%(design)s
    --output-filename-pattern=%(transcript_prefix)s
    --log=%(transcript_log)s
    --method=sleuth
    --fdr=%(sleuth_fdr)s
    --model=%(model)s
    --contrast=%(contrast)s
    --sleuth-counts-dir=%(quantifier)s.dir
    --reference-group=%(refgroup)s
    --de-test=%(detest)s
    '''
    if detest == "lrt":
        statement += '''
        --reduced-model=%(reduced_model)s
        '''
    statement += '''
    -v 0
    >%(transcript_out)s
    '''

    P.run(statement)

    if PARAMS['sleuth_genewise']:

        assert PARAMS['sleuth_gene_biomart'], (
            "Must provide a biomart (see pipeline.yml)")

        # gene-wise sleuth seems to be even more memory hungry!
        # Use 2 * transcript memory estimate
        job_memory = rnaseq.estimateSleuthMemory(
            PARAMS["%(quantifier)s_bootstrap" % locals()],
            2 * number_samples, number_transcripts)

        statement = '''
        python -m cgatpipelines.tasks.counts2table
        --design-tsv-file=%(design)s
        --output-filename-pattern=%(gene_prefix)s
        --log=%(gene_log)s
        --method=sleuth
        --fdr=%(sleuth_fdr)s
        --model=%(model)s
        --contrast=%(contrast)s
        --sleuth-genewise
        --sleuth-counts-dir=%(quantifier)s.dir
        --reference-group=%(refgroup)s
        --gene-biomart=%(sleuth_gene_biomart)s
        --de-test=%(detest)s
        '''
        if detest == "lrt":
            statement += '''
            --reduced-model=%(reduced_model)s
            '''
        statement += '''
        -v 0
        >%(transcript_out)s
        '''

        P.run(statement)

@mkdir("DEresults.dir/deseq2")
@transform(merge_tpm,
           regex("(\S+).dir/transcripts.tsv.gz"),
           [r"DEresults.dir/deseq2/\1_normalised_transcripts_expression.tsv.gz",
            r"DEresults.dir/deseq2/\1_normalised_genes_expression.tsv.gz"])
def getDESeqNormExp(infiles, outfiles):
    ''' Use the Deseq2 size factors method to obtain normalised
    expression values for summary plots '''
    # currently DESeq expression factors is not working
    # to normalise the expression values. In some workflows
    # the columns produce infinity values. Have defaulted to
    # total column until issue is identified

    transcripts_inf, genes_inf = infiles
    transcripts_outf, genes_outf = outfiles

    normalisation_method = "total-column"

    rnaseq.normaliseCounts(
        transcripts_inf, transcripts_outf, normalisation_method)
    rnaseq.normaliseCounts(
        genes_inf, genes_outf, normalisation_method)


@mkdir("DEresults.dir/sleuth")
@collate(
    run_kallisto,
    formatter(
        "(?P<QUANTIFIER>(kallisto|salmon|sailfish)).dir/(\S+)/transcripts.tsv.gz"),
    [r"DEresults.dir/sleuth/{QUANTIFIER[0]}_normalised_transcripts_expression.tsv.gz",
     r"DEresults.dir/sleuth/{QUANTIFIER[0]}_normalised_genes_expression.tsv.gz"],
    r"{QUANTIFIER[0]}")
def getSleuthNormExp(infiles, outfiles, quantifier):
    ''' get the Normalised expression from the quantification tools
    which we will run sleuth from'''

    t2gMap = infiles[0][1]
    transcript_infiles = [x[0][0] for x in infiles]
    transcripts_outf, genes_outf = outfiles

    if quantifier == "kallisto":
        basename = "abundance.h5.tsv"
        column = "tpm"
    else:
        raise ValueError("using unknown quantifier!")

    rnaseq.getAlignmentFreeNormExp(
        transcript_infiles, basename, column,
        transcripts_outf, genes_outf, t2gMap)


# Define the task for differential expression and normalisation
DETARGETS = []
NORMTARGETS = []
mapToDETargets = {'deseq2': (runDESeq2,),
                  'sleuth': (runSleuth,)}

mapToNormTargets = {'deseq2': (getDESeqNormExp, ),
                    'sleuth': (getSleuthNormExp,)}

for x in P.as_list(PARAMS["de_tools"]):
    DETARGETS.extend(mapToDETargets[x])
    if x in mapToNormTargets:
        NORMTARGETS.extend(mapToNormTargets[x])


@follows(*DETARGETS)
def differentialExpression():
    ''' dummy task to define upstream differential expression tasks'''


@follows(*NORMTARGETS)
def NormaliseExpression():
    ''' dummy task to define upstream normalisation tasks'''


###################################################
# Summary plots
###################################################

@mkdir("summary_plots")
@product(NORMTARGETS,
         formatter(
             "DEresults.dir/(?P<DETOOL>\S+)/(?P<QUANTIFIER>\S+)_normalised_transcripts_expression.tsv.gz"),
         ["design%s.tsv" % x.asFile() for x in DESIGNS],
         formatter(".*/design(?P<DESIGN>\S+).tsv$"),
         ["summary_plots/{DETOOL[0][0]}_{QUANTIFIER[0][0]}_{DESIGN[1][0]}_transcripts_plots.log",
          "summary_plots/{DETOOL[0][0]}_{QUANTIFIER[0][0]}_{DESIGN[1][0]}_genes_plots.log"])
def expressionSummaryPlots(infiles, logfiles):
    ''' make summary plots for expression values for each design file'''

    expression_infs, design_inf = infiles
    transcript_inf, gene_inf = expression_infs
    transcript_log, gene_log = logfiles

    job_memory = "10G"

    if not os.path.exists(os.path.dirname(gene_log)):
        os.mkdir(os.path.dirname(gene_log))

    rnaseq.makeExpressionSummaryPlots(
        transcript_inf, design_inf, transcript_log, submit=True,
        job_memory=job_memory)

    rnaseq.makeExpressionSummaryPlots(
        gene_inf, design_inf, gene_log, submit=True,
        job_memory=job_memory)

###################################################
# target functions for code execution             #
###################################################

def full():
    ''' dummy task for full ruffus tasks'''
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))


   