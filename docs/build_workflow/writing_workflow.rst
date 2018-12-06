.. _writing_workflow-tutorial:


=============================
Writing a pipeline - Tutorial
=============================

Before attempting this tutorial please make sure you have cgat-showcase installed correctly. please see here (see :ref:`getting_started-Installation`) for instructions.

This tutorial is intended as a walkthrough to show you the processes inviolved in writing the cgat-showcase transdiffexprs pipeline. It wont go through all of the functions in the pipeline but it should give you an idea of how to construct a pipeline using our cgat-core python library. Before attempting this tutorial, it may be best to visit cgat-core `documentation <https://cgat-core.readthedocs.io/en/latest/defining_workflow/Tutorial.html>`_  for an idea of how to construct a very simple pipeline.

Tutorial
========

1. The problem and the solution
-------------------------------

The first step in any pipeline is to think about the problem that you are trying solve. 
For example: 

**The problem:** I need to perform differential expression analysis quickly and consistently across all my samples. I also need the solution to scale and be ran on my mac or districuted across a cluster.
**The solution:** Build a cgat-core workflow that wraps pseduocounting tools (kallisto) and differential expression analysis tools (deseq2). 

Then the next step is to think about the steps needed to reach the solution.

For example:

   1. Build a reference transcriptome from a gtf
   2. Build a kallisto index of that transcriptome
   3. Perform pseudocounting of the fasta sequencing data
   4. Run differential expression analysis using deseq2
   5. Generate a report of the pseudocounting quality
   6. Generate a report displaying the differential expression output

Then the next step is to begin to code the pipeline. There are certain core modules that will help with the building of our workflows and they are
described below.

2. Setting up the pipeline
--------------------------

First navigate to a directory where you would like to start building the pipeline then create a directory called cgat-showcase::

    mkdir cgat-showcase && mkdir cgat-showcase/cgatshowcase && touch cgat-showcase/cgatshowcase/pipeline_transcriptdiffexprs.py && mkdir cgat-showcase/cgatshowcase/pipeline_transcriptdiffexprs && touch cgat-showcase/cgatshowcase/ModuleTrans.py

This will create a directory with the following input::

   |-- cgat-showcase
      |-- cgatshowcase
         -- pipeline_transdiffexprs.py
         -- ModuleTrans.py
         |-- pipeline_transdiffexprs

The layout has the following components::

pipeline_transdiffexprs.py
   This is the file that will contain all of the ruffus workflows, the file needs
   the format pipeline_<name>.py
pipeline_transdiffexprs/
   Directory containing the configuration yml file. The directory needs to be named
   the same as the pipeline_<name>.py file
ModuleTrans.py
   This file will contain functions that will be imported into the main ruffus
   workflow file, pipeline_transdiffexprs.py

3. Creating basic utility functions
-----------------------------------

Open the cgat-shocase/cgatshowcase/pipeline_transdiffexprs.py file and start populating code with basic utility functions to help
with building workflows.

The code begins with a doc
string detailing the pipeline functionality.You should use this section to document your
pipeline.::

    '''This pipeline is a test and this is where the documentation goes '''


**Config parser:** This code will help with parsing the pipeline.yml file::


   # load options from the config file
   P.get_parameters(
       ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
        "../pipeline.yml",
        "pipeline.yml"])

   PARAMS = P.PARAMS


**Commandline parser:** This bit of code allows pipeline to parse arguments and run the pipleine using the 
`cgatshowcase` command::

    def main(argv=None):
	if argv is None:
	    argv = sys.argv
	P.main(argv)


    if __name__ == "__main__":
	sys.exit(P.main(sys.argv))    

4. Create the ruffus functions
------------------------------

The first ruffus function will build a reference transcriptome so kallisto can build an index.

For this to be acomplished the following code is written::

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

Subsequent functions can then be written to form the workflow of the pipeline. More infomation on how to build ruffus pipelines can be found on the `ruffus  <www.ruffus.org.uk>`_ documentation and `cgat-core <https://cgat-core.readthedocs.io/en/latest/>`_ documentation.
Please see the source code in pipeline_tranmsdiffexprs.py for more infomation on the specific ruffus functions that have been written for the current pipeline. 

5. Creating a pipeline.yml
--------------------------

There are different ways to pass configuration values to modify the output of a pipeline. For example, the threshold of the padj filtering in deseq2 can be modified. In order to 
achieve this we modify a configuration file, pipeline.yml. This needs to be created in the pipeline_transdiffexprs/ directory.::

   touch &&  cgat-showcase/cgatshowcase/pipeline_transcriptdiffexprs/pipeline.yml

The configuration file then needs to be written so that the yml parser can pick up different options. For example,::

   kallisto:
    kmer: 31

The kmer will then have the value of 31 and can be passed into then pipeline. This is accessed within the pipeline as PARAMS['kallisto_kmer'] and will have the value 31.

6. Generating a report
----------------------

In order to present our results in a visually appealing manner, we use multiQC reports to display generic sequencing metrics and `rmarkdown <https://rmarkdown.rstudio.com/>`_ (or ipython if you prefer) for generating
bespoke reports.

In order to write a report you will need to generate an rmarkdown document folder as follows::

   mkdir cgat-showcase/cgatshowcase/pipeline_docs/pipeline_transdiffexprs/R_report

Rmarkdown reports are best constructed using rstudio as they can be written dynamically, then best tested within the build functionality within rstudio.

In order to build a Rmarkdown website we usually follow the Rmarkdown `tutorial <https://rmarkdown.rstudio.com/lesson-1.html>`_  

The basic requirements for building a website are::

   index.Rmd
   _site.yml
   R_report.Rproj

Once the report has been developed in rstudio and tested then the pipeline requires a function to copy and then render the report::

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

This function copies the R_report.dir from the main repository and then calles rmarkdown to render the website according to the _site.yml configuration file. 

