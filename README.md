# cgat-showcase

cgat-showcase is a repository containing an example pipeline constructed to demonstrate how our [cgat-core](https://github.com/cgat-developers/cgat-core) workflow management system can be used to create common workflows required in bioinformatics analysis.

Within this repository the example [pipeline](https://github.com/cgat-developers/cgat-showcase/blob/master/cgatshowcase/pipeline_transdiffexprs.py) performs pseudoalignment of fastq files
with kallisto and differential expression using deseq2. It can be run locally on your own machine or distributed across a cluster depending on your requirements.

Documentation on how to run this pipeline can be found [here](https://cgat-showcase.readthedocs.io/en/latest/) and documentation on how
to build your own custom workflow from scratch can be found [here](https://cgat-core.readthedocs.io/en/latest/defining_workflow/Tutorial.html).

Installation
------------

The following sections describe how to install the cgat-showcase pipeline.

We recommend installing using conda and the steps are described below::

   `conda install -c cgat cgatshowcase`

Alternatively, the pipeline can be installed using pip::

   `pip install cgatshowcase`

However, you will require certain software to run the pipeline. More detail on installation can be found on the [Installation](https://cgat-showcase.readthedocs.io/en/latest/getting_started/Installation.html) documentation.
