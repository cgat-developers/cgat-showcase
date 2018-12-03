# cgat-showcase

cgat-showcase is a showcase example pipeline constructed using our [cgat-core](https://github.com/cgat-developers/cgat-core) workflow management system.

Within this repository we have an example of one [pipeline](https://github.com/cgat-developers/cgat-showcase/blob/master/cgatshowcase/pipeline_transdiffexprs.py) that performs pseudoalignment
with kallisto and differential expression using deseq2. It can be ran locally or distributed across a cluster.

Documentation on how to run this pipeline can be found [here](https://cgat-showcase.readthedocs.io/en/latest/) and documentation on how
to build a workflow from scratch can be found [here](https://cgat-core.readthedocs.io/en/latest/defining_workflow/Tutorial.html).

Installation
------------

The following sections describe how to install the cgat-showcase pipeline.

We recommend installing using conda and the steps are described below::

   `conda install -c cgat cgatshowcase`

Alternatively, the pipeline can be installed using pip::

   `pip install cgatshowcase`

However, you will require certain software to run the pipeline. More detail on installation can be found on the [Installation](https://cgat-showcase.readthedocs.io/en/latest/getting_started/Installation.html) documentation.
