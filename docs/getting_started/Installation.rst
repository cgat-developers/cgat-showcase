.. _getting_started-Installation:


============
Installation
============


The following sections describe how to install the cgatshowcase framework. Since there are very minimal dependancies
the installation of the code should not take too long. The main dependancies are cgat-core, cgat-apps, kallisto, deseq2,
biomart, tximport and genomeinfodb.


Installation using conda
------------------------

The most convinient way to use our code is to install the software using conda. First you will need to install
conda using `miniconda <https://conda.io/miniconda.html>`_ or `anaconda <https://www.anaconda.com/download/#linux>`_ and following the instructions `here <https://conda.io/docs/user-guide/install/linux.html>`_. (we recomend miniconda rather than anaconda) 

Once you have conda installed then to install cgat-showcase please do the following::

   conda install -c cgat cgatshowcase

This should install all of the dependancies and you should be ready to proceded to the Writing workflows and Tutorial section.
However, conda is currently having issues with speed of installation. This is related to the solver and is something that is known and the conda developers are workin on new fixes for this.
Please see `conda issue <https://github.com/conda/conda/issues/7239>`_ for more information.

If you find that the package installation just sits at `Solving environent`, just be patient this will take a while but will install. Alternmatively install using the anaconda cloud environment ppackage, using the instructions below.

Installtion using conda environment
----------------------------------

The package distribution is quite large at the moment and as a consequence of problems with the conda solver it takes quite a long time to install the conda package.

As a temporary work around we have included the conda environemnt used to run the cgatshowcase on the anaconda cloud. Please follow the instructions below::

    conda env create cgat/cgatshowcase-env

You will then need to clone the cgatshowcase resository and run setup as follows::

    conda activate cgatshowcase-env
    git clone https://github.com/cgat-developers/cgat-showcase.git
    cd cgat-showcase
    python setup.py develop
    cgatshowcase --help

**Or** use one of our environment.yml files (see conda-envs/environment-mac.yml or conda-envs/environment-linux.yml). Then install a new
environment by::

    conda env create -f [environment-mac.yml/environment-linux.yml]


Pip installation
----------------

Alternatively, cgatshowcase can also be installed using pypi::

   pip install cgatshowcase

However, the dependancies will have to be installed seperately.

Manual installation
-------------------

To obtain the latest code, check it out from the public git repository and activate it::

   git clone https://github.com/cgat-developers/cgat-showcase.git
   cd cgat-showcase
   python setup.py develop

Once checked-out, you can get the latest changes via pulling::

   git pull 


Installing additonal software
-----------------------------

When building your own workflows we recomend using conda to install software into your environment where possible.

This can easily be performed by::

   conda search <package>
   conda install <package>
