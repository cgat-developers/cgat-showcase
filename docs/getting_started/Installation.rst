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

Once you have conda installed then to install cgat-showcase please do the following:

   conda install -c cgat cgatshowcase

This should install all of the dependancies and you should be read to proceded to the Writing workflows and Tutorial section.


Manual installation
-------------------

To obtain the latest code, check it out from the public git repository and activate it::

   git clone https://github.com/cgat-developers/cgat-core.git
   cd cgat-core
   python setup.py develop

Once checked-out, you can get the latest changes via pulling::

   git pull 


Installing additonal software
-----------------------------

When building your own workflows we recomend using conda to install software into your environment where possible.

This can easily be performed by::

   conda search <package>
   conda install <package>
