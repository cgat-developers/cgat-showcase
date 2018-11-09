.. _getting_started-Tutorial:


=============================
Running a pipeline - Tutorial
=============================

Before beginning this tutorial plase make sure you have cgat-showcase installed correctly, please see here (see :ref:`getting_started-Installation`) for instructions.

As a tutorial example we have written a showcase example of how you can build robust and scalable pipelines
written in python. For this toy example we have written a workflow that performs pseudocounting using kallisto and differential expression
analysis using deseq. Two reports are then generated, one using `MultiQC <https://multiqc.info/>` and another written in
`Rmarkdown <https://rmarkdown.rstudio.com/>`. However, because of the flexible nature of our workflow management system,
any reporting strategy can be implimented.

Tutorial
--------

**1.** First download the tutorial data::

   mkdir showcase
   cd showcase
   wget 
   tar -zxvf 

**2.** Next we will generate a configuration yml file so the pipeline output can be modified::

   cgatshowcase transdiffexpres config

or you can alternatively call the workflow file directly::

   python /path/to/file/pipeline_transdiffexpres.py config

This will generate a **pipeline.yml** file containing the configuration parameters than can be used to modify
the output of the pipleine. However, for this tutorial you do not need to modify the parameters to run the 
pipeline. In the :ref:`modify_config` section below I have detailed how you can modify the config file to
change the output of the pipeline.

**3.** Next we will run the pipleine::

   cgatshowcase transdiffexpres make full -v5 --local

This ``--local`` will run the pipeline locally if you do not have access to a cluster. Alternatively if you have a
cluster remove the ``--local`` option and the pipleine will distribute your jobs accross the cluster.

.. note::

   There are many commandline options available to run the pipeline. To see available options please run :code:`cgatshowcase --help`.

**4.** Generate a report

The final step is to generate a report to display the output of the pipeline. We have a preference for using MultiQC
for generate bioinformatics tools (such as mappers and pseudoaligners) and Rmarkdown for generating custom reports.
In order to generate these run the command::

    cgatshowcase transdiffexprs make build_report -v 5 --local

This will generate a MultiQC report in the folder MultiQC_report.dir/ and an Rmarkdown report in R_report.dir/. 



This completes the tutorial for running the transdiffexprs pipeline for cgat-showcase, hope you find it as useful as
we do for writing workflows within python. 

.. _modify_config

Modify the configuration file
-----------------------------

Having generated a default pipeline.yml file, you may require the output of the pipeline to be
modified. Foe example, if you require a different genome or a different fdr value for significance. These can
be changed modifying the yml file and adding user specific infomation.

TO BE ADDED AS I HAVENT FINISHED THE PIPELINE
