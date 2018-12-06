.. _manual-main:

============================
CGAT-showcase documentation!
============================

.. image:: https://readthedocs.org/projects/cgat-showcase/badge/?version=latest
    :target: http://cgat-showcase.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://img.shields.io/twitter/follow/CGAT_Oxford.svg?style=social&logo=twitter&label=Follow
    :target: https://twitter.com/cgat_oxford?lang=en
    :alt: Twitter Followers

.. image:: https://img.shields.io/twitter/url/http/shields.io.svg?style=social&logo=twitter
    :target: https://twitter.com/cgat_oxford?lang=en
    :alt: Twitter URL

CGAT-showcase contains a showcase example of how to build workflows suing our cgatcore  workflow management system. 
cgatcore allows you to quickly build both simple and complex reproducible pipelines completely within python. 

CGAT-core is  a set of libraries and helper functions used to enable researchers
to design and build computational workflows for the analysis of large-scale data-analysis. Access to code is
through `github <https://github.com/cgat-developers/cgat-showcase>`_ and documentation can be accessed `here <https://cgat-core.readthedocs.io/en/latest/>`_.

Used in combination with `CGAT-apps <https://github.com/cgat-developers/cgat-apps>`_, we have deomonstrated the functionality of our
flexible implementation using a set of well documented, easy to install and easy to use workflows, which can be found in our `cgat-flow <https://github.com/cgat-developers/cgat-flow>`_ repository.

However, in this repository we have developed a toy example showing the functionality of our code. For further information on the advanced functionality of our
code please refer to the `cgat-flow docuemntation <https://www.cgat.org/downloads/public/cgatpipelines/documentation/>`_. 

CGAT-core is open-sourced, powerful and user-friendly, and has been continually developed
as a Next Generation Sequencing (NGS) workflow management system over the past 10 years.


.. _manual-quick_example:

--------
Citation
--------

To be added....

.. _manual-support:

-------
Support
-------

- Please refer to our :ref:`FAQ` section
- In case of questions, please add these to `stack overflow <https://stackoverflow.com/search?q=cgat>`_
- For bugs and issues, please raise an issue on `github <https://github.com/cgat-developers/cgat-showcase>`_



.. toctree::
   :caption: Getting started
   :name: getting_started
   :maxdepth: 1
   :hidden:

   getting_started/Installation.rst
   getting_started/Writing_workflow.rst
   getting_started/Tutorial.rst


.. toctree::
   :caption: Build a workflow
   :name: build_workflow
   :maxdepth: 1
   :hidden:

   build_workflow/writing_workflow.rst
   build_workflow/packaging.rst
   build_workflow/conda_build.rst

.. toctree::
   :caption: Project Info
   :name: project-info
   :maxdepth: 1
   :hidden:

   project_info/Contributing.rst
   project_info/FAQ.rst
   project_info/Licence.rst
