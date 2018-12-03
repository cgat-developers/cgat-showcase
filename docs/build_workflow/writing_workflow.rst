.. _writing_workflow-tutorial:


=============================
Writing a pipeline - Tutorial
=============================

Before attempting this tutorial please make sure you have cgat-showcase installed correctly. please see here (see :ref:`getting_started-Installation`) for instructions.

This tutorial is intended as a walkthrough to show you the processes inviolved in writing the cgat-showcase transdiffexprs pipeline. It wont go through all of the functions in the pipeline but it should give you an idea of how to construct a pipeline using out cgat-core python library. Before attempting this tutorial, it may be best to visit cgat-core `documentation <https://cgat-core.readthedocs.io/en/latest/defining_workflow/Tutorial.html>`_  for an idea of how to construct a very simple pipeline.

Tutorial
--------
**1.** The problem and the solution

The first step in any pipeline is to think about the problem that you are trying solve. 
For example: 

**The probelm:** I need to perform differential expression analysis quickly and consistently across all my samples. I also need the solution to scale and be ran on my mac or districuted acrposs a cluster.
**The solution:** Build a cgat-core workflow that wraps pseduocounting tools differential expression analysis tools. 

Then the next step is to think about the steps needed to reach the solution.

For example:

   1.  
   2. 
