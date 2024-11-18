********
LevSeq
********


 Sequence-function data  provides valuable information about the protein functional landscape, but is rarely obtained during directed evolution campaigns. Here, we present Long-read every variant Sequencing (LevSeq), a pipeline that combines a dual barcoding strategy with nanopore sequencing to rapidly generate sequence-function data for entire protein-coding genes. LevSeq integrates into existing protein engineering workflows and comes with open-source software for data analysis and visualization. The pipeline facilitates data-driven protein engineering by consolidating sequence-function data to inform directed evolution and provide the requisite data for machine learning-guided protein engineering (MLPE). LevSeq enables quality control of mutagenesis libraries prior to screening, which reduces time and resource costs. Simulation studies demonstrate LevSeq’s ability to accurately detect variants under various experimental conditions. Finally, we show LevSeq’s utility in engineering protoglobins for new-to-nature chemistry. Widespread adoption of LevSeq and sharing of the data will enhance our understanding of protein sequence-function landscapes and empower data-driven directed evolution.


Extending LevSeq
================

1. Make a pull request on github - we have made the code extendable for the loss function etc.


Citing LevSeq
=============
LevSeq can be cited as in :ref:`references`, where we also provide citations for the used tools (e.g. matplotlib & seaborn).

.. toctree::
   :caption: Getting started
   :maxdepth: 1

   about
   installing/index


.. toctree::
   :caption: Running LevSeq
   :maxdepth: 1

   examples/Figures

.. toctree::
   :caption: About
   :maxdepth: 1

   faq
   changelog
   references
