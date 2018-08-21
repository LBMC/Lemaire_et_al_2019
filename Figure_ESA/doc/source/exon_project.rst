Exon projectr's Documentation
==========================================

Description
-----------

This script aims, for each project to show the number of regulated exons contained in the same gene.


For each project it will produce this kind of graphics :

.. Figure:: images/img2.png
  :align: center

Each bar represent the percentage of exons sharing the same gene.
The X labels says that there is 1, 2, 3, 4 or 5 and more (5+) regulated exons belonging to the same gene.


Prerequisites
-------------

This program uses `python <https://www.python.org>`_ version ``3.5`` and this following dependencies:
  * ``figure_producer`` : This script should be present in ``src`` folder of the project directory
  * `plotly v2.7.0 <https://plot.ly/python/>`_
  * `os <https://docs.python.org/3.5/library/os.html>`_


Command Line executed to create the graphics
--------------------------------------------


.. code:: bash

  python3 src/exon_project_analysis.py
