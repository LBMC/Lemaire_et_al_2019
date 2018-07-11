Description of ``correlation_iupac_gene_and_intron_size`` script
=================================================================

The goal of this script is to create correlation graphics between the relative median **gene_size/median_intron_size** and the relative median frequency of a iupac nucleotide (A, T, G, C, R, W, Y, R) for every up, down, up+down regulated exons in every splicing lore project.


Example
-------

Let's say we want to create the correlation graphics of gene_size between the A nucleotide fréquency of down-regulated exons.
To do that we :

  * calculate the median gene_size of every splicing lore projects (i.e. the median of gene_size for every down-regulated exons in each particular project)
  * calculate the median frequnecy of A nucleotide of every splicing lore projects (i.e. the median of A nucleotide frequency for every down-regulated exons in each particular project)
  * For each median previously calculated we calculate the relative median toward CCE control exons (i.e the median value of gene_size or A nucleotide frequency in CCE (consitutive control exons)
    * A relative median is calculated as follow :math:`relative\_median=\frac{median\_obs - cce\_median}{cce\_median}`
  * We then create figure.


Prerequisites
=============
This program uses `python <https://www.python.org>`_ version ``3.5`` and this following dependencies:

  * ``figure_producer`` : This script should be present in ``src`` folder of the project directory
  * ``exon_control_handler`` : This script should be present in ``src`` folder of the project directory
  * `plotly v2.7.0 <https://plot.ly/python/>`_
  * `os <https://docs.python.org/3.5/library/os.html>`_
  * `numpy <http://www.numpy.org/>`_
  * `scipy <https://www.scipy.org/>`_


Command Line executed to create the graphics
============================================

.. code:: bash

  python3 src/correlation_iupac_gene_and_intron_size.py
