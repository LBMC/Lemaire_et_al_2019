Description of ``correlation_iupac_gene_and_intron_size`` script
=================================================================

The goal of this script is :
  1. To create correlation graphics between the relative median **gene_size/median_intron_size** and the relative gene median frequency of a iupac nucleotide (A, T, G, C, R, W, Y, R) for every up, down, up+down regulated exons in every splicing lore project.
  2. To create correlation graphics between the relative median exon frequency of a iupac nucleotide (A, T, G, C, R, W, Y, R) of up and down regulated exons in every splicing lore project. (a dot corresponds to a spliing lore project)
  3. To create correlation graphics between the relative median frequency of a nucleotide (A, T, G, C, R, W, Y, R) in every group of exon (up or down) and their associated gene. (a dot corresponds to a splicing lore project)
  4. To create the correlation graphics between the relative median frequency of nucleotide (A, T, G, C, R, W, Y, R) in every group of exons (up or down) and their associated proximal intronic regions (25 nucleotide befores and after the group of exons of interest) (for every splicing factor)
  5. To create the correlation graphics between the relative median frequency of nucleotide (A, T, G, C, R, W, Y, R) in every group of exons (up or dowb) and their associated splicing force of the donor (5') or acceptor (3') splicing site (a dot correponds to a splicing lore project)
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

  # creation of correlation graphics between the median intron_size/gene_size and the frequency of a (iupac) nucleotides for every slicing lore projects.
  python3 src/correlation_iupac_gene_and_intron_size.py --type "gene_size_vs_gene_iupac"
  # creation of correlation graphics between the median exons frequency of a (iupac) nucleotides between up and down exons of every slicing lore project.
  python3 src/correlation_iupac_gene_and_intron_size.py --type "iupac_up_vs_iupac_down"
  # creation of correlation graphics between the median frequency of a nucleotide in a (up or down) group of exons and their associated genes (for every splicing lore project)
  python3 src/correlation_iupac_gene_and_intron_size.py --type "iupac_gene_vs_iupac_exon"
  # creation of correlation graphics between the median frequency of a nucleotide in a (up or down) group of exons and their associated proximal intronic sequence (for every splicing lore project)
  python3 src/correlation_iupac_gene_and_intron_size.py --type "iupac_exon_vs_iupac_intron_proxi"
  # creation of correlation graphics between the median frequency of a nucleotide in a (up or down) group of exons and their associated spling force (donor -5' and acceptor 3') (for every splicing lore project)
  python3 src/correlation_iupac_gene_and_intron_size.py --type "force_vs_iupac_exon"
