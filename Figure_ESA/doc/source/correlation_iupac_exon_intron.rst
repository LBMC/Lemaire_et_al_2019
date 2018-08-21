Description of ``correlation_iupac_exon_intron`` script
=================================================================

The goal of this script is :
  1. To create the correlation graphics between the relative median frequency of nucleotide (A, T, G, C, R, W, Y, R) in every group of exons (up or down) and their associated intronic regions (upstream and downstream introns) (for every splicing factor)

Example
-------

Let's say we want to create the correlation graphics of gene_size between the A nucleotide frequency of down-regulated exons.
To do that we :

  * calculate the median gene_size of every splicing lore projects (i.e. the median of gene_size for every down-regulated exons in each particular project)
  * calculate the median frequnecy of A nucleotide of every splicing lore projects (i.e. the median of A nucleotide frequency for every down-regulated exons in each particular project)
  * For each median previously calculated we calculate the relative median toward CCE control exons (i.e the median value of gene_size or A nucleotide frequency in CCE (consitutive control exons)
    * A relative median is calculated as follow :math:`relative\_median=\frac{median\_obs - cce\_median}{cce\_median}`
  * We then create figure.

.. note::

  This graphics can also be created at the **exon level** : it means that for each splicing factor, every exons thery up- or down regulated will be displayed in the graphics. To display the exon regulated by each splicing factor we've made an union between the exons up and down regulated by a given splicing factor in multiple cell line. The method used to create those **union exons set** are displayed below


Example:

    Let's say we have some sets of exons regulated by SRSF1 in two cell lines as follow:

    **SF:** SRSF1 - **Cell line:** Hela

    +------------+-----------+---------------+
    |  gene      | exon_pos  | regulation    |
    +------------+-----------+---------------+
    | SNRPC      |    3      |     UP        |
    +------------+-----------+---------------+
    | hnRNPK     |    4      |     DOWN      |
    +------------+-----------+---------------+
    | MBNL1      |    9      |     DOWN      |
    +------------+-----------+---------------+


    **SF:** SRSF1 - **Cell line:** 293T

    +------------+-----------+---------------+
    |  gene      | exon_pos  | regulation    |
    +------------+-----------+---------------+
    | SNRPC      |    3      |     DOWN      |
    +------------+-----------+---------------+
    | hnRNPK     |    4      |     DOWN      |
    +------------+-----------+---------------+
    | FUS        |    19     |     UP        |
    +------------+-----------+---------------+

    Then we will have the following exon set in result :

    +------------+-----------+---------------+
    |  gene      | exon_pos  | regulation    |
    +------------+-----------+---------------+
    | hnRNPK     |    4      |     DOWN      |
    +------------+-----------+---------------+
    | MBNL1      |    9      |     DOWN      |
    +------------+-----------+---------------+
    | FUS        |    19     |     UP        |
    +------------+-----------+---------------+

.. note::

        * data about ``SNRPC_3`` is lost because it show different regulation in the cell lines of interest. \
        * ``hnRNPK_4`` exons is only displayed once (we do not repeat exon information). \
        * ``MBNL1_9`` and ``FUS_19`` are both displayed because they are regulated by SRSF1 in at least 1 cell line.

.. note::

  At the **exon level** the graphics show the real value (for a given caracteristics), not relative values are calculated.


Prerequisites
-------------
This program uses `python <https://www.python.org>`_ version ``3.5`` and this following dependencies:

  * ``figure_producer`` : This script should be present in ``src`` folder of the project directory
  * ``exon_control_handler`` : This script should be present in ``src`` folder of the project directory
  * `plotly v2.7.0 <https://plot.ly/python/>`_
  * `os <https://docs.python.org/3.5/library/os.html>`_
  * `numpy <http://www.numpy.org/>`_
  * `scipy <https://www.scipy.org/>`_


Command Line executed to create the graphics
--------------------------------------------

.. code:: bash

  # Project level
  # To create the correlation graphics between the relative median frequency of nucleotide (A, T, G, C, R, W, Y, R) in every group of exons (up or down) and their associated intronic regions (upstream and downstream introns) (for every splicing factor)
  python3 correlation_iupac_exon_intron.py

  # Exon level
  python3 correlation_iupac_exon_intron.py exon_level
