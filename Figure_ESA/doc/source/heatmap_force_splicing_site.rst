relative heatmap_creator's Documentation
==========================================

Description
------------

This script aims, to cluster exons regulated by every **splicing lore project** or **splicing factor** by :
  * **relative median toward control** for ``force_donor`` and ``force_acceptor`` characteristics
  * **median** for ``relative_donor_upstream``, ``relative_donor_downstream``, ``relative_acceptor_upstream`` and ``relative_acceptor_downstream``

For a the characteristics `force_donor`` and ``force_acceptor`` denoted :math:`X` here, the relative median is calculated as follow:

.. math::

  relative\_median = \frac{median(X_{obs}) - median(X_{ctrl})}{median(X_{ctrl})}

Where:
  * :math:`X_{obs} = \{ {X_{{obs}_1}, ..., X_{{obs}_n}} \}`  where :math:`X_{{obs}_i}` is the value of :math:`X` for the exon :math:`i` of the interest sets of exons.
  * :math:`X_{ctrl} = \{{X_{{ctrl}_1}, ..., X_{{ctrl}_m}}\}`  where :math:`X_{{ctrl}_j}` is the value of :math:`X` for the exon :math:`j` of the interest sets of exons.


The median of ``relative_donor_upstream``, ``relative_donor_downstream`` is calculated as follow:

  relative\_donor\_upstream = median({\frac{F_{n}_{i} - F_{n-1}_{1}}{F_{n-1}_{1}} * 100, ..., \frac{F_{n}_{k} - F_{n-1}_{k}}{F_{n-1}_{k}})
  relative\_donor\_downstream = median({\frac{F_{n}_{i} - F_{n+1}_{1}}{F_{n+1}_{1}} * 100, ..., \frac{F_{n}_{k} - F_{n+1}_{k}}{F_{n+1}_{j}})

Where:
  * :math:`F_{n}_{i}` is the force of the donor splicing site for the :math:`i^{th}` exons regulated by a splicing factor located in the position :math:`n` in a gene
  * :math:`F_{n+1}_{i}` is the force of the donor splicing site for the :math:`i^{th}` exons regulated by a splicing factor located in the position :math:`n+1` in a gene (downstream exon than the one of interest)
  * :math:`F_{n-1}_{i}` is the force of the donor splicing site for the :math:`i^{th}` exons regulated by a splicing factor located in the position :math:`n-1` in a gene (upstream exon than the one of interest)
  * :math:`k` is the number of exons regulated by a splicing

The same formula was applied for the acceptor splicing forces.

For each project it will produce this kind of graphics :

.. Figure:: images/rhp.png
  :align: center

.. note::

  Note that a white color doest not mean the same thing for every column in the heatmap.
    * For `force_donor`` and ``force_acceptor`` a white color means that the splicing forces are equal to those of the control exons.
    * For `relative_donor_upstream``, ``relative_donor_downstream``, ``relative_acceptor_upstream`` and ``relative_acceptor_downstream`` a white color means that the splicing forces are equal to those of the upstream/doownstream exons




Each dendrogram was computed using an **Enclidean distance matrix** and a **complete linkage clustering**.
For more information see :

  * `Euclidean distance matrix <https://en.wikipedia.org/wiki/Euclidean_distance_matrix>`_
  * `Complete-linkage clustering <https://en.wikipedia.org/wiki/Complete-linkage_clustering>`_


Prerequisites:
--------------

This program uses `python <https://www.python.org>`_ version ``3.5`` and this following dependencies:
  * ``figure_producer`` : This script should be present in ``src`` folder of the project directory
  * ``exon_control_handler`` : This script should be present in ``src`` folder of the project directory
  * `plotly v2.7.0 <https://plot.ly/python/>`_
  * `os <https://docs.python.org/3.5/library/os.html>`_
  * `random <https://docs.python.org/3.5/library/random.html>`_
  * `numpy <http://www.numpy.org/>`_


Command Line executed to create the graphics
============================================


.. code:: bash

  # figure displayed for exons regulated in every projects
  python3 src/relative_heatmap_creator.py
  # figure displayed for every exons regulated by a splicing factor.
  python3 src/relative_heatmap_creator.py union
