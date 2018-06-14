heatmap_creator's Documentation
==========================================

Description
------------

This script aims, to cluster exons regulated by every **splicing lore project** or **splicing factor** by relative medians of characteristics of interest.

For a given characteristic :math:`X` the relative median is calculated as follow:

.. math::

  relative\_median = \frac{median(X_{obs}) - median(X_{ctrl})}{median(X_{ctrl})}

Where:
  * :math:`X_{obs} = \{ {X_{{obs}_1}, ..., X_{{obs}_n}} \}`  where :math:`X_{{obs}_i}` is the value of :math:`X` for the exon :math:`i` of the interest sets of exons.
  * :math:`X_{ctrl} = \{{X_{{ctrl}_1}, ..., X_{{ctrl}_m}}\}`  where :math:`X_{{ctrl}_j}` is the value of :math:`X` for the exon :math:`j` of the interest sets of exons.


For each project it will produce this kind of graphics :

.. Figure:: images/img3.png
  :align: center

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
  python3 src/heatmap_creator.py
  # figure displayed for every exons regulated by a splicing factor.
  python3 src/heatmap_creator.py union
