``relative_exon_size``'s Documentation
============================================

Description of ``control_exon_adapter`` script
----------------------------------------------

This script uses the ``control_exon_adapter`` script.

The ``control_exon_adapter`` script allows to adapt the control dictionary file by adding 3 columns to it:

  * ``rel_exon_intron_up``: this column is calculated with by dividing each exon_size by the upstream_intron_size for every control exons. Then the median of those values is computed.
  * ``rel_exon_intron_down``: this column is calculated with by dividing each exon_size by the downstream_intron_size for every control exons. Then the median of those values is computed.
  * ``rel_exon_introns`` : this column is calculated with by dividing each exon_size by the average size of the downstream_intron and the upstream intron for every control exons. Then the median of those values is computed.

Description of ``relative_exon_size`` script
----------------------------------------------

This script aims, to cluster exons regulated by every **splicing lore project** or **splicing factor** by the following caracteristics :

  *  ``rel_exon_intron_up`` (same things for control exons in ``control_exon_adapter`` script but for exons regulated in a given splicing lore project)
  * ``rel_exon_intron_down``  (same things for control exons in ``control_exon_adapter`` script but for exons regulated in a given splicing lore project)
  * ``rel_exon_introns``  (same things for control exons in ``control_exon_adapter`` script but for exons regulated in a given splicing lore project)


For a given characteristic :math:`X` the relative median is calculated as follow:

.. math::

  relative\_median = \frac{median(X_{obs}) - median(X_{ctrl})}{median(X_{ctrl})}

Where:
  * :math:`X_{obs} = \{ {X_{{obs}_1}, ..., X_{{obs}_n}} \}`  where :math:`X_{{obs}_i}` is the value of :math:`X` for the exon :math:`i` of the interest sets of exons.
  * :math:`X_{ctrl} = \{{X_{{ctrl}_1}, ..., X_{{ctrl}_m}}\}`  where :math:`X_{{ctrl}_j}` is the value of :math:`X` for the exon :math:`j` of the interest sets of exons.


Each dendrogram was computed using an **Enclidean distance matrix** and a **complete linkage clustering**.
For more information see :

  * `Euclidean distance matrix <https://en.wikipedia.org/wiki/Euclidean_distance_matrix>`_
  * `Complete-linkage clustering <https://en.wikipedia.org/wiki/Complete-linkage_clustering>`_


Prerequisites:
--------------

This program uses `python <https://www.python.org>`_ version ``3.5`` and this following dependencies:
  * ``figure_producer`` : This script should be present in ``src`` folder of the project directory
  * ``exon_control_adapter`` : This script should be present in ``src`` folder of the project directory
  * `os <https://docs.python.org/3.5/library/os.html>`_
  * `numpy <http://www.numpy.org/>`_
  * `sys <https://docs.python.org/3.5/library/sys.html>`_


Command Line executed to create the graphics
----------------------------------------------


.. code:: bash

  python3 src/relative_exon_size.py
  # whith union exon datasets
  python3 src/relative_exon_size.py union
