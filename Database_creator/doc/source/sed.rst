SedDB Creation's documentation
==============================

Description
-----------

This program contains the following scripts:
  * ``exon_class.py`` : This script aims to calculates lots of information for every exons from **FasterDB-Lite database**
  * ``exon_information_retrieverr.py``: This script will create and empty **Sed database** and then fill it (table **sed**) thanks to ``exon_class.py``. Then it will create 2 additional tables : **rnaseq_project** and  **ase_events** thanks to splicing lore database.


At the end the **Sed (Simple Exons description) database** created have the following relational schema :

.. figure:: images/sed_schema.png
  :align: center

  **Relational schema of the Sed database**

Here is described every field in the table sed:

+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|          **Field**                  |                                                                         **Description**                                                                              |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|           gene_symbol               | HGNC symbol of the gene                                                                                                                                              |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|           gene_id                   | FasterDB identifiant of a gene                                                                                                                                       |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|           exon_pos                  | Position of the exon within the gene (given in FasterDB)                                                                                                             |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|           exon_type                 | The type the exon, it can be CCE (Constitutive Coding Exon), ACE (Alternative Coding Exon), FCE (First Coding Exon), LCE (Last Coding Exon)                          |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
+          gene_size                  | The size (in nucleotide) of the gene that contains the exons                                                                                                         |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|          nb_intron_gene             |  The number of intron within the gene                                                                                                                                |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|          median_intron_size         | The median of the size of every intron in the gene that contains the exon of interest                                                                                |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|          iupac_gene                 | The frequency (%) of nucleotides A, C, G, T, S, W, R, Y, K, M respectively within the gene that contains the exon of interest                                        |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|          dnt_gene                   | The frequency (%) of di-nucleotides  AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT respectively within the gene that contains the exon of interest  |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|     upstream_exon_size              | The nucleotide size of the upstream exon from the one of interest                                                                                                    |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|         exon_size                   | The nucleotide size of the exon of interest                                                                                                                          |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|       downstream_exon_size          | The nucleotide size of the downstream exon from the one of interest                                                                                                  |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|  force_acceptor_upstream_exon       | The force of the splice acceptor site of the upstream exon                                                                                                           |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|       force_acceptor                | The force of the splice acceptor site of the interest exon                                                                                                           |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|  force_acceptor_downstream_exon     | The force of the splice acceptor site of the downstream exon                                                                                                         |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|  force_donor_upstream_exon          | The force of the splice donor site of the upstream exon                                                                                                              |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|       force_donor                   | The force of the splice donor site of the interest exon                                                                                                              |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|  force_donor_downstream_exon        | The force of the splice donor site of the downstream exon                                                                                                            |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|          iupac_exon                 | The frequency (%) of nucleotides A, C, G, T, S, W, R, Y, K, M respectively within the exon of interest                                                               |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|          dnt_exon                   | The frequency (%) of dinucleotides  AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT respectively within the exon of interest                          |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|  iupac_upstream_intron              | The frequency (%) of nucleotides A, C, G, T, S, W, R, Y, K, M respectively within the region [-100,-26] of the upstream intron                                       |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| dnt_upstream_intron                 | The frequency (%) of dinucleotides AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT respectively within the region [-100,-26] of the upstream intron   |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|  iupac_upstream_intron_ppt_area     | The frequency (%) of nucleotides A, C, G, T, S, W, R, Y, K, M + seq size respectively within the region [-75,-25] of the upstream intron (0 first nucleotide of exon)|
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|  iupac_upstream_intron_adjacent1    | The frequency (%) of nucleotides A, C, G, T, S, W, R, Y, K, M + seq size respectively within the region [-25,-0] of the upstream intron (0 first nucleotide of exon) |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|  iupac_upstream_intron_adjacent2    | The frequency (%) of nucleotides A, C, G, T, S, W, R, Y, K, M + seq size respectively within the region [-50,-0] of the upstream intron (0 first nucleotide of exon) |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|  iupac_upstream_intron_proxi        | The frequency (%) of nucleotides A, C, G, T, S, W, R, Y, K, M + seq size respectively within the region [-100,-0] of the upstream intron (0 first nucleotide of exon)|
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|  dnt_upstream_intron_proxi          | The frequency (%) of dinucleotides AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT + seq size respectively within region[-100,-0] of upstream intron  |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|  iupac_downstream_intron_adjacent1  | The frequency (%) of nucleotides A, C, G, T, S, W, R, Y, K, M + seq size respectively within the region [1,25] of the downstream intron (0 last nucleotide of exon)  |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|  iupac_downstream_intron_adjacent2  | The frequency (%) of nucleotides A, C, G, T, S, W, R, Y, K, M + seq size respectively within the region [1,50] of the downstream intron (0 last nucleotide of exon)  |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|  iupac_downstream_intron_proxi      | The frequency (%) of nucleotides A, C, G, T, S, W, R, Y, K, M + seq size respectively within the region [1,100] of the downstream intron (0 last nucleotide of exon) |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|  dnt_downstream_intron_proxi        | The frequency (%) of dinucleotides AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT + seq size respectively within region [1,100] of downstream intron |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|  iupac_upstream_intron              | The frequency (%) of nucleotides A, C, G, T, S, W, R, Y, K, M respectively within the region [26;100] of the downstream intron                                       |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|  dnt_upstream_intron                | The frequency (%) of dinucleotides AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT respectively within the region [26;100] of the downstream intron   |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|  iupac_intron_exon                  | The frequency (%) of nucleotides A, C, G, T, S, W, R, Y, K, M and sequence size, respectively within the region [-100;+50] of the upstream intron                    |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|  dnt_intron_exon                    | The frequency (%) of dinucleotides AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT seq size respectively within region [-100;50] of upstream intron   |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|  iupac_exon_intron                  | The frequency (%) of nucleotides A, C, G, T, S, W, R, Y, K, M and sequence size, respectively within the region [-100;+50] of the downstream intron                  |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|  dnt_exon_intron                    | The frequency (%) of dinucleotides AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT seq size respectively within region [-100;50] of downstream intron |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|  iupac_exon_env                     | The frequency (%) of nucleotides A, C, G, T, S, W, R, Y, K, M and sequence size, respectively within the region [-100;exon;+100] of the downstream intron            |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|  dnt_exon_env                       | The frequency (%) of dinucleotides AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT seq size respectively within region [-100;exon;100] of down-intron |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| relative_donor_upstream             | ((force_donor - force_donor_upstream_exon) / force_donor_upstream_exon) * 100                                                                                        |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| relative_donor_downstream           | ((force_donor - force_donor_downstream_exon) / force_donor_downstream_exon) * 100                                                                                    |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| relative_acceptor_upstream          | ((force_acceptor - force_acceptor_upstream_exon) / force_acceptor_upstream_exon) * 100                                                                               |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| relative_acceptor_downstream        | ((force_acceptor - force_acceptor_downstream_exon) / force_acceptor_downstream_exon) * 100                                                                           |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+

.. note::

  The columns labeled *_proxi* and *_exon_intron* or *_intron_exon* have a suplementary information corresponding to the size of the sequence used to compute the iupac or dnt frequencies

Description of the **rnaseq_projects** table:

.. note::

	A project here is an experiment where the trancriptome of a cell line depleted for a splicing factor is compared to the trancriptiome wild-type of the same cell line


+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|          **Field**                  |                                                                         **Description**                                                                              |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|            id                       | The id of the project : an unique identifier for a project on a particular cell line for a particular splicing factor                                                |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|           project_name              | The name of the project : like this : SFname_DBID_CellLine where SFname is a name of a plicing factor, DBid is an id project like GSE00000 and cell line a cell line |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|           source_db                 | The database where the project was downloaded (GEO, DRAsearch, EBI, HOME, ENCODE)                                                                                    |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|           db_id_project             | The id of the project in the database where the project is located                                                                                                   |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|           sf_name                   | the name of the splicing factor studied in the project                                                                                                               |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|           cl_name                   | The name of the cell line used in the project                                                                                                                        |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+


Description of the **ase_event** table:

This table describe the exons that are differentially skipped in each project defined in **rnaseq_projects**. Those value were obtain using farline.

+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|            **Field**                |                                                                         **Description**                                                                              |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|           id                        |  An unique identifier of a splicing event in a particular project on a particular cell line/splicing factor                                                          |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|           id_project                |  Foreing key of the field id in rnaseq_projects table                                                                                                                |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|           id_gene                   | The gene id of the gene that contains the exon differentially splicing in the project identified by id_project                                                       |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|           gene_symbol               | The HGNC symbol of the gene that contains the exon differentially splicing in the project identified by id_project                                                   |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|           exon_skipped              | The exon position skipped on the gene identified by gene_id                                                                                                          |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|           chromosome                | The chromosome where the exon differentially spliced is located                                                                                                      |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|           start                     | Chromosome coordinates where  the exon differentially spliced begins                                                                                                 |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|           stop                      | Chromosome coordinates where  the exon differentially spliced ends                                                                                                   |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|           exons_flanquant           | The position of the surrounding exons of the one differentially spliced in the gene                                                                                  |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|           deltapsi                  | The differential inclusion of the exon differentially  spliced (negative value: exon less included in the absence of a splicing factor)                              |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|           pvalue                    | The pvalue of the splicing events                                                                                                                                    |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|           pvalue_glm_cor            | The pvalue corrected (if many biological replicate are available)                                                                                                    |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+


.. note::

  This schema induce a lots of redundancy in the database. Indeed, we keep for each exon, data about its gene, so, for a gene we have the same data repeated as many times as the number of exons within the gene.
  The choice of keeping redundancy was made to improve the speed of programs that will use *Sed database*; Indeed, for one exon we have all the data we need. Despite this redundancy, the size of the database is reasonable (a little more than 100 Mo).


.. note::

  The frequencies of nucleotides and dinucleotides for a feature is only reported if 95% of the feature sequence is well defined (not unidentified N nucleotides). The frequency givene for a feature doesn't take into account the undefined nucleoitdes.

.. note::
  If a feature show a lenght below 1 nucleotide, this lenght is reported as a *NULL* value in the sed database

Issue
-----

.. warning::

        There are 4 exons in fasterDB having a length below 0 nucleotide. Those exons are present in SED database too.


Prerequisite
------------

This program uses `python <https://www.python.org>`_ version ``3.5`` and this following dependencies:
  * `numpy v1.14.0 <https://docs.scipy.org/doc/numpy-1.14.0/user/quickstart.html>`_
  * `sqlite3 v2.6.0 <https://docs.python.org/3.5/library/sqlite3.html>`_ : To create *Sed* database
  * `sys v3.5.2 <https://docs.python.org/3.5/library/sys.html>`_
  * `re v2.2.1 <https://docs.python.org/3.5/library/re.html>`_
  * `pymysql v0.8.0 <https://pymysql.readthedocs.io/en/latest/>`_

Exectuted commands to create the *Sed* database
---------------------------------------------------------

.. code-block:: bash

	python3 src/exon_information_retriever.py
