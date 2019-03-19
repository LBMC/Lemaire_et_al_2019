GC- AT- U1- U2-exons figures's documentation
============================================


The following documentation contains a brief description of what scripts inside the  `GC_rich_AT_rich_exon_list/src` folder does. For more details, read the *Materials and Methods* in the article of **Lemaire *et al***.


.. note::

  All the following scripts use the SED and/or the fasterDB database. To see how to create those databases, you can read the documentation in the `Database_creator/doc` folder. Those databases must be present into the ``data`` folder.


Root Folder
-----------

The script ``figure_creator`` will create a venn diagram showing the exons regulated by splicing factor regulating AT-rich exons (with long introns)(AT-exons) and the exons regulated by splicing factor regulating GC-rich exons (with short introns) (GC-exons).
It also create files containing the GC and AT-exons. Those files contains 2 columns : the first one correspond to the id of a FasterDB [Mal14]_ gene and the second one corresponds to the id of a FasterDB exons.

The script ``exon_set_input_creator`` create files containing the exons activated by DDX5_17, SNRP70, SNRPC, UAF1, U2AF2, SF3A3 and SF3B4. Those files contains 2 columns : the first one correspond to the id of a FasterDB gene and the second one corresponds to the id of a FasterDB exons.



``boxplot_GC_content_and_flanking_intron_size`` folder
-------------------------------------------------------

The scripts in the folder `boxplot_GC_content_and_flanking_intron_size` allow to create boxplots representing the GC content and the smallest flanking intron size in AT- and GC-exons.
They also allow to create  boxplots representing the GC content, the median intron size and the size of the genes hosting only GC- or AT-exons (not both groups).



``metaexon`` folder
--------------------

The scripts into that folder allow to create the frequency maps (for every nucleotides A, T, G, C, S=GC, W=AT) on exon-intron junctions (containing the 3' splicing site or the 5' splicing site) of different group of exons.
The exon-intron junctions corresponds to sequences with100 nucleotides within the intron and 50 nucleotides within the exons.



``minimum_free_energy`` folder
-------------------------------

The scripts into that project allow to create the violin plots representing minimum free energy of different group of exons (AT- GC- U1 and U2-exons).
The computation of the minimum free energy uses the RNAfold script of the ViennaRNA package (version 2.4.1) [Lor11]_. The RNAfold script must be present in the folder `data/RNAfold/`.



``make_control_files_bp_ppt`` folder
-------------------------------------

The scripts into this folder allow to create several kind of figures:

  * The barplots representing the proportion of exons having less or more than two branch points in the 100 nucleotides upstream two different groups of exons (AT/GC exons or U1/U2 dependent exons).
  * The weblogos representing the nucleotide sequences surrounding the predicted branch points in the 25 nucleotides upstream two different groups of exons (AT/GC exons or U1/U2 dependent exons).
  * The boxplots representing the GC content of the weblogos (described above) of the two groups of exons conciderated.
  * The boxplots representing the number of TNA in the 50 nucleotides upstream two groups of exons (AT/GC exons or U1/U2 dependent exons)

Branch points prediction is perform with SVM-BP finder [CHS10]_. To run the script you must change the path to ``svm_bpfinder.py`` in the line 19 of `make_control_files_bp_ppt/function_bp.py`.



``stretch_calculator`` folder
------------------------------

The scripts into this folder allow to create the boxplots representing the number of T-rich low complexity sequence in the 50 nucleotides upstream two groups of exons (AT/GC exons or U1/U2 dependent exons)


``GC_AT_group_regulated_U1_U2`` folder
---------------------------------------

This folder allow the creation of a scatter plot representing the V-Values of each spliceosome associated factors (see method of the paper of **Lemaire *et al*** to better understand what V-value is).



Prerequisites
--------------

Those scripts were written in **Ubuntu 16.04** with the version **3.5 of python**. They should work on other unix system with a bash shell.

Steps to launch those scripts:

  * You may have to install python 3, if you do click `here <https://www.python.org/downloads/release/python-356/>`_
  * Install the required dependencies by running ``sudo pip3 install -r requierements.txt``. The file ``requierements.txt`` is the folder ``Figure_ESA`` and contains the name and the version of every required module.
  * You must copy (or create a shortcut to) the sed and the fasterDB databases into the ``data/`` folder of ``GC_rich_AT_rich_exon_list``.
  * You must create a folder named ``RNAfold`` into the ``data`` folder and put the RNAfold script within it.
  * You must download the ViennaRNA package (2.4.1) [Lor11]_ and put the RNAduplex script in you PATH enviroment variable.
  * You must launch the following commands :

.. code:: bash

  python3 src/make_control_files_bp_ppt/control_bp_ppt.py # create the control file into the ``make_control_files_bp_ppt/`` folder. This may take a long time
  python3 src/metaexon_figure/control_dictionnary.py # create the control file into the ``metaexon_figure`` folder.
  python3 src/minimum_free_energy/control_mfe.py  # create a control mfe file.
  python3 src/stretch_calculator/stretch_calculator.py # control figures

You can download the `ViennaRNA package <https://www.tbi.univie.ac.at/RNA/>`_  and the `SVM BP Finder program <https://bitbucket.org/regulatorygenomicsupf/svm-bpfinder/downloads/>`_


To create the figures you can run the following commands (in the folder ``GC_rich_AT_rich_exon_list``):

.. code:: bash

  # creation of the boxplot GC content/intron size for AT- GC exons/genes
  python3 src/boxplot_GC_content_and_flanking_intron_size/launcher.py

  # Creation of the v-value figure
  python3 src/GC_AT_group_regulated_U1_U2/spliceosome_regulation_enrichment.py
  python3 src/GC_AT_group_regulated_U1_U2/barplot_pvalue_maker.py

  # Creation of the branch points figures
  python3 src/make_control_files_bp_ppt/bp_ppt_figure_creator.py


  # creation of the files containing GC and AT exons
  python3 src/figure_creator.py

  # creation of the files containing the exons activayed by some splicerosome activated factors.
  python3 src/exon_set_input_creator.py

  # Creation of the nucleotide frequency maps
  python3 src/metaexon_figure/launcher_metaexon.py --files result/GC_rich_exons,result/AT_rich_exons --name_files 'GC_pure,AT_pure' --nt C,S,A,T,G,W,Y,R --name_fig GC-ATgroup-legend --exon_type CCE --color '#0000FF,#00aa00' --legend True
  python3 src/metaexon_figure/launcher_metaexon.py --files result/inputs_union/input_SNRPC-union.txt,result/inputs_union/input_SNRNP70-union.txt,result/inputs_union/input_DDX5_DDX17-union.txt,result/inputs_union/input_U2AF1-union.txt,result/inputs_union/input_U2AF2-union.txt,result/inputs_union/input_SF1-union.txt,result/inputs_union/input_SF3A3-union.txt,result/inputs_union/input_SF3B4-union.txt --name_files 'SNRPC,SNRNP70,DDX5_17,U2AF1,U2AF2,SF1,SF3A3,SF3B4' --nt T,S,A,G,C,W,Y,R --name_fig spliceosome_group-legend --exon_type CCE --color 'cyan,navy,#AA00FF,#006400,olive,#55FF55,#D8EF48,#8FBC8F' --legend True

  # Creation of the minimum free energy figures
  python3 src/minimum_free_energy/mfe_figure_creator.py

  # Creation of the T-rich low complexity sequences figures
  python3 src/stretch_calculator/stretch_calculator.py


.. rubric:: References

.. [Lor11] Lorenz, R. et al. ViennaRNA Package 2.0. Algorithms for Molecular Biology 6, (2011).
.. [CHS10] Corvelo, A., Hallegger, M., Smith, C. W. & Eyras, E. Genome-wide association between branchpoint properties and alternative splicing. PLoS Comput Biol 6, e1001016, doi:10.1371/journal.pcbi.1001016 (2010).
.. [Mal14] Mallinjoud, P. et al. Endothelial, epithelial, and fibroblast cells exhibit specific splicing programs independently of their tissue of origin. Genome Res 24, 511-521, doi:10.1101/gr.162933.113 (2014)
