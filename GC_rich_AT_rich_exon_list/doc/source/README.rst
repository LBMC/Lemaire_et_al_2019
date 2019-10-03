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



``correlation_GC_exon_TAD_and_Exon`` folder
--------------------------------------------

The scripts in the folder ``correlation_GC_exon_TAD_and_Exon`` allows to create figures showing the correlation of different features (such as the GC content) between an exon and the other exons in the TAD (the mean value is display for those exons).
They also allow to build models to see if the exons having particular properties are randomly distributed between TADs or not.



``experimental_branch_point`` folder
--------------------------------------------
The scripts in the folder ``experimental_branch_point`` allow to create barplots figure with experimental branch points as for predicted branch points in ``make_control_files_bp_ppt``. It also allow to compute the intersection between experimental branch points and predicted branch points.

..note:
  The predicted branch points are predicted with svm_bpfinder [CHS10]_.
  The experimental branch points where recovered from previous publications [MCA15]_ [PiBr18]_ [TLS17]_/


``figures_2_3_creator`` folder
-------------------------------

This folder contains scritps that allows the creation of the full figure 2 of the article of Lemaire *et al* and partial parts of figures 1 and 3 for any custom list of exons given to the scripts.
They were used to produce some figures with the AT- GC exons and the exons regulated by SRSF3, SFSF2 and HNRNPC. They were also used to produce those figures with the GA and CT-exons.


``metaexon_figure`` folder
--------------------------

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


``study_microexons`` folder
----------------------------

This folder contains different scripts:
``check_enrichment_in_GC_AT_group.py``  allow to create barplot graphics showing if Irimia microexons [IWE14]_ (or fasterDB small exons with a size < 27 nucleotides) are enriched in GC/AT-like exons when compared to control exons. GC/AT-like exons are defined below:
* AT-like exons must have a GC content lower than 49.3% and flanking introns longer than 691 nucleotides
* GC-like exons must have a GC content greater than 49.3% and at least one flanking intron smaller than 691 nucleotides
The script ``variance_analysis.py`` allow to test if control exons with a size > 27 nucleotides have a GC-content variation greater than the one of Irimia microexons [IWE14]_ or small fasterdb exons (with a size below 27 nucleotides).



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

  # Creation of correlation figures of differents features between exons and the other exons in the same TAD
  ## Creation of correlation figures
  python3 src/correlation_GC_exon_TAD_and_Exon/create_GC_AT_bed_exon.py # creation of a bed file containing same data about every regulated and control exons
  grep -v "'MFE_3SS': None" result/correlation_GC-AT-exons_TAD/data_for_regulated_CCE_exons.bed > result/correlation_GC-AT-exons_TAD/data_for_regulated_CCE_exons2.bed # removing 2 exons with some values missing due to annotation error.
  python3 src/correlation_GC_exon_TAD_and_Exon/tad_exon.py # creation of correlation figures
  ## Building models to test if some features are randmly distributed between tads
  # k562 TADs
  python3 src/correlation_GC_exon_TAD_and_Exon/branch_point_homogeneity_in_tad.py -t data/K562_Lieberman-raw_TADs.hg19.nochr.bed -e result/correlation_GC-AT-exons_TAD/data_for_regulated_CCE_exons2.bed -o result/correlation_GC-AT-exons_TAD/ -n table_TAD_K562_regulated_CCE_exons
  # MCF7 TADs
  python3 src/correlation_GC_exon_TAD_and_Exon/branch_point_homogeneity_in_tad.py -t data/HiC_TADs_Stein-MCF7-WT.hg19.nochr.bed -e result/correlation_GC-AT-exons_TAD/data_for_regulated_CCE_exons2.bed -o result/correlation_GC-AT-exons_TAD/ -n table_TAD_MCF7_regulated_CCE_exons

  # experimental branch points analysis
  mkdir result/experimental_branch_point
  python3 src/experimental_branch_point/branch_poin_bed_merged.py  # Create many bed files containing experimental branch points comming from different experiment filtered and merged together
  ## Creation of barplot representing the experimental branch points found in the 100-nucleotides sequences upstream the AT/GC and control exons.
  ### experimental branch points of mercer et al.
  python3 src/experimental_branch_point/GC_AT_analysis_experimental_bp.py \
  -b result/experimental_branch_point/mercer_2015_intron_end_filter.bed \
  -n mercer_intron_filtered -a result/AT_rich_exons \
  -g result/GC_rich_exons \
  -f data/fasterDB_lite.db  \
  -s data/sed.db \
  -o result/experimental_branch_point
  ### experimental branch points of pineda et al. with at least of coverage of 5 reads
  python3 src/experimental_branch_point/GC_AT_analysis_experimental_bp.py \
  -b result/experimental_branch_point/pineda_2018_5cov_filter_intron_end_filter.bed \
  -n pineda_5cov_intron_filtered -a result/AT_rich_exons \
  -g result/GC_rich_exons \
  -f data/fasterDB_lite.db  \
  -s data/sed.db \
  -o result/experimental_branch_point
  ### experimental branch points of taggart et al. with at least of coverage of 5 reads
  python3 src/experimental_branch_point/GC_AT_analysis_experimental_bp.py \
  -b result/experimental_branch_point/taggart_2017_5cov_filter_intron_end_filter.bed \
  -n taggart_5cov_intron_filtered -a result/AT_rich_exons \
  -g result/GC_rich_exons \
  -f data/fasterDB_lite.db  \
  -s data/sed.db \
  -o result/experimental_branch_point
  ### experimental branchpoints of mercer pineda and taggart
  python3 src/experimental_branch_point/GC_AT_analysis_experimental_bp.py \
  -b result/experimental_branch_point/merged_branch_point_all_intron_filter.bed \
  -n merged_intron_filtered -a result/AT_rich_exons \
  -g result/GC_rich_exons \
  -f data/fasterDB_lite.db  \
  -s data/sed.db \
  -o result/experimental_branch_point
  ### experimental branchpoints of pineda and taggart
  python3 src/experimental_branch_point/GC_AT_analysis_experimental_bp.py \
  -b result/experimental_branch_point/pineda-taggart_filtered_merged.bed \
  -n pineda-taggart_filtered -a result/AT_rich_exons \
  -g result/GC_rich_exons \
  -f data/fasterDB_lite.db  \
  -s data/sed.db \
  -o result/experimental_branch_point
  ## Checking the intersection between experimental and predicted branch points


  # Creation of the figures 2 and 3 of the article for custom list of exons
  mkdir result/figure_2_3
  python3.5 src/figures_2_3_creator/custom_exon_list_creator.py # Create many custom list of exons
  ## Creation of some folders
  mkdir result/figure_2_3/GA_CT_figures
  mkdir result/figure_2_3/otherGC_GC_figures
  mkdir result/figure_2_3/otherGC_AT_figures
  mkdir result/figure_2_3/unregulated_GC_AT_exons
  mkdir result/figure_2_3/otherGC_AT_GC_figures
  ## Creation of the figures 2 and 3 for custom list of exons
  ### For GA and CT exons
  python3.5 src/figures_2_3_creator/figure_2_3_creator.py \
  -l result/figure_2_3/exon_list/CT_rich_exons.txt result/figure_2_3/exon_list/GA_rich_exons.txt \
  -n CT-exons GA-exons \
  -N R \
  -s data/sed.db \
  -f data/fasterDB_lite.db \
  -o result/figure_2_3/GA_CT_figures \
  -b ../Clip_analysis/data/coverage_project_selected/ \
  -r ../Clip_analysis/data/hg19.ren.chrom.sizes \
  -m /media/nicolas/DD_1/Splicing_Lore_project/FarLine_exons_results_summary/src/skipped_exon_list_results_summary/coverage_summary/metagene_coverage.r \
  -S 'ct_rich_down' 'ga_rich_down' \
  --reverse y
  ### For GC/AT exons, SRSF2, SRSF3 and  hnRNPC down-regulated exons
  python3.5 src/figures_2_3_creator/figure_2_3_creator.py \
  -l result/GC_rich_exons result/AT_rich_exons result/figure_2_3/exon_list/other_GC_exons/SRSF2_all_exons.txt result/figure_2_3/exon_list/other_GC_exons/SRSF3_all_exons.txt result/figure_2_3/exon_list/other_GC_exons/HNRNPC_all_exons.txt \
  -n GC-exons AT-exons SRSF2-down-all SRSF3-down-all hnRNPC-down-all \
  -s data/sed.db \
  -f data/fasterDB_lite.db \
  -o result/figure_2_3/otherGC_AT_GC_figures \
  -b ../Clip_analysis/data/coverage_project_selected/ \
  -r ../Clip_analysis/data/hg19.ren.chrom.sizes \
  -m /media/nicolas/DD_1/Splicing_Lore_project/FarLine_exons_results_summary/src/skipped_exon_list_results_summary/coverage_summary/metagene_coverage.r
  ### For oAT and oGC exons. oAT exons are AT-like exons (see above) not regulated by any studied splicing factors. oGC exons are GC-like exons (see above) not regulated by any studied splicing factors.
  python3.5 src/figures_2_3_creator/figure_2_3_creator.py \
  -l result/figure_2_3/exon_list/GC_unregulated_exons.txt result/figure_2_3/exon_list/AT_unregulated_exons.txt \
  -n GC-unreg-exons AT-unreg-exons \
  -s data/sed.db \
  -f data/fasterDB_lite.db \
  -o result/figure_2_3/unregulated_GC_AT_exons \
  -b ../Clip_analysis/data/coverage_project_selected/ \
  -r ../Clip_analysis/data/hg19.ren.chrom.sizes \
  -m /media/nicolas/DD_1/Splicing_Lore_project/FarLine_exons_results_summary/src/skipped_exon_list_results_summary/coverage_summary/metagene_coverage.r
  ### Creation of addition figures for GC/AT exons, SRSF2, SRSF3 and  hnRNPC down-regulated exons
  mkdir result/figure_2_3/custom_figure_oGC_exons
  python3.5 src/figures_2_3_creator/custom_figure_creator.py \
  -l result/GC_rich_exons result/AT_rich_exons result/figure_2_3/exon_list/other_GC_exons/SRSF2_all_exons.txt result/figure_2_3/exon_list/other_GC_exons/SRSF3_all_exons.txt result/figure_2_3/exon_list/other_GC_exons/HNRNPC_all_exons.txt \
  -n GC-exons AT-exons SRSF2-down-all SRSF3-down-all hnRNPC-down-all \
  -s data/sed.db \
  -f data/fasterDB_lite.db \
  -o result/figure_2_3/custom_figure_oGC_exons \
  -t force_donor force_acceptor iupac_exon upstream_intron_size downstream_intron_size exon_size


  # microexon analysis
  python3 src/study_microexons/variance_analysis.py # check if microexon have the same GC-content variation than control exons.
  python3.5 src/study_microexons/bed_file_creator_from_sed.py # creation of 2 lists of exons : one containing small exons and the other big exons
  mkdir result/variance_analysis/enrichment_barplots
  ## checking if small exons are more AT or GC-like exons.
  python3.5 src/study_microexons/check_enrichment_in_GC_AT_group.py -l result/variance_analysis/bed_file/small_sf-downregulated_exons_\(3-27nt\).bed result/variance_analysis/bed_file/small_CCE_exons_\(3-27nt\).bed result/variance_analysis/bed_file/CCE_exons_27nt+.bed -b Small_down_exons Small_cce_exons CCE_exons -o result/variance_analysis/enrichment_barplots
  ##  checking if Irimia microexons are more AT or GC-like exons.
  mkdir result/variance_analysis/enrichment_Irimia_barplots
  python3.5 src/study_microexons/check_enrichment_in_GC_AT_group.py -l result/irimia_bed/Irimia_et_al_microexons_freq.bed result/variance_analysis/bed_file/CCE_exons_27nt+.bed -b Small_Irimia_exons CCE_exons -o result/variance_analysis/enrichment_Irimia_barplots


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
.. [MCA15] Mercer TR, Clark MB, Andersen SB, et al. Genome-wide discovery of human splicing branchpoints. Genome Res. 2015;25(2):290–303. doi:10.1101/gr.182899.114
.. [PiBr18] Pineda JMB, Bradley RK. Most human introns are recognized via multiple and tissue-specific branchpoints. Genes Dev. 2018;32(7-8):577–591. doi:10.1101/gad.312058.118
.. [TLS17] Taggart AJ, Lin CL, Shrestha B, Heintzelman C, Kim S, Fairbrother WG. Large-scale analysis of branchpoint usage across species and cell lines. Genome Res. 2017;27(4):639–649. doi:10.1101/gr.202820.115
.. [IWE14] Irimia M, Weatheritt RJ, Ellis JD, et al. A highly conserved program of neuronal microexons is misregulated in autistic brains. Cell. 2014;159(7):1511–1523. doi:10.1016/j.cell.2014.11.035
