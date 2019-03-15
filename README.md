# Script of the article of **LEMAIRE *et al*** (2019)

This respository contains the script used to produce the analysis in the article of **LEMAIRE *et al*** entitled :

**<center>Splicing factors and chromatin organization enhance exon recognition
by alleviating constraints generated by gene nucleotide composition bias </center>**

It contains many subdirectories. Each subdirectories contains a `doc` folder that contains a detailled description of the scripts in their `src` folder.


Here is a short description of what every subdirectories contains:

* `Database_creator` : the scripts under its `src` directory, allow two create two `sqlite3` databases `fasterDB_lite.db` and `sed.db`. Those databases are used by the scripts in `Clip_analysis`, `Figure_ESA`, `Figures_SL` and `GC_rich_AT_rich_exon_list`.
* `Figure_ESA` :  the scripts under its `src` directory, allow two create the heatmaps and the correlation figures in the article.
* `Figures_SL` :  the scripts under its `src` directory allow you to create some suplementary figures of the article (i.e the figures representing a lot of boxplot : one boxplot for each splicing factor)
* `Clip_analysis` : the scripts under its `src` directory allow you to create the clip density figures.
* `GC_rich_AT_rich_exon_list` : the scripts under its `src` directory allow to create the boxplots, barplots and violin plots of many differents features (minimum free energy, number of branch point, number of UNA and T-rich low complexity sequences) for the AT- GC- U1- and U2-exons group. (see the method in the article for more precision)
