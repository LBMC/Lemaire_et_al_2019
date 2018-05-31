# Projects Figure Creator

## Description

This program aims to create a lots of figures that shows for each exons regulated in a project their distribution for many caracteristics described in sed database.

## Prerequisites

You must create Sed database.
The  code to create sed database is detailled here: https://gitlab.biologie.ens-lyon.fr/Auboeuf/uniform_exons_features/database_creator.
Then put the links of the sed database in the data directory of this project:

```sh
# Go to Figure_ESA folder
# Make sure to create the result directory where the result will be created.
mkdir result
# Creation of the folder data were the sed database should be located.
mkdir data/
cd data/
ln -s path_to_sed.db/sed.db .

```


This program works with `python3.5` and uses the following modules:
  * `sqlite3` : Connexion to sed database
  * `plotly` : Creation of the plots
  *  `os`
