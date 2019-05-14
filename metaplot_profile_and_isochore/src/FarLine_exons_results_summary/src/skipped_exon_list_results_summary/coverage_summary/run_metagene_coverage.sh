#!/bin/bash
#
### variable SGE
### shell du job
#$ -S /bin/bash
### nom du job
#$ -N metagene_coverage
### file d'attente
###$ -q monointeldeb48
#$ -q h48-E5-2667v2deb128,h6-E5-2667v4deb128,monointeldeb48,x5*,,E5-2670deb128*,E5-2667v2deb128nl,E5-2667v4deb256A
### charger l'environnement utilisateur pour SGE
#$ -cwd
### exporter les variables d'environnement sur tous les noeuds d'execution
#$ -V
### mails en debut et fin d'execution
###$ -m be
### change logs folder
#$ -o ppp_logs_out_dir_ppp
#$ -e ppp_logs_err_dir_ppp
#echo $dir
#
## Modif Sébastien Lemaire 15/04/2016
# crée dossier fastqc pour contenir toutes les analyses

# load module
source /usr/share/lmod/lmod/init/bash
export PYTHONPATH="${HOME}/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/src/gene_metaVisu/pyBigWig/lib/python3.6/site-packages":$PYTHONPATH
export PYTHONPATH="${HOME}/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/src/gene_metaVisu/pathos/lib/python3.6/site-packages":$PYTHONPATH
module load Python/3.6.1
umask 002

cd ${SGE_O_WORKDIR}
ls /scratch/ > /dev/null

#~ cmd="fastqc -t 1 -o ppp_output_dir_ppp/ ppp_input_ppp"
cmd="ppp_cmd_ppp"


echo "$cmd"
(>&2 echo ">>> command logs")
eval "$cmd"
(>&2 echo "<<< end command logs")

echo "END"
