
------------------------------
module add R/4.1.0-foss-2021a
module add R-bundle-Bioconductor/3.13-foss-2021a-R-4.1.0
------------------------------

#-------
Pairs=190
#-------  
scriptName=Script_DCOX
DIR=/well/singlecell/projects_Isar/scRNAseq_Parkinson/RAW_DATA/GSE140231_RAW/DCOXA/
mkdir -p $DIR'/scripts/'

cd $DIR'/scripts/'

for i in `seq 1 $Pairs`
do
    echo i in seq is $i

cp -fr 'Step02_flemix_DCOX_multiple_Generic.R' $DIR'/scripts/'$scriptName'P'$i'.R'

sed -i '5s/Pair=1/Pair='$i'/' $DIR'/scripts/'$scriptName'P'$i'.R'

done

#-------
rm -f 'Step02_flemix_DCOX_multiple_Generic.R'
#-------

#------- Slurm

cd '/well/singlecell/projects_Isar/scRNAseq_Parkinson/RAW_DATA/GSE140231_RAW/DCOXA//scripts/'

for f in $(find -type f -name '*.R'); do
echo $f

echo -e '#!/bin/bash
#SBATCH -p long
#SBATCH --chdir=./
#SBATCH --job-name=saver
#SBATCH --output=saver-%A_%a.out
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=20G

module purge
module add R/4.1.0-foss-2021a
module add R-bundle-Bioconductor/3.13-foss-2021a-R-4.1.0

cd /well/singlecell/projects_Isar/scRNAseq_Parkinson/RAW_DATA/GSE140231_RAW/DCOXA//scripts/

Rscript '${f}' ' > ${f}'.sh'

done

#-------

for line in $(find -type f -name '*.sh'); do
sbatch $line;
echo $line
done

squeue -u ani054



======== qsub ========

#----- Creat .sh script
cd '/gpfs3/well/singlecell/projects_Isar/scRNAseq_Parkinson/output/DCOXA_visium/SAVER/scripts/'

for f in $(find -type f -name '*.R'); do
echo $f

echo -e '#$ -cwd
module purge
module add R/4.1.0-foss-2021a
module add R-bundle-Bioconductor/3.13-foss-2021a-R-4.1.0

cd /gpfs3/well/singlecell/projects_Isar/scRNAseq_Parkinson/output/DCOXA_visium/SAVER/scripts/

Rscript '${f}' ' > ${f}'.sh'

done

#-------
cd $DIR'/scripts/'
# Submit .sh script
for line in $(find -type f -name '*'$scriptName'*.sh'); do
 
qsub -P SAVER $line;

echo $line
done

qstat

