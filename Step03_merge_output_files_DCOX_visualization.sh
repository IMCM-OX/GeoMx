------------------------------
module add R/4.1.0-foss-2021a
module add R-bundle-Bioconductor/3.13-foss-2021a-R-4.1.0
------------------------------

#-------------------- Slurm

PI=GSE140231_RAW
DIR=/well/singlecell/projects_Isar/scRNAseq_Parkinson/RAW_DATA/$PI/DCOXA/DCOX_Plots/
mkdir -p $DIR'/scripts/'
cd $DIR'/scripts/'

scriptName=Step03_merge_output_files_DCOX_visualization  # generate an R script including the content of "Step04_flemix_DCOX_multiple_Generic.R"
j=1

# replace 66 with the output of length(listfiles) in "Step04_flemix_DCOX_multiple_Generic.R"
for i in {1..66}
do
echo $i

cp -fr $DIR'/scripts/'$scriptName'.R' $scriptName'_'$i'.R' -f
# Replace line x
sed -i '12s/i=1/i='$i'/' $scriptName'_'$i'.R'

done

#-------
rm -f $scriptName'.R'
#-------


#------- Creat .sh script
for f in $(find -type f -name '*.R'); do
echo $f

echo -e '#!/bin/bash
#SBATCH -p short
#SBATCH --chdir=./
#SBATCH --job-name=EACOX
#SBATCH --output=EACOX-%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G

module purge
module add R/4.1.0-foss-2021a
module add R-bundle-Bioconductor/3.13-foss-2021a-R-4.1.0

cd '$DIR'/scripts/
Rscript '${f}' ' > ${f}'.sh'

done

#------- Submit .sh script
for line in $(find -type f -name '*.sh'); do
 
sbatch $line;

echo $line
done

squeue -u ani054







#------- qsub


PI=P200233
DIR=/well/singlecell/projects_Isar/scRNAseq_Parkinson/RAW_DATA/$PI/10X-Expression_Informed_Demultiplex/DCOX_Analysis/
mkdir -p $DIR'/scripts/'
cd $DIR'/scripts/'

scriptName=P200233
j=1

# replace 66 with the output of length(listfiles) in "Step04_flemix_DCOX_multiple_Generic.R"
for i in {1..65}
do
echo $i

cp -fr $DIR'/scripts/'$scriptName'.R' $scriptName'_'$i'.R' -f
# Replace line x
sed -i '2s/Pair=1/Pair='$i'/' $scriptName'_'$i'.R'

done

#-------
rm -f $scriptName'.R'
#-------


#------- Creat .sh script
for f in $(find -type f -name '*.R'); do
echo $f

echo -e '#!/bin/bash
#$ -cwd
#$ -N DCOXA -j y
#$ -P htseq.prjc.low -q htseq.qd@@htseq.hgd
#$ -pe shmem 8

module purge
module add R/4.1.0-foss-2021a
module add R-bundle-Bioconductor/3.13-foss-2021a-R-4.1.0

cd '$DIR'/scripts/
Rscript '${f}' ' > ${f}'.sh'

done

#------- Submit .sh script
for line in $(find -type f -name '*.sh'); do
 
qsub $line;

echo $line
done


qstat
qstat | grep 'Script_' | wc

