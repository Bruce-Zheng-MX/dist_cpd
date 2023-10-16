#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --time=1:00:00
#SBATCH --array=1-2
#SBATCH --output=out/slurm_%j.out
#SBATCH --error=err/slurm_%j.err


module purge
module load gcc/11.3.0
module load openblas/0.3.20
module load r/4.2.3
#Rscript -e 'install.packages(c("kerSeg","gSeg","inline"), repos="https://cloud.r-project.org")'

for i in $(seq 1 $SLURM_NTASKS | head -n  $SLURM_NTASKS)
do
    line_to_read=$((SLURM_ARRAY_TASK_ID*$SLURM_NTASKS+$i))
    var=$(awk "NR==${line_to_read}" parameter.txt)
    read -r a b <<<$(echo $var)
    srun --ntasks=1 --nodes=1 --cpus-per-task=$SLURM_CPUS_PER_TASK Rscript distributional_data.R $a $b &
done
wait


#var=$(awk "NR==${SLURM_ARRAY_TASK_ID}" parameter.txt)
#read -r a b <<<$(echo $var)
#Rscript mul_normal.R $a $b
