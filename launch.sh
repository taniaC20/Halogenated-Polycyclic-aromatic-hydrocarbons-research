#!/bin/bash

#SBATCH --account=def-jobstk-ab
#SBATCH --gres=gpu:v100:1
#SBATCH --cpus-per-task=10
#SBATCH --mem=24G
#SBATCH --time=1-00:00 # time (DD-HH:MM)

#SBATCH --job-name=PXA_RES_-%j
#SBATCH --output=OUT_%x-%j.txt

#SBATCH --output=OUT_%N-%j.txt # standard output
#SBATCH --error=ERR_%N-%j.txt  # standard error

#SBATCH --mail-user=tnchevez@mun.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL
#SBATCH --array=0-13


module load python/3.6

# virtualenv --no-download ~/pxa_env
source ~/pxa_env/bin/activate

# pip install --no-index --upgrade pip
# pip install numpy --no-index
# pip install pandas --no-index
# pip install progressbar2 --no-index
# pip install matplotlib --no-index

python ~/aromatic_batch/Aromatic_Compounds.py --index $SLURM_ARRAY_TASK_ID



