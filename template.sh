#!/bin/bash
# change into your own slurm configuration
#SBATCH --job-name=sim_eor_test_faint
#SBATCH --mail-type=END,FAIL          
#SBATCH --mail-user=zhaoting.chen@roe.ac.uk   
#SBATCH --time=144:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=7G
# do not delete this line

SECONDS=0
echo ">>> activate conda"
# change to your own environment
# activate bash
source /home/zchen/.bashrc
# activate conda environment for 21cmfast and other things
source activate sim
module load cuda/11.8

# in each realization
cd rng_$SLURM_ARRAY_TASK_ID
echo $PWD

echo ">>> run 21cmfast"

python initialize.py

#do not delete this line
do 
python run21cmfast.py faint $i HII_EFF_FACTOR
python run21cmfast.py faint $i ION_Tvir_MIN
done
#do not delete this line
python signalcov.py faint

#do not delete this line
do 
python run21cmfast.py bright $i HII_EFF_FACTOR
python run21cmfast.py bright $i ION_Tvir_MIN
done
#do not delete this line
python signalcov.py bright


echo "****ELAPSED "$SECONDS"s sim_eor"