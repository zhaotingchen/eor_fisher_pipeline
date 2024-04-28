#!/bin/bash
#SBATCH --job-name=sim_eor_test_faint
#SBATCH --mail-type=END,FAIL          
#SBATCH --mail-user=zhaoting.chen@roe.ac.uk     
#SBATCH --time=144:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=7G
#SBATCH --array=37-46

SECONDS=0
echo ">>> activate conda"
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

# the loop size should correspond to the length of the step_list
for i in {0..5..1}
do 
python run21cmfast.py faint $i HII_EFF_FACTOR
python run21cmfast.py faint $i ION_Tvir_MIN
done
# run fiducial
python run21cmfast.py faint 6 HII_EFF_FACTOR
# calculate the covariances for Fisher analysis
python signalcov.py faint

# same but for the bright model
for i in {0..5..1}
do 
python run21cmfast.py bright $i HII_EFF_FACTOR
python run21cmfast.py bright $i ION_Tvir_MIN
done
python run21cmfast.py bright 6 HII_EFF_FACTOR
python signalcov.py bright


echo "****ELAPSED "$SECONDS"s sim_eor"