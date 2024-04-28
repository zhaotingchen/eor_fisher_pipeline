#!/bin/bash
#SBATCH --job-name=sim_eor_test_faint
#SBATCH --mail-type=END,FAIL          
#SBATCH --mail-user=zhaoting.chen@roe.ac.uk     
#SBATCH --time=144:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=7G
##SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
##SBATCH --mem-per-cpu=7G
##SBATCH --partition=Main
#SBATCH --array=37-46

SECONDS=0
echo ">>> activate conda"
source /home/zchen/.bashrc
source activate sim
module load cuda/11.8

cd rng_$SLURM_ARRAY_TASK_ID
#cd rng_44

echo $PWD

echo ">>> run 21cmfast"

python initialize.py

for i in {0..5..1}
do 
python run21cmfast.py faint $i HII_EFF_FACTOR
python run21cmfast.py faint $i ION_Tvir_MIN
done
python run21cmfast.py faint 6 HII_EFF_FACTOR
python signalcov.py faint


for i in {0..5..1}
do 
python run21cmfast.py bright $i HII_EFF_FACTOR
python run21cmfast.py bright $i ION_Tvir_MIN
done
python run21cmfast.py bright 6 HII_EFF_FACTOR
python signalcov.py bright


echo "****ELAPSED "$SECONDS"s sim_eor"