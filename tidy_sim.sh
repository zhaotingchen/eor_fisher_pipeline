for i in {37..46..1}
do
cd rng_$i
sbatch run_rest.sh
cd ../
done