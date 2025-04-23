for i in x*
do
ETERNAFOLD_PATH='/home/akshay/apps/miniconda3/envs/py311/bin/eternafold-bin' ETERNAFOLD_PARAMETERS='/home/akshay/apps/EternaFold/parameters/EternaFoldParams.v1' python run_mrna_metrics.py ${i}  >> ${i}.metrics &
done
