#$ -V
#$ -N fake_pairs
#$ -b y
#$ -cwd
#$ -e job_out
#$ -o job_out
#$ -l des
#$ -l h_vmem=3G
#$ -t 1-70

python make_fake_pair_cat.py