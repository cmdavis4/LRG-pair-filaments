#$ -V
#$ -N matching
#$ -b y
#$ -cwd
#$ -e matching_out
#$ -o matching_out
#$ -l des
#$ -l h_vmem=4G
#$ -t 1-601

python neighbors_match.py