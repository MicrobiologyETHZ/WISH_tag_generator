#$ -cwd                   # run in current directory
#$ -S /bin/bash           # Interpreting shell for the job
#$ -N tnseq_tags          # Name of the job
#$ -V                     # .bashrc is read in all nodes
#$ -pe smp 32             # number of job slots to be reserved
#$ -l h_vmem=2G           # memory required
#$ -e data/error2.log           # error file
#$ -o data/out2.log             # output file

ml Python
ml Primer3
ml BWA

python recreate_reverse_primer.py
