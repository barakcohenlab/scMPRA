#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=8000
#SBATCH --output=rep1_cr_ec.out

ml miniconda3
eval "$(conda shell.bash hook)"
conda activate mascot 

python3 cell_bcs_filtering.py --quad=../scMPRA_mixed_rep1.tsv --cellBC=../cellBCs/mixed_r1_barcodes.tsv.gz  --doc_name=mixed_quint_r1 
