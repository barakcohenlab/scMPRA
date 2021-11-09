#!/bin/bash
#
#SBATCH --job-name=ct2_no_min
#SBATCH --mail-type=ALL
#SBATCH --mail-user=szhao@wustl.edu
#
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=16000
#SBATCH --output=scMPRA_k562_r1.out

ml miniconda3
eval "$(conda shell.bash hook)"
conda activate mascot 

python3 210923_quad_to_pBC_mutiple_exp_calculation_copy.py --quint=../mixed_cell_rc_2_rep1.csv --promlib=../../210903_mixed_cell_bcs/ref_data/prom_lib_master_info.tsv --out=rep1 --exp=mixed_cell
