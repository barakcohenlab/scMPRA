#!/usr/bin/env bash

#SBATCH -v -c 8 --mem=64000 -o cell-ranger.out -e cell-ranger.err

ml cellranger/6.0.1
cellranger count --id=mix_cells_r1 \
                 --fastqs=/scratch/bclab/szhao/210923_mix_txn_rep1/reads \
                 --transcriptome=/scratch/bclab/szhao/210923_mix_txn_rep1/cellranger_index/refdata-gex-GRCh38-2020-A/ \
                 --sample=Cohen_SZ-SL2-1_SI-TT-H9_AGAACTTAGA_AAAGGACTCG \
                 --localcores=8 \
                 --expect-cells=3000
