#!/bin/bash

#SBATCH --job-name="tri_aracne"
#SBATCH --export=ALL

# outdir=/home/ngrinber/quorn_grn/aracne_bootstraps/vst
# mkdir -p $outdir

# java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
# -e /home/ngrinber/quorn_grn/L2FC_filtered_data/vst_L2FC_filtered.txt \
# -o $outdir \
# --pvalue 1E-8 \
# --seed 10 \
# --calculateThreshold
# ####Run aracne on bootstraps of the input matrix
# for i in {1..100}
#     do
#         java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
#         -e /home/ngrinber/quorn_grn/L2FC_filtered_data/vst_L2FC_filtered.txt \
#         -o $outdir \
#         --pvalue 1E-8 \
#         -t /home/ngrinber/quorn_grn/L2FC_filtered_data/tfs_combined.txt \
#         --seed $i 
#     done
# ####Consolidate, i.e. combine the bootstraps into a final network file
# java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
# -o $outdir \
# --consolidate



outdir=/home/ngrinber/quorn_grn/aracne_bootstraps/vst_nobase
mkdir -p $outdir

java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-e /home/ngrinber/quorn_grn/L2FC_filtered_data/vst_L2FC_filtered_nobase.txt \
-o $outdir \
--pvalue 1E-8 \
--seed 10 \
--calculateThreshold
####Run aracne on bootstraps of the input matrix
for i in {1..100}
    do
        java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
        -e /home/ngrinber/quorn_grn/L2FC_filtered_data/vst_L2FC_filtered_nobase.txt \
        -o $outdir \
        --pvalue 1E-8 \
        -t /home/ngrinber/quorn_grn/L2FC_filtered_data/tfs_combined.txt \
        --seed $i 
    done
####Consolidate, i.e. combine the bootstraps into a final network file
java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-o $outdir \
--consolidate
