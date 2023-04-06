#!/bin/bash

#SBATCH --job-name="tri_aracne"
#SBATCH --export=ALL

vst=/home/ngrinber/quorn_grn/

for x in {1..4}; do 
	outdir=/home/ngrinber/quorn_grn/aracne_bootstraps/tri_genes_cond_$x
	mkdir -p $outdir
	java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
	-e ${vst}/L2FC_filtered_data/L2FC_filter_condition_tri_$x.txt \
	-o $outdir \
	--pvalue 1E-8 \
	--seed 123 \
	--calculateThreshold
	####Run aracne on bootstraps of the input matrix
	for i in {1..100}
        do
            java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
            -e ${vst}/L2FC_filtered_data/L2FC_filter_condition_tri_$x.txt \
            -o $outdir \
            --pvalue 1E-8 \
            -t /home/ngrinber/quorn_grn/L2FC_filtered_data/tri_regs.txt \
            --seed $i 
        done
	####Consolidate, i.e. combine the bootstraps into a final network file
	java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
	-o $outdir \
	--consolidate
done 