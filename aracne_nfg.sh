vst=/home/ngrinber/quorn_grn/

for x in {1..4}; do 
	outdir=/home/ngrinber/quorn_grn/aracne_bootstraps/cond_$x
	mkdir -p $outdir
	java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
	-e ${vst}/L2FC_filtered_data/L2FC_filter_condition_$x.txt \
	-o $outdir \
	--pvalue 1E-8 \
	--calculateThreshold
	####Run aracne on bootstraps of the input matrix
	for i in {1..100}
	do
        java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
        -e ${vst}/L2FC_filter_condition_$x.txt \
        -o $outdir \
        --pvalue 01E-8 \
        -t ${vst}/transcription_factors_conditon_$x.txt \
        --seed $i 
	done
	####Consolidate, i.e. combine the bootstraps into a final network file
	java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
	-o $outdir \
	--consolidate
	####Plot tri cluster genes  
done 

java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-e /home/ngrinber/quorn_grn//L2FC_filtered_data/L2FC_filter_condition_1.txt \
-o /home/ngrinber/quorn_grn/aracne_test/bootstrap \
--pvalue 1E-8 \
--seed 1 \
--calculateThreshold

cp $vst/{L2FC_filter_condition_,transcription_factors_conditon_}{1..4}.txt /home/ngrinber/quorn_grn/