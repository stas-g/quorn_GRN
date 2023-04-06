#!/usr/bin/env python

####Plot weights of GRN linklist genes 

####Import functions 
import argparse
import seaborn as sns 
import pandas as pd 
import matplotlib.pyplot as plt

####Parse args 
ap=argparse.ArgumentParser()
ap.add_argument('-i',type=str,required=True,help="linklist")
ap.add_argument('-t',type=str,help='plot title')
args=ap.parse_args()
title=args.t

###TRI gene names
####TRI6/10 are the regulators of the "tri cluster" 
####TRI4:g6429,TRI5:g6431,TRI6:g6430,TRI8:g6426,TRI10:g6432,TRI11:g6434,TRI12:g6435,TRI13:g6436,TRI14:g6437,
tri_regs=["g6430","g6432"]
tri_genes=["g6429","g6431","g6426","g6434","g6435","g6436","g6437"]

####Match TRI data and subset  
data=[x.replace("\n","").split("\t") for x in open(args.i)]
subset_data_1=[x[0:] for x in data for gene in tri_regs if x[0] == gene]
subset_data_2=[x[0:] for x in subset_data_1 for gene in tri_genes if x[1] == gene]

#print (subset_data_1)
#print ("***********************************************************")
#print (subset_data_2)

####Transform to dataframe and subset  
df=pd.DataFrame(subset_data_2)
df.columns=['Regulatory_gene', 'Target_gene', 'MI', 'P-value']
plot_data=df[['Target_gene', 'MI']]
plot_data=plot_data.astype({'MI':'float'})

####Plot barchart 
plt.figure(figsize=(20,10))
hst=sns.barplot(data=plot_data,
	x="Target_gene", 
	y="MI")
plt.suptitle(title, fontsize =25)
plt.xticks(rotation = 'vertical')

####Save Plot 
plt.savefig(title+"linklist_heatmap.png")
