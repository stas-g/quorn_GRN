#Input data format needs to be transformed long data with time points as rows and genes as cols, this seems to struggle with big data (>10,000 cols) 

#Set up for SWING reverse network construction 
#conda create -n swing
#ca swing
#conda install python=3.7.0
#conda install -c anaconda numpy=1.16.5
#conda install -c anaconda scipy -y
#conda install -c anaconda pandas -y
#conda install -c anaconda scikit-learn -y 

#git clone https://github.com/bagherilab/SWING.git
#cd SWING
#python setup.py install

#########################################################################################
#screen -S swing 
#srun --partition=himem --cpus-per-task=100 --mem=1500G --pty bash 
#ca swing 
#cd /home/jconnell/projects/gene_regulatory_network/SWING
#python
#########################################################################################

from tabulate import tabulate
import os 
from Swing import Swing
import pandas as pd

os.listdir()
dataset = pd.read_csv('C1_swing.txt', sep='\t')
dataset.head()
dataset.shape
gene_start_column = 1
time_label = "Time"
separator = "\t"
gene_end = None
file_path = "C1_swing.txt"
k_min = 1
k_max = 3
w = 2
method = 'RandomForest'
trees = 100

sg = Swing(file_path, gene_start_column, gene_end, time_label, separator, min_lag=k_min, max_lag=k_max, window_width=w, window_type=method)
sg.zscore_all_data()
sg.create_windows()
sg.window_list[0]
sg.window_list[0].earlier_windows
sg.optimize_params()
sg.window_list[0].n_trees
sg.fit_windows(n_trees=trees, show_progress=False, n_jobs=-1)
sg.window_list[0].edge_importance
sg.compile_edges()
sg.full_edge_list

#col_names=["Parent", "Child", "Importance", "P_window", "C_window", "adj_imp", "Rank", "Edge", "Lag"]
#with open('sg.full_edge_list.txt', 'w') as f:
#    f.write(tabulate(sg.full_edge_list, headers=col_names))

#f1 = open('sg.full_edge_list2.txt', 'w')
#f1.write("\t".join(col_names))
#for x in sg.full_edge_list:
#	str_list = [str(a) for a in x]
#	f1.write("\t".join(str_list))
# f1.close()


sg.make_static_edge_dict(self_edges=False, lag_method='mean_mean')
ranked_edges = sg.make_sort_df(sg.edge_dict)
#ranked_edges
ranked_edges[['Source', 'Target']] = ranked_edges['regulator-target'].apply(pd.Series)

col_names=["regulator-target", "mean_importance", "Source", "Target"]
with open('ranked_edges.txt', 'w') as f:
    f.write(tabulate(ranked_edges, headers=col_names))

link_data = pd.read_fwf('ranked_edges.txt')
data = link_data[["Source", "Target", "mean_importance"]]
data = data.rename(columns={'Source': 'regulatory.gene', 'Target': 'target.gene', 'mean_importance': 'weight'})
data = data.iloc[1: , :]
data.to_csv('SWING_link_list.txt', index=False, sep='\t')




