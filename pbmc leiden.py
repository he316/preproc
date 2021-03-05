# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 14:10:22 2021

@author: user

scanpy tutorial leiden clustering and marker gene vol6.1 embedding with coordinates
##pbmc
<<<<<<< Updated upstream
##移除RAB37KO字串
=======
##scv pancrease
>>>>>>> Stashed changes
"""

# !mkdir data
# !wget http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz -O data/pbmc3k_filtered_gene_bc_matrices.tar.gz
# !cd data; tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
# !mkdir write

#sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
#sc.logging.print_header()
#sc.settings.set_figure_params(dpi=80, facecolor='white')

#results_file = 'write/pbmc3k.h5ad'  # the file that will store the analysis results

import scvelo as scv
import seaborn as sns
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import os
foldername="scv_pancrease"
Clustermethod="Leiden"
os.chdir("C:/Users/user/Desktop/test/scanpy")
#adata = sc.read(    'data/GSE132188/GSE132188_adata.h5ad.h5')   # the directory with the `.mtx` file                             # write a cache file for faster subsequent reading
adata = scv.datasets.pancreas()                             # write a cache file for faster subsequent reading
if(1<0):
    adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
    
    
    sc.pl.highest_expr_genes(adata, n_top=20,)
    
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    ##RAB37KO 的MT-是小寫,要從MT-改成mt-
    adata.var['mt'] = adata.var_names.str.startswith('mt-')
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.4, multi_panel=True)
    
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
    
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    
    sc.pp.log1p(adata)
    
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    
    sc.pl.highly_variable_genes(adata)
    ### Principal component analysis
    
    
    sc.tl.pca(adata, svd_solver='arpack')
    
    
    sc.pl.pca_variance_ratio(adata, log=True)
    
    
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    
    sc.tl.umap(adata)

sc.tl.leiden(adata,resolution=0.41)

picture=sc.pl.umap(adata, color=['leiden'],size=20,return_fig=True)

try:
    os.mkdir("./"+foldername,755)
except:
    pass

# adata.obs.leiden[X] = 第X個細胞的leiden cluster
# adata.obs.leiden.index[X] = 第X個細胞的UMI
# len(adata.obs.leiden.index) = 有被分群的細胞數量
# len(set(adata.obs.leiden)) = 群數
cluster_size={'cells':[]}
cluster_size=pd.DataFrame(cluster_size)
for index_m,elemenet_i in enumerate(set(adata.obs.leiden)): #從群開始分
    tempDF=pd.DataFrame(columns=adata[0].var.index.to_list()) #建立一個空的DF,每個cluster一個,作為存檔用
    for index_n,element_j in enumerate(adata.obs.leiden):#依序查詢細胞的分群        
        if(str(adata.obs.leiden[index_n]) == str(elemenet_i) ):#cluster符合則將該筆資料加入tempDF,方便寫入檔案
            tempDF=tempDF.append(adata[index_n].to_df())
    tempDF.to_csv('./'+foldername+'/'+Clustermethod+'_cluster_'+str(elemenet_i)+'.csv') #+'_'+str(len(tempDF.index))+'_cells
    cluster_size.loc[Clustermethod+'_cluster_'+str(elemenet_i)] = str(len(tempDF.index))
    print(foldername+'_'+Clustermethod+'_cluster_'+str(elemenet_i)+'_'+str(len(tempDF.index))+'_cells')
cluster_size=cluster_size.sort_index()
cluster_size.to_csv('./'+foldername+'/'+Clustermethod+'_clustering_size.csv')

adata.to_df().to_csv('./'+foldername+'/preprocessed_cell.csv')#不分cluster的資料
tempDF=pd.DataFrame(adata.obs.leiden)
tempDF["UMAP1"]=adata.obsm['X_umap'][:,0].tolist()
tempDF["UMAP2"]=adata.obsm['X_umap'][:,1].tolist()
tempDF.to_csv('./'+foldername+'/UMAP_cell_embeddings_to_'+Clustermethod+'_clusters_and_coordinates.csv')
#細胞分群和他們的座標

