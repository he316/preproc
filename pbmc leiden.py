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
import warnings

from scipy.sparse import issparse, coo_matrix
import scipy as sp
def select_distances(dist, n_neighbors=None):
    D = dist.copy()
    n_counts = (D > 0).sum(1).A1 if issparse(D) else (D > 0).sum(1)
    n_neighbors = (
        n_counts.min() if n_neighbors is None else min(n_counts.min(), n_neighbors)
    )
    rows = np.where(n_counts > n_neighbors)[0]
    cumsum_neighs = np.insert(n_counts.cumsum(), 0, 0)
    dat = D.data

    for row in rows:
        n0, n1 = cumsum_neighs[row], cumsum_neighs[row + 1]
        rm_idx = n0 + dat[n0:n1].argsort()[n_neighbors:]
        dat[rm_idx] = 0
    D.eliminate_zeros()
    return D

def select_connectivities(connectivities, n_neighbors=None):
    C = connectivities.copy()
    n_counts = (C > 0).sum(1).A1 if issparse(C) else (C > 0).sum(1)
    n_neighbors = (
        n_counts.min() if n_neighbors is None else min(n_counts.min(), n_neighbors)
    )
    rows = np.where(n_counts > n_neighbors)[0]
    cumsum_neighs = np.insert(n_counts.cumsum(), 0, 0)
    dat = C.data

    for row in rows:
        n0, n1 = cumsum_neighs[row], cumsum_neighs[row + 1]
        rm_idx = n0 + dat[n0:n1].argsort()[::-1][n_neighbors:]
        dat[rm_idx] = 0
    C.eliminate_zeros()
    return C

def get_neighs(adata, mode="distances"):
    if hasattr(adata, "obsp") and mode in adata.obsp.keys():
        return adata.obsp[mode]
    elif "neighbors" in adata.uns.keys() and mode in adata.uns["neighbors"]:
        return adata.uns["neighbors"][mode]
    else:
        raise ValueError("The selected mode is not valid.")


def get_n_neighs(adata):
    return adata.uns.get("neighbors", {}).get("params", {}).get("n_neighbors", 0)


def get_connectivities(adata, mode="connectivities", n_neighbors=None):
    if "neighbors" in adata.uns.keys():
        C = get_neighs(adata, mode)
        if n_neighbors is not None and n_neighbors < get_n_neighs(adata):
            if mode == "connectivities":
                C = select_connectivities(C, n_neighbors)
            else:
                C = select_distances(C, n_neighbors)
        connectivities = C > 0
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            connectivities.setdiag(1)
            connectivities = connectivities.multiply(1.0 / connectivities.sum(1))
        return connectivities.tocsr().astype(np.float32)
    else:
        return None
foldername="scv_pancreas_impute"
Clustermethod="celltype"
impute=True

os.chdir("C:/Users/user/Desktop/test/scanpy")
#adata = sc.read(    'data/GSE132188/GSE132188_adata.h5ad.h5')   # the directory with the `.mtx` file                             # write a cache file for faster subsequent reading
adata = scv.datasets.pancreas()

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
imputed_adata=adata.copy()
if impute is True:
    
    adata_conn=get_connectivities(adata)
    imputed_adata=adata.copy()
    imputed_adata.X=sp.dot(adata_conn,adata.X)
    sc.tl.umap(imputed_adata)
# write a cache file for faster subsequent reading
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

#sc.pp.neighbors(imputed_adata, n_neighbors=10, n_pcs=40)


sc.tl.leiden(imputed_adata,resolution=1)

picture=sc.pl.umap(imputed_adata, color=['leiden'],size=30,return_fig=True)
scv.tl.velocity(imputed_adata)
scv.tl.velocity_graph(imputed_adata)
scv.pl.velocity_embedding_stream(imputed_adata, basis='umap')
scv.pl.velocity_embedding(imputed_adata, basis='umap', arrow_length=2, arrow_size=1.5, dpi=150)
scv.tl.recover_dynamics(imputed_adata)
scv.tl.velocity(imputed_adata, mode='dynamical')
scv.tl.velocity_graph(imputed_adata)
scv.tl.latent_time(imputed_adata)
scv.pl.scatter(imputed_adata, color='latent_time', color_map='gnuplot', size=80, colorbar=True)

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
for index_m,element_i in enumerate(set(adata.obs['leiden'])): #從群開始分
    tempDF=pd.DataFrame(columns=adata[0].var.index.to_list()) #建立一個空的DF,每個cluster一個,作為存檔用
    for index_n,element_j in enumerate(adata.obs['leiden']): #依序查詢細胞的分群        
        if(str(adata.obs.leiden[index_n]) == str(element_i) ):#cluster符合則將該筆資料加入tempDF,方便寫入檔案
            tempDF=tempDF.append(adata[index_n].to_df())
    tempDF.to_csv('./'+foldername+'/'+Clustermethod+'_cluster_'+str(element_i)+'.csv') #+'_'+str(len(tempDF.index))+'_cells
    cluster_size.loc[Clustermethod+'_cluster_'+str(element_i)] = str(len(tempDF.index))
    print(foldername+'_'+Clustermethod+'_cluster_'+str(element_i)+'_'+str(len(tempDF.index))+'_cells')
cluster_size=cluster_size.sort_index()
cluster_size.to_csv('./'+foldername+'/'+Clustermethod+'_clustering_size.csv')

adata.to_df().to_csv('./'+foldername+'/preprocessed_cell.csv')#不分cluster的資料
tempDF=pd.DataFrame(adata.obs.leiden)
tempDF["UMAP1"]=adata.obsm['X_umap'][:,0].tolist()
tempDF["UMAP2"]=adata.obsm['X_umap'][:,1].tolist()
tempDF.to_csv('./'+foldername+'/UMAP_cell_embeddings_to_'+Clustermethod+'_clusters_and_coordinates.csv')
#細胞分群和他們的座標



#ifcluster的名稱是字串就跑這一段
###<===start
clusters='clusters'

cluster_size={'cells':[]}
cluster_size=pd.DataFrame(cluster_size)

save_cluster_size={'cells':[]}
save_cluster_size=pd.DataFrame(save_cluster_size)
    
for index_m,element_i in enumerate(set(imputed_adata.obs['clusters'])):
    cluster_size.loc[str(element_i)] = len(imputed_adata.obs['clusters'][imputed_adata.obs['clusters']==element_i])
cluster_size=cluster_size.sort_values(by='cells',ascending=False)

for index_m,element_i in enumerate(cluster_size.index): #從群開始分
    tempDF=pd.DataFrame(columns=imputed_adata[0].var.index.to_list()) #建立一個空的DF,每個cluster一個,作為存檔用
    for index_n,element_j in enumerate(imputed_adata.obs['clusters']): #依序查詢細胞的分群        
        if(str(imputed_adata.obs.clusters[index_n]) == str(element_i) ):#cluster符合則將該筆資料加入tempDF,方便寫入檔案
            tempDF=tempDF.append(imputed_adata[index_n].to_df())
    tempDF.to_csv('./'+foldername+'/'+Clustermethod+'_cluster_'+str(index_m)+'.csv') #+'_'+str(len(tempDF.index))+'_cells
    save_cluster_size.loc[Clustermethod+'_cluster_'+str(index_m)] = str(len(tempDF.index))
    print(foldername+'_'+Clustermethod+'_cluster_'+str(index_m)+'_'+str(element_i)+'_'+str(len(tempDF.index))+'_cells')
cluster_size.to_csv('./'+foldername+'/'+Clustermethod+'_clustering_table.csv')
save_cluster_size.to_csv('./'+foldername+'/'+Clustermethod+'_clustering_size.csv')

imputed_adata.obs['clusters2num']=imputed_adata.obs['clusters']
#for index_m,element_i in enumerate(cluster_size.index):
##這邊要把替換要替換的cluster帶進來    
imputed_adata.rename_categories('clusters2num',['0','5','1','2','3','4','7','6'])

imputed_adata.to_df().to_csv('./'+foldername+'/preprocessed_cell.csv')#不分cluster的資料
tempDF=pd.DataFrame(imputed_adata.obs.clusters2num)
tempDF["UMAP1"]=imputed_adata.obsm['X_umap'][:,0].tolist()
tempDF["UMAP2"]=imputed_adata.obsm['X_umap'][:,1].tolist()
tempDF.to_csv('./'+foldername+'/UMAP_cell_embeddings_to_'+Clustermethod+'_clusters_and_coordinates.csv')
imputed_adata.write_h5ad('./'+foldername+'/'+foldername+'.h5ad')
###<===end
