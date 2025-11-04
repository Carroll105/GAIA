import scanpy as sc
import copy
import numpy as np

def module(adata, components = 1,**kwargs):
    
    gene_module1 = adata.uns['Fisher_gene_module'][str(components-1)]['up']
    gene_module2 = adata.uns['Fisher_gene_module'][str(components-1)]['down']

    sc.tl.score_genes(adata, gene_list=gene_module1)
    gene_module1_score = adata.obs['score'].copy()
    sc.tl.score_genes(adata, gene_list=gene_module2)
    gene_module2_score = adata.obs['score'].copy()

    gene_module_score = np.vstack([np.array(gene_module1_score),np.array(gene_module2_score)])
    # print(gene_module_score)
    adata.obsm['X_gene_module'] = gene_module_score.T

    sc.pl.embedding(adata, basis='X_gene_module' , **kwargs)