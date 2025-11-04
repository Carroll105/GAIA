import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture
from geomstats.learning.frechet_mean import FrechetMean
from geomstats.learning.pca import TangentPCA
from geomstats.geometry.hypersphere import Hypersphere
# from scipy.spatial import distance_matrix
# import torch
# import torch.nn as nn
# import torch.optim as optim


def fisher_PCA(adata, components = 10):
    
    data = adata.copy()
    sc.pp.normalize_total(data, target_sum=1)
    data = data.to_df()
    sqrt_df = np.sqrt(data)
    sqrt_df.dropna(inplace=True)
    sphere = Hypersphere(dim=2)
    mean = FrechetMean(metric=sphere.metric)
    mean.fit(np.array(sqrt_df))
    mean_estimate = mean.estimate_

    tpca = TangentPCA(sphere.metric, n_components=components)
    tpca = tpca.fit(sqrt_df, base_point=mean_estimate)
    tangent_projected_data = tpca.transform(sqrt_df)

    adata.obsm['X_scFisher'] = tangent_projected_data
    adata.varm['loading_scFisher'] = tpca.components_.T


def gene_module(adata, components = 10):
    adata.uns['Fisher_gene_module'] = {}
    if 'loading_scFisher' in adata.varm:
        gene_weight = adata.varm['loading_scFisher'][:,:].T
    
    else:
        print("please first run fisher_PCA function")

    if(components > adata.varm['loading_scFisher'].shape[1]):

        print("there is not enough scFisher loadings, use the max loadings")
        components = adata.varm['loading_scFisher'].shape[1]

    for i in range(components):
        com_vec = gene_weight[i,:]
        gm = GaussianMixture(n_components=3, random_state=0,n_init = 10).fit_predict(com_vec.reshape(-1,1))        
        # print(gm)
        gm_df = pd.DataFrame(index = np.sort(np.unique(gm)))
        gm_df['score'] = 0
        # print(np.unique(gm))
        gm_df['score'].iloc[0] = com_vec[np.where(gm==0)[0]].mean()
        gm_df['score'].iloc[1] = com_vec[np.where(gm==1)[0]].mean()
        gm_df['score'].iloc[2] = com_vec[np.where(gm==2)[0]].mean()

        gm_df.sort_values(by = 'score', inplace = True)

        adata.uns['Fisher_gene_module'][str(i)] = {}
        adata.uns['Fisher_gene_module'][str(i)]['down'] = adata.var_names[np.where(gm == gm_df.index[0])]
        adata.uns['Fisher_gene_module'][str(i)]['up'] = adata.var_names[np.where(gm == gm_df.index[-1])]
    
    print("finished gene module")
# import scipy.sparse as sp
# def learn_heights(adata, gene_program, lambda_smooth=1.0, k_neighbors=5, n_iter=1000, lr=0.01, use_cuda=False):
#     """
#     """
#     device = torch.device("cuda" if use_cuda and torch.cuda.is_available() else "cpu")

#     coords = adata.obsm['X_spatial']
#     n_cells = adata.n_obs

#     # 2. 计算 gene program 激活值
#     genes_idx = [i for i, gene in enumerate(adata.var_names) if gene in gene_program]
#     if len(genes_idx) == 0:
#         raise ValueError("gene_program 中的基因均未在 adata.var_names 中找到！")
#     if hasattr(adata.X, "toarray"):
#         expr = adata.X.toarray()
#     else:
#         expr = adata.X
#     gene_expr = np.mean(expr[:, genes_idx], axis=1)

#     # 3. 构造空间 kNN 图
#     dist = distance_matrix(coords, coords)
#     knn_matrix = np.zeros((n_cells, n_cells), dtype=np.float32)
#     for i in range(n_cells):
#         neighbors = np.argsort(dist[i])[1:k_neighbors+1]
#         knn_matrix[i, neighbors] = 1
#     knn_matrix = sp.csr_matrix(knn_matrix)

#     # 转换为 PyTorch 稀疏矩阵
#     indices = torch.tensor(np.vstack((knn_matrix.nonzero())), dtype=torch.long, device=device)
#     values = torch.tensor(knn_matrix.data, dtype=torch.float32, device=device)
#     sparse_knn = torch.sparse_coo_tensor(indices, values, (n_cells, n_cells), device=device)

#     # 4. 训练模型
#     heights = torch.nn.Parameter(torch.randn(n_cells, device=device, requires_grad=True))
#     a = torch.nn.Parameter(torch.randn(1, device=device, requires_grad=True))
#     b = torch.nn.Parameter(torch.randn(1, device=device, requires_grad=True))

#     optimizer = optim.Adam([heights, a, b], lr=lr)
#     gene_expr_tensor = torch.tensor(gene_expr, dtype=torch.float32, device=device)

#     for it in range(n_iter):
#         optimizer.zero_grad()

#         # 预测 gene program 值
#         pred_gene = a * heights + b

#         # 计算 gene program 拟合损失
#         gene_loss = torch.mean((pred_gene - gene_expr_tensor) ** 2)

#         # 计算空间平滑损失（利用稀疏矩阵批量计算）
#         height_diff = heights.unsqueeze(1) - heights.unsqueeze(0)
#         smooth_loss = torch.sum((sparse_knn.to_dense() * height_diff) ** 2) / (n_cells * k_neighbors)

#         # 计算总损失并优化
#         loss = gene_loss + lambda_smooth * smooth_loss
#         loss.backward()
#         optimizer.step()

#         if it % 100 == 0:
#             print(f"Iteration {it}, Total Loss: {loss.item():.4f}, Gene Loss: {gene_loss.item():.4f}, Smooth Loss: {smooth_loss.item():.4f}")

#     adata.obs['height'] = heights.detach().cpu().numpy()
#     return adata

