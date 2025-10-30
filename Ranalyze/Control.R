#!/usr/bin/env Rscript
# =========================================
# B-cell 单细胞处理主流程
# =========================================

## 0) 环境与依赖 --------------------------------------------------------------
setwd("~/NC/Ranalyze")
set.seed(20)

library(Seurat)
library(Matrix)
library(dplyr)
library(biomaRt)
data_dir <- "~/NC/Bcell_counts/outs/filtered_feature_bc_matrix"
# ---------------------------
# Step 1. 读取10x矩阵
# ---------------------------
bcell_data <- Read10X(data.dir = data_dir)
names(bcell_data)                      # 预期: "Gene Expression" "Antibody Capture"

rna_mat <- bcell_data[["Gene Expression"]]
ab_mat  <- bcell_data[["Antibody Capture"]]

cat("RNA矩阵：", nrow(rna_mat), "基因 ×", ncol(rna_mat), "细胞\n")
cat("抗体矩阵：", nrow(ab_mat), "标签 ×", ncol(ab_mat), "细胞\n")

# ---------------------------
# Step 2. 根据名称拆分 HTO 与 ADT
# ---------------------------
adt_features <- grep("^Hashtag", rownames(ab_mat), invert = TRUE, value = TRUE)
hto_features <- grep("^Hashtag", rownames(ab_mat), value = TRUE)

adt_mat <- ab_mat[adt_features, , drop = FALSE]
hto_mat <- ab_mat[hto_features, , drop = FALSE]

rownames(adt_mat) <- gsub("_", "-", rownames(adt_mat))
rownames(hto_mat) <- gsub("_", "-", rownames(hto_mat))

cat("ADT特征：", length(adt_features), "；HTO特征：", length(hto_features), "\n")

# Step 3. 创建Seurat对象
bcell <- CreateSeuratObject(
  counts       = rna_mat,
  project      = "Bcell",
  min.cells    = 3,
  min.features = 200
)



# 添加 ADT 与 HTO assay）
bcell[["ADT"]] <- CreateAssayObject(counts = adt_mat)
bcell[["HTO"]] <- CreateAssayObject(counts = hto_mat)

bcell <- PercentageFeatureSet(bcell, pattern = "^mt-", col.name = "percent_mito")
m1 <- quantile(as.numeric(bcell$percent_mito), 0.995)
bcell <- bcell[, bcell$percent_mito <= m1]

## 4) HTO 初次 CLR + Demux --------------------------
DefaultAssay(bcell) <- "HTO"

# 4.0 剔除 HTO 全 0 细胞
hto_counts <- tryCatch(
  GetAssayData(bcell, assay = "HTO", layer = "counts"),
  error = function(e) GetAssayData(bcell, assay = "HTO", slot = "counts")
)
nz <- Matrix::colSums(hto_counts) > 0
bcell <- subset(bcell, cells = colnames(bcell)[nz])

# 4.1 CLR 归一化
bcell <- NormalizeData(bcell, assay = "HTO",
                       normalization.method = "CLR",
                       margin = 2, verbose = FALSE)

# 4.2 设定 init = HTO 标签数
k_init <- nrow(GetAssayData(bcell, assay = "HTO", layer = "counts"))

# 4.3 Demux
bcell <- HTODemux(bcell, assay = "HTO",
                  positive.quantile = 0.99,
                  kfunc = "clara", init = k_init)


## 5) ADT 归一化（CLR + Scale） ----------------------------------------------
if ("ADT" %in% Assays(bcell)) {
  bcell <- NormalizeData(bcell, assay = "ADT", normalization.method = "CLR", verbose = FALSE)
  bcell <- ScaleData(bcell, assay = "ADT", verbose = FALSE)
}

## 6) SCT 归一化（回归 mt） + PCA/UMAP/TSNE + 邻居/聚类 -----------------------
options(future.globals.maxSize = 3 * 1024^3)  
bcell <- SCTransform(bcell, assay = "RNA", vars.to.regress = "percent_mito", verbose = FALSE)

bcell <- RunPCA(bcell, assay = "SCT", verbose = FALSE)
DefaultAssay(bcell) <- "SCT"
bcell <- RunUMAP(bcell, dims = 1:30, verbose = FALSE)
bcell <- RunTSNE(bcell, dims = 1:30, verbose = FALSE)

bcell <- FindNeighbors(bcell, reduction = "pca", dims = 1:30)
bcell <- FindClusters(bcell, resolution = 0.4)



## 7) 细胞周期打分 + 细胞周期 PCA ---------------------------------
# 使用 Ensembl 存档镜像来建立连接
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

cell_cycle_human <- cc.genes

s_map <- getLDS(attributes = "hgnc_symbol", 
                filters = "hgnc_symbol", 
                values = cell_cycle_human$s.genes, 
                mart = human, 
                attributesL = "mgi_symbol", 
                martL = mouse, 
                uniqueRows = TRUE)

g2m_map <- getLDS(attributes = "hgnc_symbol", 
                  filters = "hgnc_symbol", 
                  values = cell_cycle_human$g2m.genes, 
                  mart = human, 
                  attributesL = "mgi_symbol", 
                  martL = mouse, 
                  uniqueRows = TRUE)

# 取小鼠符号列
s_genes_mouse <- unique(s_map$MGI.symbol)
g2m_genes_mouse <- unique(g2m_map$MGI.symbol)

# 打分并做 PCA
bcell <- CellCycleScoring(bcell, 
                          s.features = s_genes_mouse, 
                          g2m.features = g2m_genes_mouse, 
                          set.ident = FALSE)

bcell <- RunPCA(bcell, 
                assay = "SCT", 
                features = c(s_genes_mouse, g2m_genes_mouse), 
                reduction.name = "pca_cell_cycle", 
                reduction.key = "cellcyclePC_")

saveRDS(bcell, "bcell.rds")
