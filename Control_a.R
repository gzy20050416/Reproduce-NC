# ---- 0) 环境 ----
library(Seurat)
library(Matrix)
library(data.table)
library(biomaRt)

set.seed(20)

raw_dir <- "~/NC/Ranalyze/raw"

# 简洁读表：第一列为feature名，其他列是细胞；转稀疏矩阵
read_count_matrix <- function(fp) {
  dt <- fread(fp)
  rn <- dt[[1]]
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- rn
  storage.mode(mat) <- "numeric"
  Matrix(mat, sparse = TRUE)
}

# 给列名统一追加样本后缀
add_suffix <- function(m, suf) {
  colnames(m) <- paste0(colnames(m), ":", suf)
  m
}

# ---- 1) 读取 RAW counts ----
# GEX (RNA)
gex_wt1 <- read_count_matrix(file.path(raw_dir, "GSM5130034_WT_1_singlecell_gex_raw_counts.txt.gz"))
gex_wt2 <- read_count_matrix(file.path(raw_dir, "GSM5130034_WT_2_singlecell_gex_raw_counts.txt.gz"))

# ADT
adt_wt1 <- read_count_matrix(file.path(raw_dir, "GSM5130035_WT_1_singlecell_adt_raw_counts.txt.gz"))
adt_wt2 <- read_count_matrix(file.path(raw_dir, "GSM5130035_WT_2_singlecell_adt_raw_counts.txt.gz"))

# HTO（作者就2个：HTO1、HTO2，通常是两行）
hto_wt1 <- read_count_matrix(file.path(raw_dir, "GSM5130036_WT_1_singlecell_hto_raw_counts.txt.gz"))
hto_wt2 <- read_count_matrix(file.path(raw_dir, "GSM5130036_WT_2_singlecell_hto_raw_counts.txt.gz"))

# ---- 2) 统一列名后缀，避免重名 ----
gex_wt1 <- add_suffix(gex_wt1, "WT1"); gex_wt2 <- add_suffix(gex_wt2, "WT2")
adt_wt1 <- add_suffix(adt_wt1, "WT1"); adt_wt2 <- add_suffix(adt_wt2, "WT2")
hto_wt1 <- add_suffix(hto_wt1, "WT1"); hto_wt2 <- add_suffix(hto_wt2, "WT2")

## ---- 3) 按样本建对象，并挂载 ADT/HTO（用列名交集对齐） ----

# WT1：取 gex/adt/hto 三者共同的细胞条形码，并统一顺序
cells_wt1 <- Reduce(intersect, list(colnames(gex_wt1), colnames(adt_wt1), colnames(hto_wt1)))
stopifnot(length(cells_wt1) > 0)
gex_wt1_s <- gex_wt1[, cells_wt1, drop = FALSE]
adt_wt1_s <- adt_wt1[, cells_wt1, drop = FALSE]
hto_wt1_s <- hto_wt1[, cells_wt1, drop = FALSE]

obj1 <- CreateSeuratObject(gex_wt1_s, assay = "RNA", project = "Bcell")
obj1[["ADT"]] <- CreateAssayObject(adt_wt1_s)
obj1[["HTO"]] <- CreateAssayObject(hto_wt1_s)
obj1$sample <- "WT1"

# WT2：同理
cells_wt2 <- Reduce(intersect, list(colnames(gex_wt2), colnames(adt_wt2), colnames(hto_wt2)))
stopifnot(length(cells_wt2) > 0)
gex_wt2_s <- gex_wt2[, cells_wt2, drop = FALSE]
adt_wt2_s <- adt_wt2[, cells_wt2, drop = FALSE]
hto_wt2_s <- hto_wt2[, cells_wt2, drop = FALSE]

obj2 <- CreateSeuratObject(gex_wt2_s, assay = "RNA", project = "Bcell")
obj2[["ADT"]] <- CreateAssayObject(adt_wt2_s)
obj2[["HTO"]] <- CreateAssayObject(hto_wt2_s)
obj2$sample <- "WT2"

# 合并
bcell <- merge(obj1, obj2)

# 可选：HTO重命名为 HTO1/HTO2（若作者文件行名不是这两个）
if (nrow(bcell[["HTO"]]) == 2 && !all(rownames(bcell[["HTO"]]) %in% c("HTO1","HTO2"))) {
  rownames(bcell[["HTO"]]) <- c("HTO1","HTO2")
}

cat("RNA:", nrow(bcell[["RNA"]]), "genes x", ncol(bcell), "cells\n")
cat("ADT features:", nrow(bcell[["ADT"]]), "\n")
cat("HTO features:", nrow(bcell[["HTO"]]), "\n")


# ---- 4) 合并为一个对象（保留 RNA/ADT/HTO 三个 assays）----
bcell <- merge(obj1, obj2)


# 简单回显
cat("RNA:", nrow(bcell[["RNA"]]), "genes x", ncol(bcell), "cells\n")
cat("ADT features:", nrow(bcell[["ADT"]]), "\n")
cat("HTO features:", nrow(bcell[["HTO"]]), "\n")



bcell <- PercentageFeatureSet(bcell, pattern = "^mt-", col.name = "percent_mito")
m1 <- quantile(as.numeric(bcell$percent_mito), 0.995)
bcell <- bcell[, bcell$percent_mito <= m1]

## 4) HTO 初次 CLR + Demux（用“标签个数”做 init） --------------------------
DefaultAssay(bcell) <- "HTO"



# 4.1 CLR 归一化（按细胞，margin = 2）
bcell <- NormalizeData(bcell, assay = "HTO",
                       normalization.method = "CLR",
                       margin = 2, verbose = FALSE)



## 5) ADT 归一化（CLR + Scale） ----------------------------------------------
bcell <- NormalizeData(bcell, assay = "ADT", normalization.method = "CLR", verbose = FALSE)
bcell <- ScaleData(bcell, assay = "ADT", verbose = FALSE)


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




