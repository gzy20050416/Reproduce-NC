# =========================
# B-cell 复现 · 精简统一版
# =========================
setwd("/data1/xty/zff/NC/Ranalyze")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(clusterProfiler)
  library(msigdbr)
})

# -------------------------
# 0) 数据与身份映射
# -------------------------
bcell <- readRDS("bcell.rds")
bcell  <- UpdateSeuratObject(bcell)

map_vec <- c(
  "0"="Kappa Pre B","1"="Immature B","2"="PreBCRi I","3"="Mature B",
  "4"="Pro B VDJ","5"="Lambda Pre B","6"="Cycling Pro B",
  "7"="PreBCRi II S phase","8"="PreBCRd","9"="Prepro B",
  "10"="PreBCRi II G2/M phase","11"="High Mitochondrial B",
  "12"="Plasma Cells","13"="Cycling Immature B"
)
ordered_levels <- c(
  "Prepro B","Cycling Pro B","Pro B VDJ","PreBCRd","PreBCRi I",
  "PreBCRi II S phase","PreBCRi II G2/M phase","Kappa Pre B","Lambda Pre B",
  "Immature B","Cycling Immature B","Mature B","Plasma Cells","High Mitochondrial B"
)
bcell$celltype <- unname(map_vec[as.character(bcell$seurat_clusters)])
Idents(bcell)  <- factor(bcell$celltype, levels = ordered_levels)

# -------------------------
# 1) 统一风格 & 工具函数
# -------------------------
theme_base <- theme_classic(base_size = 12) +
  theme(panel.grid = element_blank(), axis.line = element_line(color="black"))

# 通用保存（避免尺寸限制 + 嵌入字体更稳）
save_pdf <- function(p, file, w, h){
  ggsave(file, p, width = w, height = h, units = "in",
         limitsize = FALSE, useDingbats = FALSE)
}

# 计算 UMAP/TSNE/PCA 的合适宽高（保持坐标比例不变，避免“被压扁”）
.calc_dim_size <- function(obj, reduction, base = 6.8){
  emb <- Embeddings(obj, reduction = reduction)
  xr  <- diff(range(emb[,1])); yr <- diff(range(emb[,2]))
  asp <- ifelse(yr == 0, 1, xr / yr)
  if (is.na(asp) || is.infinite(asp) || asp <= 0) asp <- 1
  if (asp >= 1) c(w = base * asp, h = base) else c(w = base, h = base / asp)
}

# 降维单屏（自动尺寸 + 保持比例）
dim_basic <- function(obj, file, reduction="umap", group.by=NULL, title=NULL){
  sz <- .calc_dim_size(obj, reduction)
  p <- DimPlot(obj, reduction = reduction, label = FALSE, group.by = group.by) +
    ggtitle(title %||% "") + theme_base + coord_fixed()
  save_pdf(p, file, w = sz["w"], h = sz["h"]); p
}

# 高亮（自动尺寸 + 保持比例）
dim_highlight <- function(obj, colors_named, file, title=NULL){
  sz <- .calc_dim_size(obj, "umap")
  cols_all <- setNames(rep("#D9D9D9", length(levels(Idents(obj)))), levels(Idents(obj)))
  cols_all[names(colors_named)] <- unname(colors_named)
  p <- DimPlot(obj, reduction="umap", cols=cols_all, pt.size=0.55) +
    theme_base + theme(legend.position = "none") +
    ggtitle(title %||% "") + coord_fixed()
  save_pdf(p, file, w = sz["w"], h = sz["h"]); p
}

# FeaturePlot 网格（关键修复：combine=FALSE + wrap_plots + & theme）
# 自动根据特征个数和列数计算 PDF 宽高，统一标题对齐并收拢图例
fp_cols <- c("grey92", "#6A51A3")

# --- 固定尺寸 & 每面板自带色标 & 标题不跑位 ---
fp_grid <- function(obj, feats, file,
                    ncol, nrow = NULL,             
                    assay = "SCT",
                    pt = 0.45, minc = "q05", maxc = "q95",
                    legend_side = "right",          # 每个子图的色标位置
                    w, h) {                         
  
  DefaultAssay(obj) <- assay
  pp <- FeaturePlot(
    obj, features = feats, reduction = "umap", cols = fp_cols,
    order = TRUE, pt.size = pt, min.cutoff = minc, max.cutoff = maxc,
    combine = FALSE
  )
  
  # 逐子图统一样式；每个子图各自保留色标，标题居中；固定轴比例防止变形
  pp <- lapply(pp, function(p)
    p + theme_classic(base_size = 12) +
      theme(panel.grid = element_blank(),
            axis.line = element_line(color = "black"),
            plot.title = element_text(hjust = .5, margin = margin(b = 3)),
            plot.margin = margin(4, 4, 4, 4),
            legend.position = legend_side) +
      coord_fixed()
  )
  
  p <- patchwork::wrap_plots(pp, ncol = ncol, nrow = nrow)
  ggsave(file, p, width = w, height = h, units = "in",
         limitsize = FALSE, useDingbats = FALSE)
  invisible(p)
}


# 小提琴自动宽度（按组数×特征数线性放大）
vln_by_celltype <- function(obj, feats, file, assay="SCT", ncol=2, pt=0.2,
                            base_h = 4.8, per_group_w = 0.5){
  DefaultAssay(obj) <- assay
  ng <- length(levels(factor(obj$celltype)))
  w  <- max(7.0, ng * per_group_w * length(feats))  # 组越多、特征越多，自动变宽
  p <- VlnPlot(
    obj, features = feats, group.by = "celltype", pt.size = pt, ncol = ncol
  ) +
    labs(x = NULL, y = "Expression") +
    theme_base +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_pdf(p, file, w = w, h = base_h + 1.2); p
}


# 二组热图：
heatmap_two_groups <- function(obj, g1, g2, topN = 30, file,
                               assay = "SCT", w = 10, h = 8) {
  DefaultAssay(obj) <- assay
  sub <- subset(obj, idents = c(g1, g2))
  Idents(sub) <- factor(Idents(sub), c(g1, g2))
  
  m <- FindMarkers(sub, ident.1 = g1, ident.2 = g2,
                   logfc.threshold = 0.25, min.pct = 0.25)
  
  m_up   <- m[m$avg_log2FC > 0, , drop = FALSE]
  m_down <- m[m$avg_log2FC < 0, , drop = FALSE]
  
  up   <- head(rownames(m_up  [order(m_up$avg_log2FC,  decreasing = TRUE), , drop = FALSE]), topN)
  down <- head(rownames(m_down[order(m_down$avg_log2FC, decreasing = FALSE), , drop = FALSE]), topN)
  genes <- unique(c(up, down))
  genes <- genes[!is.na(genes)]
  genes <- intersect(genes, rownames(sub))
  sub <- ScaleData(sub, features = genes, verbose = FALSE)
  p <- DoHeatmap(sub, features = genes, assay = assay, slot = "scale.data",
                 group.by = "ident", group.bar = TRUE, draw.lines = TRUE, raster = TRUE) +
    scale_fill_gradientn(colors = c("#b44699", "#000000", "#ebe231")) +
    theme(axis.text.y = element_text(size = 6))
  ggsave(file, p, width = w, height = h)
  p
}

# 火山：5.3x4.1 -> 7x5.5
volcano_simple <- function(obj, ident1, ident2, file, assay="RNA",
                           xlim=NULL, ylim=NULL, labels_left=c(), labels_right=c(),
                           title=NULL, w=7, h=5.5){
  DefaultAssay(obj) <- assay
  de <- FindMarkers(obj, ident.1=ident1, ident.2=ident2, test.use="wilcox", logfc.threshold=0, min.pct=0.1)
  de$gene <- rownames(de); de$logP <- -log10(pmax(de$p_val, .Machine$double.eps))
  p <- ggplot(de, aes(avg_log2FC, logP)) +
    geom_point(size=1.6) + geom_vline(xintercept=0, linetype="dashed") +
    labs(x="avg_log2FC", y="-log10(p)", title=title %||% "") + theme_base
  if(!is.null(xlim)) p <- p + coord_cartesian(xlim=xlim, ylim=ylim)
  if(length(labels_left))  p <- p + ggrepel::geom_text_repel(data=subset(de, gene %in% labels_left),
                                                             aes(label=gene), size=3, max.overlaps=Inf)
  if(length(labels_right)) p <- p + ggrepel::geom_text_repel(data=subset(de, gene %in% labels_right),
                                                             aes(label=gene), size=3, max.overlaps=Inf)
  save_pdf(p, file, w, h); p
}

# GSEA 四联图：18x4.5 -> 22x6
gsea_4panel <- function(gsea_res, geneList, terms, file, w=22, h=6){
  m_df <- msigdbr(species="Mus musculus", category="H")[,c("gs_name","gene_symbol")]
  pathways <- split(m_df$gene_symbol, m_df$gs_name)
  run_es <- function(stats, genes){
    stats <- sort(stats, decreasing=TRUE)
    hit   <- names(stats) %in% genes
    Nm <- sum(!hit); if(!any(hit) || Nm==0) return(NULL)
    step <- ifelse(hit, abs(stats[hit])/sum(abs(stats[hit])), -1/Nm)
    es   <- cumsum(replace(rep(-1/Nm, length(stats)), which(hit), step[hit]))
    data.frame(rank=seq_along(es), es=es, hit=hit, stats=stats)
  }
  plot_one <- function(term){
    d <- run_es(-geneList, pathways[[term]]); if(is.null(d)) return(NULL)
    res <- as.data.frame(gsea_res@result); row <- res[match(term, res$ID),]
    NES <- row$NES; qv <- if("qvalues"%in%names(row)) row$qvalues else row$p.adjust
    top <- ggplot(d, aes(rank, stats)) + geom_area(fill="grey70") +
      labs(y="Ranked list metric", x=NULL) + theme_base +
      theme(axis.text.x=element_blank(), plot.margin=margin(5,5,0,5))
    bot <- ggplot(d, aes(rank, es)) + geom_line() +
      geom_vline(xintercept=which.max(abs(d$es)), linetype="dashed", color="red") +
      geom_segment(data=data.frame(rank=which(d$hit)),
                   aes(x=rank, xend=rank, y=0, yend=0.06), linewidth=.25) +
      labs(y="Running Enrichment Score",
           x="Position in the Ranked List of Genes",
           title=sprintf("%s\nNES=%.3f  q=%.2e", term, NES, qv)) +
      theme_base + theme(plot.margin=margin(0,5,5,5),
                         plot.title=element_text(hjust=.5, size=10))
    top / bot + plot_layout(heights=c(1,3))
  }
  p <- wrap_plots(lapply(terms, plot_one), nrow=1, guides="collect")
  save_pdf(p, file, w, h); p
}

`%||%` <- function(a,b) if(!is.null(a)) a else b

# -------------------------
# 2) Fig1：全局降维与基础 Feature/ADT
# -------------------------
dim_basic(bcell, "01_umap_clusters.pdf",      reduction="umap", title="UMAP clusters")
dim_basic(bcell, "01_umap_phase.pdf",         reduction="umap", group.by="Phase", title="UMAP Phase")
dim_basic(bcell, "01_tsne_clusters.pdf",      reduction="tsne", title="tSNE clusters")
dim_basic(bcell, "01_pca_phase.pdf",          reduction="pca",  group.by="Phase", title="PCA Phase")

DefaultAssay(bcell) <- "SCT"
p_mki67 <- FeaturePlot(bcell, features="Mki67", reduction="umap", cols=fp_cols,
                       order=TRUE, pt.size=0.5) + theme_base + ggtitle("Mki67")
save_pdf(p_mki67, "01_feature_Mki67.pdf", 7.5, 7.5)


DefaultAssay(bcell) <- "ADT"
fp_grid(bcell, c("adt_B220","adt_CD19","adt_CD93","adt_CD25","adt_IgM","adt_CD43"),
        "01_feature_ADT_panel.pdf" , ncol = 6, nrow = 1, assay = "ADT",w = 30, h = 6)

DefaultAssay(bcell) <- "SCT"
fp_grid(bcell, c("Ptprc","Cd19","Cd93","Il2ra","Ighm","Spn"),
        "01_feature_RNA_panel.pdf", ncol = 6, nrow = 1, assay = "SCT",w = 30, h = 6)

# -------------------------
# 3) Fig2：Prepro vs Cycling Pro · 高亮 / 早期标志 / 二组热图 / 轻链小提琴
# -------------------------
dim_highlight(bcell, c("Prepro B"="#d73027","Cycling Pro B"="#3b4cc0"),
              "02_highlight_prepro_vs_cyclingPro.pdf",
              title="pre-pro B / cycling pro B")

fp_grid(bcell, c("Flt3","Il7r","Cd79a"),
        "02_feature_earlyB_markers.pdf",ncol = 3, nrow = 1, assay = "SCT",w = 12, h = 4.2)

heatmap_two_groups(bcell, "Prepro B","Cycling Pro B", topN=30,
                   "02_heatmap_prepro_vs_cyclingPro.pdf", assay="SCT")

vln_by_celltype(bcell, c("Vpreb1","Igll1"),
                "02_violin_Vpreb1_Igll1_by_celltype.pdf")

# -------------------------
# 4) Fig3：PreBCRd 高亮 + PreBCR 相关 Feature +GSEA 四联图
# -------------------------
dim_highlight(bcell, c("PreBCRd"="#b2182b"),
              "03_highlight_preBCRd.pdf", title="PreBCR-dependent expansion")

fp_grid(
  bcell,
  c("Nrgn","Ybx1","Ybx3","Slc7a5","Slc3a2","Myc"),
  "03_feature_preBCR_markers.pdf",
  ncol = 3, nrow = 2, assay = "RNA",  
  w = 12, h = 8
)

if (exists("gsea_res") && exists("geneList")) {
  gsea_4panel(gsea_res, geneList,
              c("HALLMARK_MYC_TARGETS_V2","HALLMARK_GLYCOLYSIS",
                "HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_OXIDATIVE_PHOSPHORYLATION"),
              "03_hallmark_GSEA_4paths.pdf")
}

# -------------------------
# 5) Fig4：PreBCRi I/II 高亮 + ADT Violin + 单基因 UMAP/Violin + 四群热图 + 两张火山
# -------------------------
dim_highlight(bcell, c("PreBCRi I"="#b2182b","PreBCRi II S phase"="#2166ac","PreBCRi II G2/M phase"="#2166ac"),
              "04_highlight_preBCRi_groups.pdf", title="PreBCR independent proliferation")

DefaultAssay(bcell) <- "ADT"
vln_by_celltype(bcell, "adt_CD43", "04_violin_ADT_CD43.pdf", assay="ADT", ncol=1)
vln_by_celltype(bcell, "adt_CD25", "04_violin_ADT_CD25.pdf", assay="ADT", ncol=1)

DefaultAssay(bcell) <- "SCT"
fp_grid(bcell, c("Bach2","Cd74","Cxcr4","Cd44"),
        "04_feature_Bach2_Cd74_Cxcr4_Cd44.pdf" ,ncol = 2, nrow = 2, assay = "SCT",w = 9, h = 6)

# 四群热图：8x10 -> 10x12
Idents(bcell) <- factor(bcell$celltype,
                        levels=c("PreBCRd","PreBCRi I","PreBCRi II S phase","PreBCRi II G2/M phase"))
sub4 <- subset(bcell, idents=levels(Idents(bcell)))
get_top <- function(obj, grp, n=35){
  m <- FindMarkers(obj, ident.1=grp, ident.2=setdiff(levels(Idents(obj)), grp),
                   logfc.threshold=0.25, min.pct=0.1)
  m$gene <- rownames(m); m <- m[m$avg_log2FC>0 & m$p_val_adj<0.05,]
  head(m[order(m$p_val_adj, -m$avg_log2FC),"gene"], n)
}
genes4 <- unique(c(get_top(sub4,"PreBCRd"), get_top(sub4,"PreBCRi I"),
                   get_top(sub4,"PreBCRi II G2/M phase"), get_top(sub4,"PreBCRi II S phase")))
sub4 <- ScaleData(sub4, features=genes4, verbose=FALSE)
p4 <- DoHeatmap(sub4, features=genes4, assay="SCT", slot="scale.data",
                group.by="celltype", group.bar=TRUE, draw.lines=TRUE, raster=TRUE) +
  scale_fill_gradientn(colors=c("#742374","black","#F9ED37")) +
  theme(axis.text.y = element_text(size=6))
save_pdf(p4, "04_heatmap_four_preBCR_states.pdf", 10, 12)

# 火山图：两张都 7x5.5
Idents(bcell) <- factor(bcell$celltype, levels = ordered_levels)
volcano_simple(
  bcell, "PreBCRd", "PreBCRi I", "04_volcano_PreBCRd_vs_PreBCRiI.pdf",
  assay="RNA", xlim=c(-3,2), ylim=c(0,220),
  labels_left=c("Mif","C1qbp","Srm","Ccnd2","Nrgn","Mt1","Vpreb1","Igll1"),
  labels_right=c("Ebf1","Il7r","Igkc","Vpreb3","H2afx","H2afv","Hist1h2ae"),
  title="PreBCRd vs PreBCRi I", w=7, h=5.5
)
volcano_simple(
  bcell, "PreBCRi I", c("PreBCRi II S phase","PreBCRi II G2/M phase"),
  "04_volcano_PreBCRiI_vs_PreBCRiII.pdf",
  assay="RNA", xlim=c(-3,1), ylim=c(0,260),
  labels_left=c("Top2a","Pclaf","Fbxo5","Pdia4","Birc5","Cdk1"),
  labels_right=c("Tmsb4x","Tmsb10","Vim","Actg1","Dnajc7","Dynll1","Cd74","Ifi27l2a"),
  title="PreBCRi I vs PreBCRi II", w=7, h=5.5
)

# -------------------------
# 6) Fig5：VDJ 阶段 & 晚期成熟 + 对应基因 Feature
# -------------------------
bcell$VDJ_stage <- "Other"
bcell$VDJ_stage[Cells(bcell) %in% WhichCells(bcell, idents="Lambda Pre B")] <- "Lambda Pre B"
bcell$VDJ_stage[Cells(bcell) %in% WhichCells(bcell, idents="Kappa Pre B")]  <- "Kappa Pre B"
bcell$VDJ_stage[Cells(bcell) %in% WhichCells(bcell, idents="Pro B VDJ")]    <- "Pro B VDJ"
bcell$VDJ_stage <- factor(bcell$VDJ_stage, levels=c("Other","Lambda Pre B","Kappa Pre B","Pro B VDJ"))
p_vdj <- DimPlot(bcell, reduction="umap", group.by="VDJ_stage",
                 cols=c("grey85","#3953A4","#C23B31","#2E7D32")) + theme_base +
  ggtitle("V(D)J Recombination")
save_pdf(p_vdj, "05_highlight_VDJ.pdf", 7.5, 6.2)

bcell$Mat_stage <- "Other"
bcell$Mat_stage[Cells(bcell) %in% WhichCells(bcell, idents="Immature B")]         <- "Immature B"
bcell$Mat_stage[Cells(bcell) %in% WhichCells(bcell, idents="Cycling Immature B")] <- "Cycling Immature B"
bcell$Mat_stage[Cells(bcell) %in% WhichCells(bcell, idents="Mature B")]           <- "Mature B"
bcell$Mat_stage <- factor(bcell$Mat_stage, levels=c("Other","Immature B","Cycling Immature B","Mature B"))
p_mat <- DimPlot(bcell, reduction="umap", group.by="Mat_stage",
                 cols=c("grey85","#2E7D32","#D33F6A","#3B4CC0")) + theme_base +
  ggtitle("Immature + Mature B cells")
save_pdf(p_mat, "05_highlight_maturation.pdf", 7.5, 6.2)

DefaultAssay(bcell) <- "SCT"
fp_grid(bcell, c("Rag1","Rag2","Igkc","Iglc1","Iglc2","Iglc3"),
        "05_feature_VDJ_genes.pdf", ncol = 3, nrow = 2, assay = "SCT",w = 12, h = 8)
fp_grid(bcell, c("Ms4a1","Ms4a4c","H2-Aa","Sell","Apoe","Ltb"),
        "05_feature_maturation_genes.pdf", ncol = 3, nrow = 2, assay = "SCT",w = 12, h = 8)
