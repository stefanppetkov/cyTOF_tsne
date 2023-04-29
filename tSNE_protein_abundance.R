
# load('CD4_PR_RT.RData')
library(readxl)
metadata_filename <- "CD4_metadata_additional.xlsx"
md <- read_excel(metadata_filename)

## Make sure condition variables are factors with the right levels 
md$condition <- factor(md$condition, levels = c("RT", "PR")) 
head(data.frame(md))

## Define colors for conditions
color_conditions <- c("#FF7F00","#42f495") 
names(color_conditions) <- levels(md$condition)


## Load data
library(flowCore) 
fcs_raw <- read.flowSet(md$file_name, transformation = FALSE, truncate_max_range = FALSE) 
fcs_raw

## rename channel desc
for (x in 1:length(fcs_raw)) {
  pData(parameters(fcs_raw[[x]]))$desc[6] <- 'CD4'
  pData(parameters(fcs_raw[[x]]))$desc[7] <- 'CD8'
  pData(parameters(fcs_raw[[x]]))$desc[8] <- 'CCR6'
  pData(parameters(fcs_raw[[x]]))$desc[9] <- 'CCR10'
  pData(parameters(fcs_raw[[x]]))$desc[10] <- 'CD45RA'
  pData(parameters(fcs_raw[[x]]))$desc[11] <- 'CXCR5'
  pData(parameters(fcs_raw[[x]]))$desc[12] <- 'CXCR3'
  pData(parameters(fcs_raw[[x]]))$desc[13] <- 'CCR4'
  pData(parameters(fcs_raw[[x]]))$desc[14] <- 'VIABILITY'
}

panel_filename <- "SPLENOCYTE_panel.xlsx" 
panel <- read_excel(panel_filename) 
head(data.frame(panel))

# Replace problematic characters
panel$Antigen <- gsub("-", "_", panel$Antigen)

panel_fcs <- pData(parameters(fcs_raw[[1]]))
head(panel_fcs)

# Replace problematic characters
panel_fcs$desc <- gsub("-", "_", panel_fcs$desc)

# Lineage markers
(lineage_markers <- panel$Antigen[panel$Lineage == 1])

# Functional markers
(functional_markers <- panel$Antigen[panel$Functional == 1])

# Spot checks
all(lineage_markers %in% panel_fcs$desc)
all(functional_markers %in% panel_fcs$desc)

## arcsinh transformation and column subsetting
fcs <- fsApply(fcs_raw, function(x, cofactor=150){
  colnames(x) <- panel_fcs$desc
  expr <- exprs(x)
  expr <- asinh(expr[ , c(lineage_markers, functional_markers)] / cofactor)
  exprs(x) <- expr
  x
})
fcs

## Extract expression
expr <- fsApply(fcs, exprs)
dim(expr)

library(matrixStats)
rng <- colQuantiles(expr, probs = c(0.01, 0.99))
expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
expr01[expr01 < 0] <- 0
expr01[expr01 > 1] <- 1

## Generate sample IDs corresponding to each cell in the 'expr' matrix
sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))

library(ggplot2)
library(reshape2)

ggdf <- data.frame(sample_id = sample_ids, expr)
ggdf <- melt(ggdf, id.var = "sample_id",
             value.name = "expression", variable.name = "antigen")
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- md$condition[mm]
###########################################################################
ggplot(ggdf, aes(x = expression, color = condition,
                 group = sample_id)) +
  geom_density() +
  facet_wrap(~ antigen, nrow = 4, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(size = 7), axis.text = element_text(size = 5)) +
  scale_color_manual(values = color_conditions)
ggsave('expression.pdf', plot = last_plot(), device = 'pdf')
## Number of cells per sample
##########################################################################
cell_table <- table(sample_ids)

ggdf <- data.frame(sample_id = names(cell_table),
                   cell_counts = as.numeric(cell_table))
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- md$condition[mm]

ggplot(ggdf, aes(x = sample_id, y = cell_counts, fill = condition)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = cell_counts), hjust=0.5, vjust=-0.5, size = 2.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_manual(values = color_conditions, drop = FALSE) +
  scale_x_discrete(drop = FALSE)

ggsave('cell_number.pdf', plot = last_plot(), device = 'pdf')

library(dplyr)
# Get the median marker expression per sample

expr_median_sample_tbl <- data.frame(sample_id = sample_ids, expr) %>%
  group_by(sample_id) %>%  summarize_all(funs(median))

expr_median_sample <- t(expr_median_sample_tbl[, -1])
colnames(expr_median_sample) <- expr_median_sample_tbl$sample_id

library(limma)
mds <- plotMDS(expr_median_sample, plot = FALSE)

library(ggrepel)
ggdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y,
                   sample_id = colnames(expr_median_sample))
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- md$condition[mm]

ggplot(ggdf, aes(x = MDS1, y = MDS2, color = condition)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_label_repel(aes(label = sample_id)) +
  theme_bw() +
  scale_color_manual(values = color_conditions) +
  coord_fixed()
ggsave('MDA.pdf', plot = last_plot(), device = 'pdf')

##heatmap

library(RColorBrewer)
library(pheatmap)

# Column annotation for the heatmap
mm <- match(colnames(expr_median_sample), md$sample_id)
annotation_col <- data.frame(condition = md$condition[mm],
                             row.names = colnames(expr_median_sample))
annotation_colors <- list(condition = color_conditions)

# Colors for the heatmap
color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)

pheatmap(expr_median_sample, color = color, display_numbers = TRUE,
         number_color = "black", fontsize_number = 5, annotation_col = annotation_col,
         annotation_colors = annotation_colors, clustering_method = "average", filename = 'pheatmap.pdf')

## Define a function that calculates the NRS (non-redundancy scores, explaining the variability) per sample
NRS <- function(x, ncomp = 3){
  pr <- prcomp(x, center = TRUE, scale. = FALSE)
  score <- rowSums(outer(rep(1, ncol(x)),
                         pr$sdev[1:ncomp]^2) * abs(pr$rotation[,1:ncomp]))
  return(score)
}

## Calculate the score
nrs_sample <- fsApply(fcs[, lineage_markers], NRS, use.exprs = TRUE)
rownames(nrs_sample) <- md$sample_id
nrs <- colMeans(nrs_sample, na.rm = TRUE)

## Plot the NRS for ordered markers
lineage_markers_ord <- names(sort(nrs, decreasing = TRUE))
nrs_sample <- data.frame(nrs_sample)
nrs_sample$sample_id <- rownames(nrs_sample)

ggdf <- melt(nrs_sample, id.var = "sample_id",
             value.name = "nrs", variable.name = "antigen")

ggdf$antigen <- factor(ggdf$antigen, levels = lineage_markers_ord)
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- md$condition[mm]

ggplot(ggdf, aes(x = antigen, y = nrs)) +
  geom_point(aes(color = condition), alpha = 0.9,
             position = position_jitter(width = 0.3, height = 0)) +
  geom_boxplot(outlier.color = NA, fill = NA) +
  stat_summary(fun.y = "mean", geom = "point", shape = 21, fill = "white") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(values = color_conditions)
ggsave('NRS.pdf', plot = last_plot(), device = 'pdf')

## SOM clustering!!!
library(FlowSOM)
fsom <- ReadInput(fcs, transform = FALSE, scale = FALSE)
set.seed(1234)
som <- BuildSOM(fsom, colsToUse = lineage_markers)
## Get the cell clustering into 100 SOM codes
cell_clustering_som <- som$map$mapping[,1]

## METACLUSTERING
## Metaclustering into 20 clusters with ConsensusClusterPlus
library(ConsensusClusterPlus)

codes <- som$map$codes
plot_outdir <- "consensus_plots"
nmc <- 20

mc <- ConsensusClusterPlus(t(codes), maxK = nmc, reps = 100,
                           pItem = 0.9, pFeature = 1, title = plot_outdir, plot = "png",
                           clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average",
                           distance = "euclidean", seed = 1234)

## Get cluster ids for each cell
code_clustering1 <- mc[[nmc]]$consensusClass
cell_clustering1 <- code_clustering1[cell_clustering_som]

# Define cluster colors (here there are 30 colors)
color_clusters <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

plot_clustering_heatmap_wrapper <- function(expr, expr01,
                                            cell_clustering, color_clusters, cluster_merging = NULL){
  
  # Calculate the median expression
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  expr01_median <- data.frame(expr01, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  
  # Calculate cluster frequencies
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  # Sort the cell clusters with hierarchical clustering
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_median[, colnames(expr01)])
  rownames(expr_heat) <- expr01_median$cell_clustering
  
  # Colors for the heatmap
  color_heat <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
  legend_breaks = seq(from = 0, to = 1, by = 0.2)
  labels_row <- paste0(expr01_median$cell_clustering, " (", clustering_prop ,
                       "%)")
  
  # Annotation for the original clusters
  annotation_row <- data.frame(Cluster = factor(expr01_median$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  # Annotation for the merged clusters
  if(!is.null(cluster_merging)){
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$Cluster_merging <- cluster_merging$new_cluster
    color_clusters2 <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters2) <- levels(cluster_merging$new_cluster)
    annotation_colors$Cluster_merging <- color_clusters2
  }
  
  pheatmap(expr_heat, color = color_heat, cluster_cols = FALSE,
           cluster_rows = cluster_rows, labels_row = labels_row,
           display_numbers = TRUE, number_color = "black",
           fontsize = 8, fontsize_number = 6,  legend_breaks = legend_breaks,
           annotation_row = annotation_row, annotation_colors = annotation_colors)
  
}
######################################################################################
## HEATMAP: 100 SOM CLUSTERS INTO 20 METACLUSTERS
dev.off(which = dev.cur())
plot_clustering_heatmap_wrapper(expr = expr[, lineage_markers_ord],
                                expr01 = expr01[, lineage_markers_ord],
                                cell_clustering = cell_clustering1, color_clusters = color_clusters)
######################################################################################

library(ggridges)
plot_clustering_distr_wrapper <- function(expr, cell_clustering){
  # Calculate the median expression
  cell_clustering <- factor(cell_clustering)
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  # Sort the cell clusters with hierarchical clustering
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  # Calculate cluster frequencies
  freq_clust <- table(cell_clustering)
  freq_clust <- round(as.numeric(freq_clust)/sum(freq_clust)*100, 2)
  cell_clustering <- factor(cell_clustering,
                            labels = paste0(levels(cell_clustering), "  (", freq_clust, "%)"))
  ### Data organized per cluster
  ggd <- melt(data.frame(cluster = cell_clustering, expr),
              id.vars = "cluster", value.name = "expression",
              variable.name = "antigen")
  ggd$antigen <- factor(ggd$antigen, levels = colnames(expr))
  ggd$reference <- "no"
  ### The reference data
  ggd_bg <- ggd
  ggd_bg$cluster <- "reference"
  ggd_bg$reference <- "yes"
  
  ggd_plot <- rbind(ggd, ggd_bg)
  ggd_plot$cluster <- factor(ggd_plot$cluster,
                             levels = c(levels(cell_clustering)[rev(cluster_rows$order)], "reference"))
  
  ggplot() +
    geom_density_ridges(data = ggd_plot, aes(x = expression, y = cluster,
                                             color = reference, fill = reference), alpha = 0.3) +
    facet_wrap( ~ antigen, scales = "free_x", nrow = 2) +
    theme_ridges() +
    theme(axis.text = element_text(size = 7),  
          strip.text = element_text(size = 7), legend.position = "none")
  
}
plot_clustering_distr_wrapper(expr = expr[, lineage_markers_ord],
                              cell_clustering = cell_clustering1)

library(ComplexHeatmap)

plot_clustering_heatmap_wrapper2 <- function(expr, expr01,
                                             lineage_markers, functional_markers = NULL, sample_ids = NULL,
                                             cell_clustering, color_clusters, cluster_merging = NULL,
                                             plot_cluster_annotation = TRUE){
  
  # Calculate the median expression of lineage markers
  expr_median <- data.frame(expr[, lineage_markers],
                            cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  expr01_median <- data.frame(expr01[, lineage_markers],
                              cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  
  # Calculate cluster frequencies
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  # Sort the cell clusters with hierarchical clustering
  d <- dist(expr_median[, lineage_markers], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_median[, lineage_markers])
  
  # Median expression of functional markers in each sample per cluster
  expr_median_sample_cluster_tbl <- data.frame(expr01[, functional_markers,
                                                      drop = FALSE], sample_id = sample_ids, cluster = cell_clustering) %>%
    group_by(sample_id, cluster) %>% summarize_all(funs(median))
  
  # Colors for the heatmap
  color_heat <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
  legend_breaks = seq(from = 0, to = 1, by = 0.2)
  labels_row <- paste0(expr01_median$cell_clustering, " (", clustering_prop , "%)") 
  
  ### Annotation for the original clusters
  annotation_row1 <- data.frame(Cluster = factor(expr01_median$cell_clustering))
  color_clusters1 <- color_clusters[1:nlevels(annotation_row1$Cluster)]
  names(color_clusters1) <- levels(annotation_row1$Cluster)
  
  
  ### Annotation for the merged clusters
  if(!is.null(cluster_merging)){
    mm <- match(annotation_row1$Cluster, cluster_merging$original_cluster)
    annotation_row2 <- data.frame(Cluster_merging =
                                    factor(cluster_merging$new_cluster[mm]))
    color_clusters2 <- color_clusters[1:nlevels(annotation_row2$Cluster_merging)]
    names(color_clusters2) <- levels(annotation_row2$Cluster_merging)
  }
  
  
  ### Heatmap annotation for the original clusters
  ha1 <- Heatmap(annotation_row1, name = "Cluster",
                 col = color_clusters1, cluster_columns = FALSE,
                 cluster_rows = cluster_rows, row_dend_reorder = FALSE,
                 show_row_names = FALSE, width = unit(0.5, "cm"),
                 rect_gp = gpar(col = "grey"))
  ### Heatmap annotation for the merged clusters
  if(!is.null(cluster_merging)){
    ha2 <- Heatmap(annotation_row2, name = "Cluster \nmerging",
                   col = color_clusters2, cluster_columns = FALSE,
                   cluster_rows = cluster_rows, row_dend_reorder = FALSE,
                   show_row_names = FALSE, width = unit(0.5, "cm"),
                   rect_gp = gpar(col = "grey"))
  }
  ### Cluster names and sizes - text
  ha_text <- rowAnnotation(text = row_anno_text(labels_row,
                                                gp = gpar(fontsize = 6)), width = max_text_width(labels_row))
  ### Cluster sizes - barplot
  ha_bar <- rowAnnotation("Frequency (%)" = row_anno_barplot(
    x = clustering_prop, border = FALSE, axis = TRUE,
    axis_gp = gpar(fontsize = 5), gp = gpar(fill = "#696969", col = "#696969"),
    bar_width = 0.9), width = unit(0.7, "cm"), show_annotation_name = TRUE,
    annotation_name_rot = 0, annotation_name_offset = unit(5, "mm"),
    annotation_name_gp = gpar(fontsize = 7))
  ### Heatmap for the lineage markers
  ht1 <- Heatmap(expr_heat, name = "Expr", column_title = "Lineage markers",
                 col = color_heat, cluster_columns = FALSE, cluster_rows = cluster_rows,
                 row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks,
                                                                       labels = legend_breaks, color_bar = "continuous"),
                 show_row_names = FALSE, row_dend_width = unit(2, "cm"),
                 rect_gp = gpar(col = "grey"), column_names_gp = gpar(fontsize = 8))
  
  if(plot_cluster_annotation){
    draw_out <- ha1
  }else{
    draw_out <- NULL
  }
  if(!is.null(cluster_merging)){
    draw_out <- draw_out + ha2 + ht1 + ha_bar + ha_text
  }else{
    draw_out <- draw_out + ht1 + ha_bar + ha_text
  }
  ### Heatmaps for the signaling markers
  if(!is.null(functional_markers)){
    for(i in 1:length(functional_markers)){
      ## Rearange so the rows represent clusters
      expr_heat_fun <- as.matrix(dcast(expr_median_sample_cluster_tbl[,
                                                                      c("sample_id", "cluster", functional_markers[i])],
                                       cluster ~ sample_id, value.var = functional_markers[i])[, -1])
      draw_out <- draw_out + Heatmap(expr_heat_fun,
                                     column_title = functional_markers[i], col = color_heat,
                                     cluster_columns = FALSE, cluster_rows = cluster_rows,
                                     row_dend_reorder = FALSE, show_heatmap_legend = FALSE,
                                     show_row_names = FALSE, rect_gp = gpar(col = "grey"),
                                     column_names_gp = gpar(fontsize = 8))
    }
  }
  draw(draw_out, row_dend_side = "left")
}
#######################################################################################

plot_clustering_heatmap_wrapper2(expr = expr, expr01 = expr01,
                                 lineage_markers = lineage_markers, functional_markers = "CD8",
                                 sample_ids = sample_ids, cell_clustering = cell_clustering1,
                                 color_clusters = color_clusters, cluster_merging = NULL)
#######################################################################################

## Find and skip duplicates
dups <- which(!duplicated(expr[, lineage_markers]))

## Data subsampling: create indices by sample
inds <- split(1:length(sample_ids), sample_ids)

## How many cells to downsample per-sample
tsne_ncells <- pmin(table(sample_ids), 5000)

## Get subsampled indices
set.seed(1234)
tsne_inds <- lapply(names(inds), function(i){
  s <- sample(inds[[i]], tsne_ncells[i], replace = FALSE)
  intersect(s, dups)
})

tsne_inds <- unlist(tsne_inds)

tsne_expr <- expr[tsne_inds, lineage_markers]
## Run t-SNE
library(Rtsne)

set.seed(1234)
tsne_out <- Rtsne(tsne_expr, check_duplicates = FALSE, pca = FALSE)

## Plot t-SNE colored by CD4 expression
dr <- data.frame(tSNE1 = tsne_out$Y[, 1], tSNE2 = tsne_out$Y[, 2],
                 expr[tsne_inds, lineage_markers])
######################################################################################
ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = CXCR3)) +
  geom_point(size = 0.8) +
  theme_bw() +
  scale_color_gradientn("CD8",
                        colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))
ggsave('tSNE_CXCR3.pdf', plot = last_plot(), device = 'pdf')
dev.off(which = dev.cur())
######################################################################################

dr$sample_id <- sample_ids[tsne_inds]
mm <- match(dr$sample_id, md$sample_id)
dr$condition <- md$condition[mm]
dr$cell_clustering1 <- factor(cell_clustering1[tsne_inds], levels = 1:nmc)

## Plot t-SNE colored by clusters
ggp <- ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = cell_clustering1)) +
  geom_point(size = 0.8) +
  theme_bw() +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
######################################################################################
ggp
ggsave('tSNE_20_clusters.pdf', plot = last_plot(), device = 'pdf')
## Facet per sample
ggp + facet_wrap(~ sample_id)
ggsave('tSNE_20_clusters_by_sample.pdf', plot = last_plot(), device = 'pdf')
## Facet per condition
ggp + facet_wrap(~ condition)
ggsave('tSNE_20_clusters_by_condition.pdf', plot = last_plot(), width= 12, height = 4,
       device = 'pdf')
dev.off(which = dev.cur())
######################################################################################
# Density plots faceter per condition
ggplot(dr,  aes(x = tSNE1, y = tSNE2)) +  
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_distiller(palette = "Spectral") +
  theme_bw() +
  facet_wrap(~ condition)
ggsave('tSNE_20_clusters_by_condition_density_polygon.pdf', plot = last_plot(), 
       width= 12, height = 4, device = 'pdf')
#or

ggplot(dr,  aes(x = tSNE1, y = tSNE2)) +  
  stat_density_2d(aes(color = ..level..), geom = "path") +
  scale_color_distiller(palette = "Spectral") +
  theme_bw() +
  facet_wrap(~ condition)
ggsave('tSNE_20_clusters_by_condition_density_path.pdf', plot = last_plot(), 
       width= 12, height = 4, device = 'pdf')
dev.off(which = dev.cur())
######################################################################################


# --------------------------------END----------------------------------------
# ## Get code sizes; sometimes not all the codes have mapped cells so they will have size 0
# code_sizes <- table(factor(som$map$mapping[, 1], levels = 1:nrow(codes)))
# code_sizes <- as.numeric(code_sizes)
# ## Run t-SNE on codes
# set.seed(1234)
# tsne_out <- Rtsne(codes, perplexity = 5, pca = FALSE)
# ## Run PCA on codes
# pca_out <- prcomp(codes, center = TRUE, scale. = FALSE)
# codes_dr <- data.frame(tSNE1 = tsne_out$Y[, 1], tSNE2 = tsne_out$Y[, 2],
#                        PCA1 = pca_out$x[, 1], PCA2 = pca_out$x[, 2])
# codes_dr$code_clustering1 <- factor(code_clustering1)
# codes_dr$size <- code_sizes
# 
# ## Plot t-SNE on codes
# gg_tsne_codes <- ggplot(codes_dr,  aes(x = tSNE1, y = tSNE2,
#                                        color = code_clustering1, size = size)) +
#   geom_point(alpha = 0.9) +
#   theme_bw() +
#   scale_color_manual(values = color_clusters) +
#   guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
# ## Plot PCA on codes
# gg_pca_codes <- ggplot(codes_dr,  aes(x = PCA1, y = PCA2,
#                                       color = code_clustering1, size = size)) +
#   geom_point(alpha = 0.9) +
#   theme_bw() +
#   scale_color_manual(values = color_clusters) +
#   guides(color = guide_legend(override.aes = list(size = 4), ncol = 2)) +  
#   theme(legend.position = "right", legend.box = "vertical")
# 
# library(cowplot)
# 
# legend <- get_legend(gg_tsne_codes)
# ggdraw() +
#   draw_plot(gg_tsne_codes + theme(legend.position = "none"), 0, .5, .7, .5) +
#   draw_plot(gg_pca_codes + theme(legend.position = "none"), 0, 0, .7, .5) +
#   draw_plot(legend, .7, .0, .2, 1) +
#   draw_plot_label(c("A", "B", ""), c(0, 0, .7), c(1, .5, 1), size = 15)
# 
# ## SOM clusters
# plot_clustering_heatmap_wrapper2(expr = expr, expr01 = expr01,
#                                  lineage_markers = lineage_markers, functional_markers = "CD8",
#                                  sample_ids = sample_ids, cell_clustering = cell_clustering_som,
#                                  color_clusters = rep(color_clusters, length.out = 100),
#                                  cluster_merging = data.frame(original_cluster = 1:100,
#                                                               new_cluster = code_clustering1), plot_cluster_annotation = FALSE)
# 
# 
# ----------------BEGIN-----------------------
### CLUSTER MERGING
cluster_merging1_filename <- "CD4_cluster_merging.xlsx"
# download.file(paste0(url, "/", cluster_merging1_filename),
#               destfile = cluster_merging1_filename, mode = "wb")
cluster_merging1 <- read_excel(cluster_merging1_filename)
data.frame(cluster_merging1)

## Convert to factor with merged clusters in desired order
levels_clusters_merged <- c("Th1", "Th2", "Th9", "Tfh", "CCR10+", "undefined")
cluster_merging1$new_cluster <- factor(cluster_merging1$new_cluster,
                                       levels = levels_clusters_merged)
## MANUAL CLUSTERS
## New clustering1m
mm <- match(cell_clustering1, cluster_merging1$original_cluster)
cell_clustering1m <- cluster_merging1$new_cluster[mm]

mm <- match(code_clustering1, cluster_merging1$original_cluster)
code_clustering1m <- cluster_merging1$new_cluster[mm]

dr$cell_clustering1m <- cell_clustering1m[tsne_inds]

## PLOT WITH NAMED CLUSTERS
ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = cell_clustering1m)) +
  geom_point(size = 0.8) +
  theme_bw() +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4)))
ggsave('tSNE_named_clusters.pdf', plot = last_plot(), device = 'pdf')

## Named clusters by condition
ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = cell_clustering1m)) +
  geom_point(size = 0.8) +
  theme_bw() +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  facet_wrap(~ condition)
ggsave('tSNE_named_clusters_by_condition.pdf', plot = last_plot(),
       width = 12, height = 4, device = 'pdf')

######################################################################################
ggplot(dr,  aes(x = tSNE1, y = tSNE2)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_distiller(palette = "Spectral") +
  theme_bw() +
  facet_wrap(~ sample_id)
ggsave('tSNE_named_clusters_by_sample.pdf', plot = last_plot(),
       width = 12, height = 4, device = 'pdf')
dev.off(which = dev.cur())

## COMPARISON BETWEEN 20 CLUSTER METACLUSTERING AND MERGING
######################################################################################
plot_clustering_heatmap_wrapper(expr = expr[, lineage_markers_ord],
                                expr01 = expr01[, lineage_markers_ord], cell_clustering = cell_clustering1,
                                color_clusters = color_clusters, cluster_merging = cluster_merging1)

plot_clustering_heatmap_wrapper(expr = expr[, lineage_markers_ord],
                                expr01 = expr01[, lineage_markers_ord], cell_clustering = cell_clustering1m,
                                color_clusters = color_clusters)
######################################################################################

## Get cluster ids for each cell
nmc2 <- 10
code_clustering2 <- mc[[nmc2]]$consensusClass
cell_clustering2 <- code_clustering2[som$map$mapping[, 1]]

dr$cell_clustering2 <- factor(cell_clustering2[tsne_inds], levels = 1:nmc2)
ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = cell_clustering2)) +
  geom_point(size = 0.8) +
  theme_bw() +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
plot_clustering_heatmap_wrapper(expr = expr[, lineage_markers_ord],
                                expr01 = expr01[, lineage_markers_ord], cell_clustering = cell_clustering2,
                                color_clusters = color_clusters)

cluster_merging2_filename <- "CD4_cluster_merging2.xlsx"
# download.file(paste0(url, "/", cluster_merging2_filename),
#               destfile = cluster_merging2_filename, mode = "wb")
cluster_merging2 <- read_excel(cluster_merging2_filename)
data.frame(cluster_merging2)

## Convert to factor with merged clusters in correct order
cluster_merging2$new_cluster <- factor(cluster_merging2$new_cluster,
                                       levels = levels_clusters_merged)
## New clustering2m
mm <- match(cell_clustering2, cluster_merging2$original_cluster)
cell_clustering2m <- cluster_merging2$new_cluster[mm]
dr$cell_clustering2m <- cell_clustering2m[tsne_inds]
gg_tsne_cl2m <- ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = cell_clustering2m)) +
  geom_point(size = 0.8) +
  theme_bw() +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4)))
gg_tsne_cl2m
plot_clustering_heatmap_wrapper(expr = expr[, lineage_markers_ord],
                                expr01 = expr01[, lineage_markers_ord], cell_clustering = cell_clustering2,
                                color_clusters = color_clusters, cluster_merging = cluster_merging2)
plot_clustering_heatmap_wrapper(expr = expr[, lineage_markers_ord],
                                expr01 = expr01[, lineage_markers_ord], cell_clustering = cell_clustering2m,
                                color_clusters = color_clusters)

## Get cluster ids for each cell
nmc3 <- 8
code_clustering3 <- mc[[nmc3]]$consensusClass
cell_clustering3 <- code_clustering3[som$map$mapping[, 1]]
# tabular comparison of cell_clustering3 and cell_clustering2m
table(algorithm=cell_clustering3, manual=cell_clustering2m)

dr$cell_clustering3 <- factor(cell_clustering3[tsne_inds], levels = 1:nmc3)
gg_tsne_cl3 <- ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = cell_clustering3)) +
  geom_point(size = 0.8) +
  theme_bw() +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4)))

plot_grid(gg_tsne_cl2m, gg_tsne_cl3, ncol = 1, labels = c("A", "B"))
plot_clustering_heatmap_wrapper(expr = expr[, lineage_markers_ord],
                                expr01 = expr01[, lineage_markers_ord], cell_clustering = cell_clustering3,
                                color_clusters = color_clusters)


## DIFFERENTIAL ABUNDANCE ANALYSIS
library(lme4)
library(multcomp)
## Model formula without random effects
model.matrix(~ condition, data = md)

## Create contrasts
contrast_names <- c("RTvsPR")
k1 <- c(0, 1)
K <- matrix(k1, nrow = 1, byrow = TRUE, dimnames = list(contrast_names))
K

FDR_cutoff <- 0.05

## Differential cell population abundance

counts_table <- table(cell_clustering1m, sample_ids)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)

ggdf <- melt(data.frame(cluster = rownames(props), props),
             id.vars = "cluster", value.name = "proportion", variable.name = "sample_id")
ggdf$cluster <- factor(ggdf$cluster, levels = levels_clusters_merged)
## Add condition info
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- factor(md$condition[mm])

##############################################
ggplot(ggdf, aes(x = sample_id, y = proportion, fill = cluster)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ condition, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = color_clusters)
ggsave('cluster_proportions.pdf', plot = last_plot(), device = 'pdf')
dev.off(which = dev.cur())
#############################################
ggdf$patient_id <- factor(md$patient_id[mm])

###########################################
ggplot(ggdf) +
  geom_boxplot(aes(x = condition, y = proportion, color = condition,
                   fill = condition), position = position_dodge(), alpha = 0.5,
               outlier.color = NA) +
  geom_point(aes(x = condition, y = proportion, color = condition,
                 shape = patient_id), alpha = 0.8, position = position_jitterdodge()) +
  facet_wrap(~ cluster, scales = "free", nrow = 2) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(), strip.text = element_text(size = 6)) +
  scale_color_manual(values = color_conditions) +
  scale_fill_manual(values = color_conditions) +
  scale_shape_manual(values = c(16, 17, 8, 3, 12, 0, 1, 2, 4, 5, 9, 10, 11, 7))
ggsave('cluster_boxplots.pdf', plot = last_plot(), device = 'pdf')
dev.off(which = dev.cur())
############################################
formula_glmer_binomial1 <- y/total ~ condition + (1|sample_id)
formula_glmer_binomial2 <- y/total ~ condition + (1|patient_id) + (1|sample_id)

differential_abundance_wrapper <- function(counts, md, formula, K){
  ## Fit the GLMM for each cluster separately
  ntot <- colSums(counts)
  fit_binomial <- lapply(1:nrow(counts), function(i){
    
    data_tmp <- data.frame(y = as.numeric(counts[i, md$sample_id]),
                           total = ntot[md$sample_id], md)
    
    fit_tmp <- glmer(formula, weights = total, family = binomial,
                     data = data_tmp)
    
    ## Fit contrasts one by one
    out <- apply(K, 1, function(k){
      contr_tmp <- glht(fit_tmp, linfct = matrix(k, 1))
      summ_tmp <- summary(contr_tmp)
      pval <- summ_tmp$test$pvalues
      return(pval)
    })
    return(out)
  })
  pvals <- do.call(rbind, fit_binomial)
  colnames(pvals) <- paste0("pval_", contrast_names)
  rownames(pvals) <- rownames(counts)
  ## Adjust the p-values
  adjp <- apply(pvals, 2, p.adjust, method = "BH")
  colnames(adjp) <- paste0("adjp_", contrast_names)
  return(list(pvals = pvals, adjp = adjp))
}

da_out1 <- differential_abundance_wrapper(counts, md = md,
                                          formula = formula_glmer_binomial1, K = K)
apply(da_out1$adjp < FDR_cutoff, 2, table)

da_out2 <- differential_abundance_wrapper(counts, md = md,
                                          formula = formula_glmer_binomial2, K = K)
apply(da_out2$adjp < FDR_cutoff, 2, table)

da_output2 <- data.frame(cluster = rownames(props), props,
                         da_out2$pvals, da_out2$adjp, row.names = NULL)
print(head(da_output2), digits = 2)
write.csv(da_output2, file = 'stats.csv')

normalization_wrapper <- function(expr, th = 2.5){
  expr_norm <- apply(expr, 1, function(x){
    sdx <- sd(x, na.rm = TRUE)
    if(sdx == 0){
      x <- (x - mean(x, na.rm = TRUE))
    }else{
      x <- (x - mean(x, na.rm = TRUE)) / sdx
    }
    x[x > th] <- th
    x[x < -th] <- -th
    return(x)
  })
  expr_norm <- t(expr_norm)
}
plot_differential_heatmap_wrapper <- function(expr_norm, sign_adjp,
                                              condition, color_conditions, th = 2.5){
  ## Order samples by condition
  oo <- order(condition)
  condition <- condition[oo]
  expr_norm <- expr_norm[, oo, drop = FALSE]

  ## Create the row labels with adj p-values and other objects for pheatmap
  labels_row <- paste0(rownames(expr_norm), " (",
                       sprintf( "%.02e", sign_adjp), ")")
  labels_col <- colnames(expr_norm)
  annotation_col <- data.frame(condition = factor(condition))
  rownames(annotation_col) <- colnames(expr_norm)
  annotation_colors <- list(condition = color_conditions)
  color <- colorRampPalette(c("#87CEFA", "#56B4E9", "#56B4E9", "#0072B2",
                              "#000000", "#D55E00", "#E69F00", "#E69F00", "#FFD700"))(100)
  breaks = seq(from = -th, to = th, length.out = 101)
  legend_breaks = seq(from = -round(th), to = round(th), by = 1)
  gaps_col <- as.numeric(table(annotation_col$condition))

  pheatmap(expr_norm, color = color, breaks = breaks,
           legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE,
           labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col,
           annotation_col = annotation_col, annotation_colors = annotation_colors,
           fontsize = 8)
}
## Apply the arcsine-square-root transformation to the proportions
asin_table <- asin(sqrt((t(t(counts_table) / colSums(counts_table)))))
asin <- as.data.frame.matrix(asin_table)
## Get significant clusters and sort them by significance
sign_clusters <- names(which(sort(da_out2$adjp[, "adjp_PRvspVax"]) < FDR_cutoff))
## Get the adjusted p-values for the significant clusters
sign_adjp <- da_out2$adjp[sign_clusters , "adjp_PRvspVax", drop=FALSE]
## Normalize the transformed proportions to mean = 0 and sd = 1
asin_norm <- normalization_wrapper(asin[sign_clusters, ])
mm <- match(colnames(asin_norm), md$sample_id)
plot_differential_heatmap_wrapper(expr_norm = asin_norm, sign_adjp = sign_adjp,
                                  condition = md$condition[mm], color_conditions = color_conditions)

## Get median marker expression per sample and cluster
expr_median_sample_cluster_tbl <- data.frame(expr[, functional_markers],
                                             sample_id = sample_ids, cluster = cell_clustering1m) %>%
  group_by(sample_id, cluster) %>%
  summarize_all(funs(median))
## Melt
expr_median_sample_cluster_melt <- melt(expr_median_sample_cluster_tbl,
                                        id.vars = c("sample_id", "cluster"), value.name = "median_expression",
                                        variable.name = "antigen")
## Rearange so the rows represent clusters and markers
expr_median_sample_cluster <- dcast(expr_median_sample_cluster_melt,
                                    cluster + antigen ~ sample_id,  value.var = "median_expression")
rownames(expr_median_sample_cluster) <- paste0(expr_median_sample_cluster$cluster,
                                               "_", expr_median_sample_cluster$antigen)
## Eliminate clusters with low frequency
clusters_keep <- names(which((rowSums(counts < 5) == 0)))
keepLF <- expr_median_sample_cluster$cluster %in% clusters_keep
expr_median_sample_cluster <- expr_median_sample_cluster[keepLF, ]
## Eliminate cases with zero expression in all samples
keep0 <- rowSums(expr_median_sample_cluster[, md$sample_id]) > 0
expr_median_sample_cluster <- expr_median_sample_cluster[keep0, ]

ggdf <- expr_median_sample_cluster_melt[expr_median_sample_cluster_melt$cluster
                                        %in% clusters_keep, ]
## Add info about samples
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- factor(md$condition[mm])
ggdf$patient_id <- factor(md$patient_id[mm])
ggplot(ggdf) +
  geom_boxplot(aes(x = antigen, y = median_expression,
                   color = condition, fill = condition),
               position = position_dodge(), alpha = 0.5, outlier.color = NA) +
  geom_point(aes(x = antigen, y = median_expression, color = condition,
                 shape = patient_id), alpha = 0.8, position = position_jitterdodge(),
             size = 0.7) +
  facet_wrap(~ cluster, scales = "free_y", ncol=2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(values = color_conditions) +
  scale_fill_manual(values = color_conditions) +
  scale_shape_manual(values = c(16, 17, 8, 3, 12, 0, 1, 2)) +
  guides(shape = guide_legend(override.aes = list(size = 2)))
## No need to run further



differential_expression_wrapper <- function(expr_median, md, model = "lmer", formula, K){
  ## Fit LMM or LM for each marker separately
  fit_gaussian <- lapply(1:nrow(expr_median), function(i){
    data_tmp <- data.frame(y = as.numeric(expr_median[i, md$sample_id]), md)
    switch(model,
           lmer = {
             fit_tmp <- lmer(formula, data = data_tmp)
           },
           lm = {
             fit_tmp <- lm(formula, data = data_tmp)
           })
    ## Fit contrasts one by one
    out <- apply(K, 1, function(k){
      contr_tmp <- glht(fit_tmp, linfct = matrix(k, 1))
      summ_tmp <- summary(contr_tmp)
      pval <- summ_tmp$test$pvalues
      return(pval)
    })
    return(out)
  })
  pvals <- do.call(rbind, fit_gaussian)
  colnames(pvals) <- paste0("pval_", contrast_names)
  rownames(pvals) <- rownames(expr_median)
  ## Adjust the p-values
  adjp <- apply(pvals, 2, p.adjust, method = "BH")
  colnames(adjp) <- paste0("adjp_", contrast_names)
  return(list(pvals = pvals, adjp = adjp))
}

formula_lm <- y ~ condition
formula_lmer <- y ~ condition + (1|patient_id)

de_out1 <- differential_expression_wrapper(expr_median = expr_median_sample_cluster,
                                           md = md, model = "lm", formula = formula_lm, K = K)
apply(de_out1$adjp < FDR_cutoff, 2, table)

de_out2 <- differential_expression_wrapper(expr_median = expr_median_sample_cluster,
                                           md = md, model = "lmer", formula = formula_lmer, K = K)
apply(de_out2$adjp < FDR_cutoff, 2, table)

de_output2 <- data.frame(expr_median_sample_cluster,
                         de_out2$pvals, de_out2$adjp, row.names = NULL)
print(head(de_output2), digits = 2)

## Keep the significant markers, sort them by significance and group by cluster
sign_clusters_markers <- names(which(de_out2$adjp[, "adjp_PRvspVax"] < FDR_cutoff))
oo <- order(expr_median_sample_cluster[sign_clusters_markers, "cluster"],
            de_out2$adjp[sign_clusters_markers, "adjp_PRvspVax"])
sign_clusters_markers <- sign_clusters_markers[oo]
## Get the significant adjusted p-values
sign_adjp <- de_out2$adjp[sign_clusters_markers , "adjp_PRvspVax"]
## Normalize expression to mean = 0 and sd = 1
expr_s <- expr_median_sample_cluster[sign_clusters_markers,md$sample_id]
expr_median_sample_cluster_norm <- normalization_wrapper(expr_s)
mm <- match(colnames(expr_median_sample_cluster_norm), md$sample_id)
plot_differential_heatmap_wrapper(expr_norm = expr_median_sample_cluster_norm,
                                  sign_adjp = sign_adjp, condition = md$condition[mm],
                                  color_conditions = color_conditions)

ggdf <- melt(data.frame(expr_median_sample[functional_markers, ],
                        antigen = functional_markers), id.vars = "antigen",
             value.name = "median_expression", variable.name = "sample_id")
## Add condition info
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- factor(md$condition[mm])
ggdf$patient_id <- factor(md$patient_id[mm])
ggplot(ggdf) +
  geom_boxplot(aes(x = condition, y = median_expression, color = condition,
                   fill = condition),  position = position_dodge(), alpha = 0.5,
               outlier.color = NA) +
  geom_point(aes(x = condition, y = median_expression, color = condition,
                 shape = patient_id), alpha = 0.8, position = position_jitterdodge()) +
  facet_wrap(~ antigen, scales = "free", nrow = 5) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_color_manual(values = color_conditions) +
  scale_fill_manual(values = color_conditions) +
  scale_shape_manual(values = c(16, 17, 8, 3, 12, 0, 1, 2))

## Fit a linear model
de_out3 <- differential_expression_wrapper(expr_median =
                                             expr_median_sample[functional_markers, ],
                                           md = md, model = "lm", formula = formula_lm, K = K)
apply(de_out3$adjp < FDR_cutoff, 2, table)

## Fit a linear mixed model with patient ID as a random effect
de_out4 <- differential_expression_wrapper(expr_median =
                                             expr_median_sample[functional_markers, ],
                                           md = md, model = "lmer", formula = formula_lmer, K = K)
apply(de_out4$adjp < FDR_cutoff, 2, table)

de_output4 <- data.frame(antigen = functional_markers,
                         expr_median_sample[functional_markers, ], de_out4$pvals, de_out4$adjp)
print(head(de_output4), digits=2)

## Keep the significant markers and sort them by significance
sign_markers <- names(which(sort(de_out4$adjp[, "adjp_PRvspVax"]) < FDR_cutoff))
## Get the adjusted p-values
sign_adjp <- de_out4$adjp[sign_markers , "adjp_PRvspVax"]
## Normalize expression to mean = 0 and sd = 1
expr_median_sample_norm <- normalization_wrapper(expr_median_sample[sign_markers, ])

mm <- match(colnames(expr_median_sample_norm), md$sample_id)
plot_differential_heatmap_wrapper(expr_norm = expr_median_sample_norm,
                                  sign_adjp = sign_adjp, condition = md$condition[mm],
                                  color_conditions = color_conditions)

