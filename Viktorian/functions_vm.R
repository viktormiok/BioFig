#####################################################################################################
#
#       PCA plot
#
####################################################################################################

pca_plot <- function(dds,
                     condition,
                     intgroup,
                     title,
                     point.size,
                     col.legend){
            pcaData <- plotPCA(vst(dds, fitType = "local"),
                               intgroup = intgroup,
                               returnData = TRUE
            )
            percentVar <- round(100 * attr(pcaData, "percentVar"))
            ggplot(pcaData, 
                   aes(PC1, PC2, color = condition)) +
                   geom_point(size = point.size) +
                   geom_label(aes(label = colnames(dds))) +
                   xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
                   ylab(paste0("PC2: ", percentVar[2], "% variance")) +
                   guides(color = guide_legend(ncol = col.legend)) +
                   ggtitle(title)
}


#####################################################################################################
#
#       sample dendogram with sample distance matrix heatmap
#
####################################################################################################

sample_distance <- function(dds,
                            variables,
                            colors){
  
                    vsd <- assay(vst(dds, 
                                     fitType="local",
                                     blind = TRUE)
                    )
                    sd_mat <- as.matrix(dist(t(vsd)))
                    
                    # Data frame with column annotations.
                    col_ann <- as.data.frame(colData(dds)[, variables])
                    rownames(col_ann) <- colnames(sd_mat)
                    
                    # List with colors for each annotation.
                    ann_colors <- list(condition = brewer.pal(12, "Set3"))
                    #names(ann_colors$condition) <- unique(colData(dds)$pers)
                    
                    
                    sdm <- pheatmap(sd_mat,
                                    clustering_distance_rows = dist(t(vsd)),
                                    clustering_distance_cols = dist(t(vsd)),
                                    show_colnames = TRUE,
                                    show_rownames = FALSE,
                                    annotation = col_ann,
                                    annotation_colors = ann_colors,
                                    col = colors
                    )
}

#####################################################################################################
#
#       Volcano plot
#
####################################################################################################

volcano_plot <- function(results,
                         title,
                         legend.name,
                         num_sig_gene,
                         col,
                         xlim, 
                         ylim){
                ggplot(results, aes(log2FoldChange, -log10(padj))) +
                       geom_point(aes(col = sig)) + 
                       scale_color_manual(values = col) +
                       guides(col = guide_legend(nrow=2)) + 
                       ggtitle(title) + 
                       xlim(xlim) +
                       geom_text_repel(data = results[1:20, ], 
                                       aes(label = rownames(results[1:20, ])))
  
}

