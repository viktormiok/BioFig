library(ShrinkBayes)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)


source("functions_vm.R")
data(CAGEdata10000)
data(design_brain)


id <- colnames(CAGEdata10000)
pers <- design_brain$pers
batch <- design_brain$batch
groupfac <- design_brain$groupfac

# meta data
metaData <- data.frame(id, pers, batch, groupfac)

# DESeq object
dds <- DESeqDataSetFromMatrix(countData = CAGEdata10000,
                              colData = metaData,
                              design =~ batch
)


dds <- estimateSizeFactors(dds)

#####################################################################################################
#
#       PCA plot
#
####################################################################################################

pca_plot(dds = dds,
         condition = groupfac,
         intgroup = "groupfac",
         title = "PCA plot",
         point.size = 3,
         col.legend = 1
)

#####################################################################################################
#
#       sample dendogram with sample distance matrix heatmap
#
####################################################################################################

sample_distance(dds = dds,
                variables = c('pers', 'batch', 'groupfac'),
                colors = colorRampPalette(rev(brewer.pal(9, "Spectral")))(n=255)
)

#####################################################################################################
#
#       Volcano plot
#
####################################################################################################


dds <- DESeq(dds, fit='local')

res <- results(dds)#, contrast = c("condition","Input_CCK","Input_Veh"))
res = res[order(res$padj),]

results = as.data.frame(mutate(as.data.frame(res), 
                               sig = ifelse(res$padj<0.05, "FDR<0.05", "Not Sig")),
                        row.names=rownames(res)
)
results = results[(!is.na(results$padj)),]


volcano_plot(results = results,
             legend.name = significant,
             title = "Volcano plot",
             col = c("red", "black"),
             xlim = c(-10, 10)
)



             
  
# )
# ggplot(results, aes(log2FoldChange, -log10(padj))) +
#            geom_point(aes(col = sig)) + 
#            scale_color_manual(values = c("red", "black")) +
#            guides(col = guide_legend(nrow=2)) + 
#            ggtitle("Input_CCK_vs_Input_Veh") + 
#            xlim(-10,10)
# 
# p + geom_text_repel(data=results[1:20, ], aes(label=rownames(results[1:20, ])))
# 
