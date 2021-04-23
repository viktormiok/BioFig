library(dplyr)
library(factoextra)
library(FactoMineR)
library(ggplot2)
library(ggrepel)
library(dendextend)


## PC plots 
plot_PCA <-function(exp.data,
                    group = NULL,
                    colors=NULL,
                    shape=NULL,
                    label =NULL,..){
  df_pca <-prcomp(t(exp.data))
  percentage <- round(df_pca$sdev^2 / sum(df_pca$sdev^2)*100,1)
  ## screeplots 
  screeplot <- fviz_eig(df_pca, barcolor = "grey", 
           linecolor = "red", 
           addlabels = TRUE, 
           hjust = 0, 
           ncp=10) +
    labs(title = "PCA Coverage", X = "Principal Components",
         y = "% of variances") +
    theme(axis.title=element_text(size= 14, face="bold"),
            axis.text = element_text(size = 14, face="bold"),
            plot.title = element_text(hjust = 0.5,size=14,face="bold"))
  
  ## PCA plot
  ggplot_data <- data.frame(sample =rownames(df_pca$x), 
                            x = df_pca$x[,1],
                            y= df_pca$x[,2])
  
  group = sample.data$condition
  p <-ggplot(ggplot_data,aes(x=x,y=y,label = sample,color=group))+
     geom_point(size=4)+
    geom_vline(xintercept = 0,linetype="dotted", color = "black", size=1) + 
    geom_hline(yintercept=0,linetype="dotted", color = "black", size=1)+
    xlab(paste("PC1:",percentage[1], "%", sep = " ")) +
    ylab(paste("PC2:",percentage[2], "%", sep = " "))+
    theme(axis.title=element_text(size= 14, face="bold"),
          axis.text = element_text(size = 14, face="bold"),
          plot.title = element_text(hjust = 0.5,size=14,face="bold"),
          legend.title = element_text(color = "Black", size = 10, face = "bold"),
          legend.text=element_text(color = "blue", size = 10, face = "bold"),
          legend.position = "bottom",
          legend.margin = margin(1, 1, 1, 1),
          panel.background = element_rect(fill = "white"),
          plot.margin = margin(2, 2, 2, 2, "cm"),
          plot.background = element_rect(fill = "grey90",colour = "black",size = 1))+
    ggtitle("PCA plot")
  
  ## plot both 
  cowplot::plot_grid(screeplot, p,ncol = 2,
                     labels=LETTERS[1:2])  
  
}

## Sample deprogram
plot_sample_den <-function(exp.data,colors, ...){
  dend <- as.dendrogram(hclust(dist(t(exp.data))))
  colors <- ifelse(sample.data$condition == "treated", "skyblue",
                   "orange")
  par(mar=c(9,4,3,3))
  dend %>%
    #set("labels_col", value = c("skyblue", "orange"), k= 2) %>%
    #set("branches_k_color", value = c("skyblue", "orange"), k = 2) %>%
    set("branches_lwd", 2) %>%
    set("labels_cex", 1) %>%
    set("nodes_cex", 0.7) %>%
    plot(axes=TRUE, main = "Sample Dendogram", 
         cex.main = 1.5, font.main= 4)
  # Add the colored bar
  colored_bars(colors = colors, dend = dend, rowLabels = "Label")
  
}


## Volcano plots 
## 1st I did DEA
library(DESeq2)
sample <- sample.data %>% mutate_if(is.character, factor)
rownames(sample) <- rownames(sample.data)
sample.data <- sample
dds <- DESeqDataSetFromMatrix(countData = exp.data,
                              colData = sample.data,
                              design = ~ condition + type)


## keep genes with 50 or higher count in each sample  
keep <- rowSums(counts(dds)) >= 350   
dds <- dds[keep,]
nrow(dds) ## 6929
dds <- DESeq(dds)
resultsNames(dds)


## transformation 
vsd <- vst(dds, blind = FALSE) 
write.csv(as.data.frame(assay(vsd)), file = "vsd_data.csv")

## contrast treated vs untreated 
data <- results(dds, alpha = 0.05, contrast = c("condition","treated", "untreated"))
summary(data) 

## Order by padj
data <- as.data.frame(data[order(exp.data$padj),])
write.csv(data,  file = "DEG.csv")

## volcano plot 
plot_volcano <- function(expression,
                         Significant,
                         title = "Volcano Plot",
                         subtitle = "Treated vs Non-Treated",
                         legend.name,
                         num_sig_gene,
                         color,..){
  
  expression$Significant <- ifelse(expression$padj < 0.05, "Sig", "Not Sig")
  expression$col <- ifelse(expression$log2FoldChange < 0 & expression$padj < 0.05, 'Downregulated',
                                    ifelse(expression$log2FoldChange > 0 & expression$padj < 0.05, "Upregulated",
                                           "none Sig & expressed"))
  ggplot(expression, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = col), size=4) + 
    scale_color_manual(values = c("royalblue", "grey", "red")) +
    geom_vline(xintercept = 0, linetype = 2, colour='black', size = 1) +
    geom_hline(yintercept = 1.3, linetype = 2, colour='black', size = 1)+
    theme_bw(base_size = 12)+ 
    theme(axis.title=element_text(size = 12,face="bold", colour = "black"),
          axis.text = element_text(size = 12,face="bold", colour = "black"),
          axis.line = element_line(colour='black'),
          axis.ticks = element_line(colour='black'),
          plot.title = element_text(hjust = 0.5,size=12,face="bold"),
          plot.subtitle = element_text(hjust = 0.5,size=10,face="bold"),
          legend.position = "bottom",
          legend.title = element_text(color = "Black", size = 10, face = "bold"),
          legend.text=element_text(color = "Black", size = 10, face = "bold")) +
    labs(x = "log2FoldChange",
         y = "-log10(padj)",
         color = "Sig.Level")+
    ggtitle(title) +
    geom_text_repel(data = expression[1:20, ],aes(label = rownames(expression[1:20, ])))

  
  }

## Test Functions
## Read expression and sample data 
exp.data <- read.csv("cts.csv", row.names = 1)
sample.data <- read.csv("Sample.csv", row.names = 1)

## PC plot 
plot_PCA(exp.data)  

## sample distance   
plot_sample_den(exp.data)

## volcano plot for DE genes 
data <- read.csv("DEG.csv", row.names = 1)
plot_volcano(expression = data, subtitle = "Treated vs Non-Treated")
####################################################################################################

