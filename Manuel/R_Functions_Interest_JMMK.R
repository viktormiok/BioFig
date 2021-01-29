#Creator: José Manuel Monroy Kuhn
#Day of creation: 25.01.2021
#Description: This file contains functions of interest for the package BioFig

#Function added 29.01.2020
#'Sample dendrogram
#'@description This function creates a sample dendrogram. It is a wrap up of two functions, plotDendroAndColors from the package WGCNA that creates a sample dendrogram and the function flashClust of the package flashClust that implements optimal hierarchical clustering.
#'@param expression An expression matrix with samples as columns (e.g gene expression)
#'@param colors A vector with colors for every sample according to a specific variable (e.g. "Sex"), or a matrix with colors when more than one variable is used ("Sex", "Age", "Condition", etc.)
#'@param method The agglomeration method to use for clustering ("ward", "single", "complete", "average", "mcquitty", "median" or "centroid"), see help documentation of flashClust::flashClust for more information
#'@param title Title of the plot
#'@param labels A list that contains one or more labels according to the variables used ("Sex","Age","Condition", etc.)
#'@examples
#'
#'@importFrom WGCNA plotDendroAndColors
#'@importFrom flashClust flashClust
#'@importFrom gplots col2hex
sample_dendrogram<-function(expression,colors, method="ward",title="Sample dendrogram",labels=c("Condition")){
  colors<-as.matrix(colors)
  if(ncol(colors)>1){
    traitColors<- apply(colors,MARGIN = 2,FUN = col2hex)
  }else{
    traitColors<-col2hex(colors)
  }
  sample_dendrogram<- flashClust(dist(t(expression)), method = method)
  par(cex =0.6)
  par(mar =c(0,4,2,0))
  plot<-plotDendroAndColors(sample_dendrogram,
                                   traitColors,
                                   groupLabels = labels,
                                   main=title,
                                   marAll = c(1,5,3,1))
  return(plot)
}

#Function added 29.01.2020
#'Volcano plot
#'@description This functions is a modification of the enhanced volcano function from the package enhanced volcano. 
#'I added an option to include the feature names (e.g. genes) with ggrepel according to the log2fold change (logFC) and adjusted p value (adj.P.val). I added an option to change the size of the dots according to a given column (e.g. average expression).Rest of parameters are like EnhancedVolcano.
#'@param keepLab1 logical. If TRUE labels for the features (e.g. genes) with a padj<0.05 and an abs(log2FC)>=2 are shown
#'@param keepLab2 logical. If TRUE labels for the features (e.g. genes) with a padj<0.05 and an 2>abs(log2FC)>1 are shown
#'@param features Name for the features from the input table ("genes", "proteins" or other). Default is genes.
#'@param SizeDots logical. If TRUE dots are resized based on the values of a selected input column
#'@param col_SizeDots Number or name of the column in the input table to use for the size of the plot dots
#'@param legend_SizeDots Optional legend name for dot colors legend
#'@param legend_significance Optional legend name for size legend 
EnhancedVolcano2<-function (toptable, lab, x="logFC", y="adj.P.val", selectLab = NULL, xlim = c(min(toptable[,x], na.rm = TRUE), max(toptable[, x], na.rm = TRUE)),
                            ylim = c(0,max(-log10(toptable[, y]), na.rm = TRUE) + 5), 
                            xlab = bquote(~Log[2] ~"fold change"), 
                            ylab = bquote(~-Log[10] ~ italic(P)), axisLabSize = 18, 
                            title = "Volcano plot", subtitle = "Bioconductor package EnhancedVolcano", 
                            features='genes', #Name of the features in the expression table
                            keepLab1=TRUE, #Add ggrepel labels if padj<0.05 and abs(log2FC)>=2
                            keepLab2=TRUE, #Add ggrepel labels if padj<0.05 and 2>abs(log2FC)>1
                            caption = paste0("Total = ", nrow(toptable), " ",features), 
                            titleLabSize = 18, subtitleLabSize = 14, captionLabSize = 10, 
                            pCutoff = 1e-05, pLabellingCutoff = pCutoff, FCcutoff = 1, 
                            cutoffLineType = "longdash", cutoffLineCol = "black", cutoffLineWidth = 0.4, 
                            transcriptPointSize = 1, transcriptLabSize = 3, transcriptLabCol = "black", 
                            transcriptLabFace = "plain", transcriptLabhjust = 0, transcriptLabvjust = 1.5, 
                            boxedlabels = FALSE, shape = 19, shapeCustom = NULL, col = c("grey30","forestgreen", "royalblue", "red2"), colCustom = NULL, 
                            colAlpha = 1/2, legend = c("NS", "Log2 FC", "P", "P & Log2 FC"), 
                            legendPosition = "top", legendLabSize = 16, legendIconSize = 6, 
                            legendVisible = TRUE, shade = NULL, shadeLabel = NULL, shadeAlpha = 1/2, 
                            shadeFill = "grey", shadeSize = 0.01, shadeBins = 2, drawConnectors = FALSE, 
                            widthConnectors = 0.5, typeConnectors = "closed", endsConnectors = "first", 
                            lengthConnectors = unit(0.01, "npc"), colConnectors = "grey10", 
                            hline = NULL, hlineType = "longdash", hlineCol = "black", 
                            hlineWidth = 0.4, vline = NULL, vlineType = "longdash", vlineCol = "black", 
                            vlineWidth = 0.4, gridlines.major = TRUE, gridlines.minor = TRUE, 
                            border = "partial", borderWidth = 0.8, borderColour = "black",
                            SizeDots = TRUE,
                            col_SizeDots=NULL, #specify column to set the size of the plot dots
                            legend_SizeDots="", #legend name for the size of the dots
                            legend_significance=""){#legend name for the significance
  #Add TRUE/FALSE columns for ggrepel labels
  if(keepLab1){
    toptable$padjLFC2<-ifelse(toptable[,y]<0.05 & (toptable[,x]>=2 | toptable[,x]<=-2),TRUE,FALSE) #It is TRUE if p-adj is less than 0.05 and the abs(log2FC)>=2
    Lab1<-toptable$padjLFC2
  }
  if(keepLab2){
    toptable$significant_padj<-ifelse(toptable[,y]<0.05 & ((toptable[,x]<2 & toptable[,x]>1) | (toptable[,x]>-2 & toptable[,x]< -1)),TRUE,FALSE) #It is TRUE if p-adj is less than 0.05 and 2 > abs(log2FC) > 1
    Lab2<-toptable$significant_padj
  }
  
  if (!requireNamespace("ggplot2")) {
    stop("Please install ggplot2 first.", call. = FALSE)
  }
  if (!requireNamespace("ggrepel")) {
    stop("Please install ggrepel first.", call. = FALSE)
  }
  if (!is.numeric(toptable[, x])) {
    stop(paste(x, " is not numeric!", sep = ""))
  }
  if (!is.numeric(toptable[, y])) {
    stop(paste(y, " is not numeric!", sep = ""))
  }
  if(SizeDots){
    if(is.null(col_SizeDots)){
      stop("Please specify the input column to resize plot dots: col_SizeDots = column_number or 'column_name'")
    }
  }
  i <- xvals <- yvals <- Sig <- NULL
  toptable <- as.data.frame(toptable)
  toptable$Sig <- "NS"
  toptable$Sig[(abs(toptable[, x]) > FCcutoff)] <- "FC"
  toptable$Sig[(toptable[, y] < pCutoff)] <- "P"
  toptable$Sig[(toptable[, y] < pCutoff) & (abs(toptable[, x]) > FCcutoff)] <- "FC_P"
  toptable$Sig <- factor(toptable$Sig, levels = c("NS", "FC", "P", "FC_P"))
  if (min(toptable[, y], na.rm = TRUE) == 0) {
    warning(paste("One or more P values is 0.", "Converting to minimum possible value..."), 
            call. = FALSE)
    toptable[which(toptable[, y] == 0), y] <- .Machine$double.xmin
  }
  toptable$lab <- lab
  toptable$xvals <- toptable[, x]
  toptable$yvals <- toptable[, y]
  if (!is.null(selectLab)) {
    names.new <- rep(NA, length(toptable$lab))
    indices <- which(toptable$lab %in% selectLab)
    names.new[indices] <- toptable$lab[indices]
    toptable$lab <- names.new
  }
  th <- theme_bw(base_size = 24) + theme(legend.background = element_rect(), 
                                         plot.title = element_text(angle = 0, size = titleLabSize, 
                                                                   face = "bold", vjust = 1,hjust = 0.5), 
                                         plot.subtitle = element_text(angle = 0, 
                                                                      size = subtitleLabSize, 
                                                                      face = "plain", vjust = 1), 
                                         plot.caption = element_text(angle = 0, 
                                                                     size = captionLabSize, 
                                                                     face = "plain", vjust = 1), 
                                         axis.text.x = element_text(angle = 0, size = axisLabSize, vjust = 1), 
                                         axis.text.y = element_text(angle = 0, size = axisLabSize, vjust = 1), 
                                         axis.title = element_text(size = axisLabSize), 
                                         legend.position = legendPosition, 
                                         legend.key = element_blank(), 
                                         legend.key.size = unit(0.5, "cm"), 
                                         legend.text = element_text(size = legendLabSize), 
                                         title = element_text(size = legendLabSize),
                                         legend.title = element_text(face="bold")) #+ theme(plot.margin = unit(c(1,1,1,1), "cm"))
  if (!is.null(colCustom) & !is.null(shapeCustom)) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      labs(size=legend_SizeDots,color=legend_significance) +
      th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)), 
                  shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
      geom_point(aes(color = factor(names(colCustom)),size=sapply(SizeDots,function(x){if(x){toptable[,col_averageExp]}else{1}}), shape = factor(names(shapeCustom))), alpha = colAlpha, 
                 na.rm = TRUE) + scale_color_manual(values = colCustom) + 
      #geom_point(aes(color = factor(names(colCustom)), shape = factor(names(shapeCustom))), alpha = colAlpha, 
      #           size = transcriptPointSize, na.rm = TRUE) + scale_color_manual(values = colCustom) + 
      scale_shape_manual(values = shapeCustom)
  }
  else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 1) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      labs(size=legend_SizeDots,color=legend_significance) +
      th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)), 
                  shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
      geom_point(aes(color = factor(names(colCustom)),size=sapply(SizeDots,function(x){if(x){toptable[,col_averageExp]}else{1}})), alpha = colAlpha, shape = shape, size = transcriptPointSize, 
                 na.rm = TRUE) +
      #geom_point(aes(color = factor(names(colCustom))), alpha = colAlpha, shape = shape, size = transcriptPointSize, 
      #           na.rm = TRUE) + 
      scale_color_manual(values = colCustom) + scale_shape_manual(guide = TRUE)
  }
  else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 4) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      labs(size=legend_SizeDots,color=legend_significance) +
      th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)), 
                  shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
      geom_point(aes(color = factor(names(colCustom)),size=sapply(SizeDots,function(x){if(x){toptable[,col_averageExp]}else{1}}), shape = factor(Sig)), alpha = colAlpha, na.rm = TRUE) + 
      #geom_point(aes(color = factor(names(colCustom)), shape = factor(Sig)), alpha = colAlpha, size = transcriptPointSize, na.rm = TRUE) + 
      scale_color_manual(values = colCustom) + 
      scale_shape_manual(values = c(NS = shape[1], FC = shape[2], P = shape[3], FC_P = shape[4]), 
                         labels = c(NS = legend[1],FC = paste(legend[2], sep = ""), P = paste(legend[3], sep = ""), FC_P = paste(legend[4], sep = "")), 
                         guide = TRUE)
  }
  else if (is.null(colCustom) & !is.null(shapeCustom)) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      labs(size=legend_SizeDots,color=legend_significance) +
      th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)), 
                  shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
      geom_point(aes(color = factor(Sig), shape = factor(names(shapeCustom)),size=sapply(SizeDots,function(x){if(x){toptable[,col_averageExp]}else{1}})), alpha = colAlpha, na.rm = TRUE) + 
      #geom_point(aes(color = factor(Sig), shape = factor(names(shapeCustom))), alpha = colAlpha, size = transcriptPointSize, na.rm = TRUE) + 
      scale_color_manual(values = c(NS = col[1],
                                    FC = col[2], P = col[3], 
                                    FC_P = col[4]), 
                         labels = c(NS = legend[1], FC = paste(legend[2], sep = ""), 
                                    P = paste(legend[3], sep = ""), 
                                    FC_P = paste(legend[4], sep = ""))) + 
      scale_shape_manual(values = shapeCustom)
  }
  else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) == 1) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      labs(size=legend_SizeDots,color=legend_significance) +
      th + guides(colour = guide_legend(order = 1, override.aes = list(shape = shape, size = legendIconSize))) + 
      geom_point(aes(color = factor(Sig),size=sapply(SizeDots,function(x){if(x){toptable[,col_averageExp]}else{1}})), alpha = colAlpha, shape = shape, na.rm = TRUE, show.legend = legendVisible) + 
      #geom_point(aes(color = factor(Sig)), alpha = colAlpha, shape = shape, size = transcriptPointSize, na.rm = TRUE, show.legend = legendVisible) + 
      scale_color_manual(values = c(NS = col[1], FC = col[2], P = col[3], FC_P = col[4]), labels = c(NS = legend[1],FC = paste(legend[2], sep = ""), P = paste(legend[3], sep = ""), FC_P = paste(legend[4], sep = "")))
  }
  else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) == 4) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      labs(size=legend_SizeDots,color=legend_significance) +
      th + guides(colour = guide_legend(order = 1, override.aes = list(shape = c(NS = shape[1], FC = shape[2], P = shape[3], FC_P = shape[4]), size = legendIconSize))) + 
      geom_point(aes(color = factor(Sig), shape = factor(Sig),size=sapply(SizeDots,function(x){if(x){toptable[,col_averageExp]}else{1}})), 
                 alpha = colAlpha, size = transcriptPointSize, 
                 na.rm = TRUE, show.legend = legendVisible) +
      geom_point(aes(color = factor(Sig), shape = factor(Sig),size=sapply(SizeDots,function(x){if(x){toptable[,col_averageExp]}else{1}})), 
                 alpha = colAlpha, 
                 na.rm = TRUE, show.legend = legendVisible) + 
      #geom_point(aes(color = factor(Sig), shape = factor(Sig)), 
      #           alpha = colAlpha, size = transcriptPointSize, 
      #           na.rm = TRUE, show.legend = legendVisible) + 
      scale_color_manual(values = c(NS = col[1], FC = col[2], 
                                    P = col[3], FC_P = col[4]), 
                         labels = c(NS = legend[1], 
                                    FC = paste(legend[2], sep = ""), 
                                    P = paste(legend[3], sep = ""), 
                                    FC_P = paste(legend[4], sep = ""))) + 
      scale_shape_manual(values = c(NS = shape[1], FC = shape[2], 
                                    P = shape[3], FC_P = shape[4]), guide = FALSE)
  }
  
  plot <- plot + xlab(xlab) + ylab(ylab) + xlim(xlim[1], xlim[2]) + 
    ylim(ylim[1], ylim[2]) + 
    geom_vline(xintercept = c(-FCcutoff, FCcutoff), linetype = cutoffLineType, 
               colour = cutoffLineCol, size = cutoffLineWidth) + 
    geom_hline(yintercept = -log10(pCutoff), linetype = cutoffLineType, 
               colour = cutoffLineCol, size = cutoffLineWidth)
  
  if(keepLab1){
    plot<-plot+ggrepel::geom_text_repel(data = dplyr::filter(toptable, Lab1), 
                                        aes(label = lab), size = transcriptLabSize, 
                                        box.padding = unit(0.1, "lines"), 
                                        point.padding = unit(0.1, "lines"), 
                                        segment.size = 0.5,
                                        segment.alpha=0.4,
                                        colour='red')
  }
  
  if(keepLab2){
    plot<-plot+ggrepel::geom_text_repel(data = dplyr::filter(toptable, Lab2), 
                                        aes(label = lab), size = 4, 
                                        box.padding = unit(0.1, "lines"), 
                                        point.padding = unit(0.1, "lines"), 
                                        segment.size = 0.5,
                                        segment.alpha=0.5,
                                        colour='grey',
                                        alpha=0.7)
  }
  
  
  plot <- plot + labs(title = title, caption = caption)
  if(!SizeDots){
    plot<-plot + scale_size(guide="none")
  }
  if (!is.null(vline)) {
    plot <- plot + geom_vline(xintercept = vline, linetype = vlineType, 
                              colour = vlineCol, size = vlineWidth)
  }
  if (!is.null(hline)) {
    plot <- plot + geom_hline(yintercept = -log10(hline), 
                              linetype = hlineType, colour = hlineCol, size = hlineWidth)
  }
  if (border == "full") {
    plot <- plot + theme(panel.border = element_rect(colour = borderColour, 
                                                     fill = NA, size = borderWidth))
  }
  else if (border == "partial") {
    plot <- plot + theme(axis.line = element_line(size = borderWidth, 
                                                  colour = borderColour), panel.border = element_blank(), 
                         panel.background = element_blank())
  }
  else {
    stop("Unrecognised value passed to 'border'. Must be 'full' or 'partial'")
  }
  if (gridlines.major == TRUE) {
    plot <- plot + theme(panel.grid.major = element_line())
  }
  else {
    plot <- plot + theme(panel.grid.major = element_blank())
  }
  if (gridlines.minor == TRUE) {
    plot <- plot + theme(panel.grid.minor = element_line())
  }
  else {
    plot <- plot + theme(panel.grid.minor = element_blank())
  }
  if (!is.null(shade)) {
    plot <- plot + stat_density2d(data = subset(toptable, 
                                                rownames(toptable) %in% shade), fill = shadeFill, 
                                  alpha = shadeAlpha, geom = "polygon", contour = TRUE, 
                                  size = shadeSize, bins = shadeBins, show.legend = FALSE, 
                                  na.rm = TRUE) + scale_fill_identity(name = shadeLabel, 
                                                                      labels = shadeLabel, guide = "legend")
  }
  return(plot)
}

#Function added 29.01.2020
#'PCA plot
#'@description This function plots a PCA
#'@param expression An expression matrix with samples as columns (e.g gene expression)
#'@param group Grouping for color (factor variable, e.g. treatment)
#'@param colors Optional color vector for chosen group variable
#'@param shape Optional grouping for shape (factor variable, e.g. sex)
#'@param samplenames Names of samples
#'@param LegendName_Color Legend name for color, default is 'group'
#'@param LegendName_Shape Legend name for shape, default is 'shape'
#'@param ggrepelLab logical. If TRUE labels for every sample in the PCA plot are displayed
#'@param size_gglab Size for the ggrepel labels that show sample names
plot_2DPCA<-function(expression,group,colors=NULL,shape=NULL,samplenames,title="",LegendName_Color="group",LegendName_Shape="shape",LegendName="group",ggrepelLab=TRUE,size_gglab=5){
  df_pca<-prcomp(t(expression))
  df_out <- as.data.frame(df_pca$x)
  df_out$group<-group
  df_out$sample_name<-samplenames
  
  percentage <- round(df_pca$sdev^2 / sum(df_pca$sdev^2)*100,2)
  percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )
  print(percentage)
  
  if(is.null(shape)){
    p<-ggplot(df_out,aes(x=PC1,y=PC2,color=group))
    if(!is.null(colors)){p<-p+scale_color_manual(values=colors,name=LegendName)}
  }else{
    p<-ggplot(df_out,aes(x=PC1,y=PC2,color=group,shape=shape))
    if(!is.null(colors)){p<-p+scale_color_manual(values=colors,name=LegendName)}
  }
  
  
  p<-p+geom_point(size=4)+ xlab(percentage[1]) + ylab(percentage[2])+
    #geom_text(label=row.names(df_out_raw),show.legend = FALSE,hjust=0,vjust=0.2)
    theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5),
          legend.title=element_text(size=14), legend.text=element_text(size=12,),
          axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),
          panel.background = element_blank(),
          panel.grid.major =  element_line(colour = "grey90", size = 0.2),
          panel.grid.minor =  element_line(colour = "grey98", size = 0.5),
          panel.border = element_rect(color='black',fill=NA))
  if(ggrepelLab){
    p<-p + geom_text_repel(aes(label=sample_name),show.legend = FALSE,size=size_gglab,force = 2) +
      labs(shape=LegendName_Shape, col=LegendName_Color)+
      ggtitle(title)
  }
  
  return(p)
}

#Function added 29.01.2020
#'3D PCA plot
#' @description This function creates a 3D PCA plot using the first three principal components
#' @param expression Expression matrix
#' @param color Vector of colors as text (e.g.'blue'), one color for each sample
#' @importFrom rgl plot3d grid3d text3d
plot_3DPCA <- function(expression, color) {
  dat<-as.matrix(t(expression))
  dat<-as.data.frame(dat)
  
  #Compute PCA
  pca.data<-prcomp(dat)
  scores = as.data.frame(pca.data$x)
  x<-pca.data$x[,c(1:3)]
  
  n <- ncol(x) 
  if(!(n %in% c(3))) { # only 3D PCA
    stop("x must have 3 columns")
  }
  
  if(n == 3) { # 3d plot
    plot3d(x, col=color, type="s", size=1, axes=F)
    axes3d(edges=c("x--", "y--", "z"), lwd=3, axes.len=2, labels=FALSE)
    grid3d("x")
    grid3d("y")
    grid3d("z")
  } 
  
  # text3d(scores[,1]+0.7, scores[,2]+0.7, scores[,3]+0.7,
  #        texts=c(rownames(scores)), cex= 0.8, pos=3)
  text3d(scores[,1], scores[,2], scores[,3],
         texts=c(rownames(scores)), cex= 0.8, pos=3)
}



#Other possible functions to add:
#GO enrichment bubble plot
#Alluvial plot
#Other 3D plots using plotly, modify PCA plot and do it with plotly?
