#' Plot exploratory plots
#' @param expression The expression table to use
#' @param group Factor specifying the groups for the samples or the name of the column containing the factor variable if the input is a deseq2 object
#' @param samplenames Vector specifying the names of the samples or the name of the column containing the sample names if the input is a deseq2 object
#' @param plottype Type of plot to generate, options are 'PCA','MDS','Distance','Expression' and 'all'
#' @param transform Apply a transformation to the data, available transformation are 'vst' and 'rlog'
#' @param colors Optional vector of colors for the factor levels
#' @param shape Optional vector of shapes for the factor levels
#' @param LegendName_Color Title for the color legend
#' @param LegendName_Shape Title for the shape legend (for PCA)
#' @param ggrepelLab logical. If TRUE sample names are included as ggrepel labels
#' @param size_gglab Size of the ggrepel labels
#' @param size_title Size of plots titles
#' @param point.size Size of the dots for the PCA and MDS plots
#' @param scale logical. If explicitly set to TRUE, the data is scaled for the PCA plot even if a transformation was applied 
#' @param ntop Use n features that explain most of the variance
#' @param MahalanobisEllips logical. If TRUE mahalanobis ellipses are included in the PCA plot
#' @return This function returns one or four plots for data exploration. The available plots are Principal Component Analysis (PCA), Multi Dimensional Scaling (MDS), Heatmap of sample distances, and a Heatmap of the expression of all features.
#' @examples
#' @import ComplexHeatmap
#' @import dplyr
#' @import ggplot2
#' @import ggrepel
#' @export
plot_ExploratoryPlots<-function(expression,group,samplenames,
                                plottype = c("all","PCA", "MDS", "Distance","Expression"),
                                transform="no",
                                colors=NULL,shape=NULL,
                                LegendName_Color="group",
                                LegendName_Shape="shape",
                                ggrepelLab=TRUE,size_gglab=5,size_title=1,
                                point.size=4,scale=T,
                                ntop=NULL,
                                MahalanobisEllips=F,...){
  
  #Internal Functions
  matchLs<-function(L1,L2){
    Idx<-match(L1,L2)
    IdxOut<-Idx[!is.na(Idx)]
    IdxOut
  }
  '%ni%' <- Negate('%in%')
  scl<-T
  
  #Obtain parameters given by the user
  ParaList<-as.list(match.call())
  
  #Check parameters that should be characters
  ParamChNames<-  c("group","colors","samplenames","title","LegendName_Color","LegendName_Shape","LegendName","transform", "plottype")
  ChParaList<-ParaList[matchLs(ParamChNames,names(ParaList))]   #Obtain the parameter values given by the user
  ChNameList <-ParamChNames[matchLs(names(ParaList),ParamChNames)] #and names 
  for(i in 1:length(ChParaList)){#Loop parameters to make sure they are of class character
    param<-eval(ChParaList[[i]])
    if(is.null(param) && ChNameList[i] == "colors"){
      next
    }else if(is.null(param)){
      stop(paste0("Input ",ChNameList[i]," cannot be null"))
    }
    if(ChNameList[i]=="group" && is.factor(param)){
      param<-as.character(param)
    }
    if (!is(param, "character")) {
      stop(paste0("Input ",ChNameList[i]," is of a wrong class. Character is expected and ",class(param)," was given"))
    }
    if(ChNameList[i]=="transform" && param %ni% c("no","vst","rlog")){
      stop(paste0("Input ",ChNameList[i]," is not part of the transformations. Only '",paste(c("no","vst","rlog"),collapse=", "),"' are possible parameters"))
    }
  }
  
  #Check parameters that should be nummeric
  ParamNumNames<-  c("shape","size_gglab","size_title","point.size","ntop")
  NumParaList<-ParaList[matchLs(ParamNumNames,names(ParaList))]
  NumNameList <-ParamNumNames[matchLs(names(ParaList),ParamNumNames)]
  for(i in 1:length(NumParaList)){#Loop parameters to make sure they are of class numeric
    param<-eval(NumParaList[[i]])
    if(is.null(param) && NumNameList[i] %in% c("shape","ntop")){
      next
    }else if(is.null(param)){
      stop(paste0("Input ",NumNameList[i]," cannot be null"))
    }
    
    if (!is(param, "numeric")) {
      stop(paste0("Input ",NumNameList[i]," is of a wrong class. Numeric is expected and ",class(param)," was given"))
    }
    if(param<0){
      stop(paste0("Input ",NumNameList[i]," is not a positive integer"))
    }
  }
  
  #Check parameters that should be logical
  ParamLogNames<-  c("scale","MahalanobisEllips")
  LogParaList<-ParaList[matchLs(ParamLogNames,names(ParaList))]
  LogNameList <-ParamLogNames[matchLs(names(ParaList),ParamLogNames)]
  for(i in 1:length(LogParaList)){#Loop parameters to make sure they are of class numeric
    param<-eval(LogParaList[[i]])
    if (!is(param, "logical")) {
      stop(paste0("Input ",LogNameList[i]," is of a wrong class. Logical is expected and ",class(param)," was given"))
    }
  }
  
  
  if(class(expression)=="DESeqDataSet"){ #Deseq2 option
    group<-colData(expression)[,group]
    samplenames<-colData(expression)[,samplenames]
    expression<-assay(expression)
  }else if(class(expression)=="DGEList" && attr(attributes(expression)[[1]],"package")=="edgeR"){  #EdgeR option
    #group <- expression[[2]][1]
    expression<-edgeR::getCounts(expression)
    samplenames<-colnames(expression)
  }
  
  #Add option transform data (could be packed in another function)
  transform<-match.arg(transform) #MAybe change it later to a simple if(is.na())
  if(!transform=="no"){
    eval_r<-switch(transform, "vst" = vst(expression),"rlog" = rlog(expression))
    expression<-eval_r
    scl<-F
    if(scale){scl<-T}
  }else{
    if(scale){scl<-T}
  }
  
  #Add option to use only the top n genes that explain most of the variance
  if(!is.null(ntop)){
    rv <- rowVars(expression)
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    expression<-expression[select,]
  }
  
  
  
  switch(plottype,
         all={
           #PCA
           p_PCA<-plot_2DPCA(expression=expression,scl=scl,
                             group=group,shape=shape,samplenames=samplenames,
                             colors=colors,ggrepelLab = ggrepelLab,
                             size_gglab = size_gglab,size_title = size_title,
                             point.size = point.size,LegendName_Color = LegendName_Color,
                             LegendName_Shape = LegendName_Shape,
                             MahalanobisEllips = MahalanobisEllips,...)
           #MDS
           p_MDS<-plot_MDS(expression=expression,group=group,point.size=point.size,
                           LegendName_Color = LegendName_Color,
                           LegendName_Shape = LegendName_Shape)
           
           #Distance
           pList<-plot_SampleDistance(expression=expression,...)
           p_Similarity<-pList[[2]]$gtable
           
           #Expression
           p_heatmap<-plot_ExpressionHeatmap(expTable=expression,sampleDGroup=group)
           grob_expression = grid.grabExpr(draw(p_heatmap)) 
           
           cowplot::plot_grid(p_PCA,p_MDS,p_Similarity,grob_expression)
           
         },
         PCA={
           p<-plot_2DPCA(expression=expression,scl=scl,
                         group=group,shape=shape,samplenames=samplenames,
                         colors=colors,ggrepelLab = ggrepelLab,
                         size_gglab = size_gglab,size_title = size_title,
                         point.size = point.size,
                         LegendName_Color = LegendName_Color,
                         LegendName_Shape = LegendName_Shape,
                         MahalanobisEllips = MahalanobisEllips,...)
           plot(p)
           ggsave(file="PCA.svg", plot=p, width=10, height=8)
           ggsave(file="PCA.pdf", plot=p, width=10, height=8)
         },
         MDS={
           p<-plot_MDS(expression=expression,group=group,point.size=point.size,
                       LegendName_Color = LegendName_Color,
                       LegendName_Shape = LegendName_Shape)
           ggsave(file="MDS.svg", plot=p, width=10, height=8)
           ggsave(file="MDS.pdf", plot=p, width=10, height=8)
           plot(p)
         },
         Distance={
           pList<-plot_SampleDistance(expression=expression,...)
           plotsSimilarityDisimilarity<-cowplot::plot_grid(pList[[1]]$gtable,pList[[2]]$gtable, ncol=2)
           ggsave(file="SampleDistances.svg", plot=plotsSimilarityDisimilarity, width=15, height=15)
           ggsave(file="SampleDistances.pdf", plot=plotsSimilarityDisimilarity, width=15, height=15)
         },
         Expression={
           p<-plot_ExpressionHeatmap(expTable=expression,sampleDGroup=group)
           draw(p)
           # return(p)
         },
         stop(paste0("plottype: ",plottype," is not supported. The available options are 'all','PCA','MDS','Distance','Expression'"))
         
  )
  
  
  
}



#' Plot 2D PCA
#' @param expression The expression table to use to create the PCA
#' @param group Factor specifying the groups for the samples
#' @param samplenames Vector specifying the names of the samples
#' @param scl logical. If TRUE the data is scaled
#' @param colors Optional vector of colors for the factor levels
#' @param shape Optional vector of shapes for the factor levels
#' @param ggrepelLab logical. If TRUE sample names are included as ggrepel labels
#' @param size_gglab Size of the ggrepel labels
#' @param size_title Size of plot title
#' @param point.size Size of the dots
#' @param MahalanobisEllips logical. If TRUE mahalanobis ellipses are included in the plot
#' @param LegendName_Color Title for the color legend
#' @param LegendName_Shape Title for the shape legend
#' @return This function return a PCA plot
#' @examples
#' @import ggplot2
#' @import ggrepel
#' @export
plot_2DPCA<-function(expression, group,samplenames,
                     scl,colors=NULL,shape=NULL,
                     ggrepelLab=TRUE,size_gglab=5,
                     size_title=1,point.size=4,
                     MahalanobisEllips=F,
                     LegendName_Color="group",LegendName_Shape="shape",...){
  df_pca<-prcomp(t(expression),scale=scl)
  df_out <- as.data.frame(df_pca$x)
  df_out$group<-group
  df_out$sample_name<-samplenames
  
  percentage <- round(df_pca$sdev^2 / sum(df_pca$sdev^2)*100,2)
  percentage <- paste0(colnames(df_out)[grep("^PC",colnames(df_out))], " (", paste0(as.character(percentage), "% variance", ")") )
  print(percentage)
  
  if(is.null(shape)){
    p<-ggplot(df_out,aes(x=PC1,y=PC2,color=group))
    if(!is.null(colors)){p<-p+scale_color_manual(values=colors,name=LegendName_Color)}
  }else{
    p<-ggplot(df_out,aes(x=PC1,y=PC2,color=group,shape=shape))
    if(!is.null(colors)){p<-p+scale_color_manual(values=colors,name=LegendName_Color)}
  }
  
  
  p<-p+geom_point(size=point.size) + 
    xlab(percentage[1]) + 
    ylab(percentage[2]) +
    ggtitle(label = "PCA plot") +
    labs(shape=LegendName_Shape, col=LegendName_Color)+
    ggTheme(1)
  if(ggrepelLab){
    p<-p + geom_text_repel(aes(label=sample_name),show.legend = FALSE,
                           size=size_gglab,
                           force = 2,
                           max.overlaps = Inf)
  }
  
  if(MahalanobisEllips){
    p<-p+stat_ellipse()
  }
  return(p)
  
  
}


#' Expression Heatmap
#' @param expTable The expression table to use to create the heatmap
#' @param Title Title of the plot
#' @param sampleDGroup A data frame specifying groups/factors for every sample
#' @param clColumns logical. If TRUE columns will be clustered
#' @param FontSRow Font size for the row labels
#' @param FontSColumn Font size for the column labels
#' @param FontSTitle Font size for the title
#' @param setseed Set seed to obtain the same colors
#' @return It returns a heatmap of the expression of all features
#' @import ComplexHeatmap
#' @export
plot_ExpressionHeatmap<-function(expTable,Title="Heatmap",sampleDGroup,clColumns=FALSE,FontSRow=10,FontSColumn=10,FontSTitle=16,setseed=NULL,...){
  #Scale data using columns...
  z <- t(scale(t(expTable)))
  z<-z[complete.cases(z), ] #Keep only those that do not have NA
  zlim <- c(-3,3)
  z <- pmin(pmax(z, zlim[1]), zlim[2])
  
  
  
  getColorsNamedVector<-function(FactorLevels,nLev){
    #Get colors
    n<-sample(0:20,1)
    qual_col_pals = RColorBrewer::brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    start<-1+n
    end<-nLev+n
    colorLevels<-col_vector[start:end]
    return(colorLevels)
  }
  
  sampleDGroup<-as.data.frame(sampleDGroup)
  
  n<-0
  colorDF<-apply(sampleDGroup,MARGIN = 2,function(x){
    n<<-n+20
    FactorLevels<-unique(x) #get unique factor levels
    nLev<-length(FactorLevels) #Number of levels
    if(!is.null(setseed)){
      set.seed(seed=setseed+n)
    }
    colorLevels<-getColorsNamedVector(FactorLevels=FactorLevels,nLev=nLev)
    names(colorLevels)<-FactorLevels
    colorLevels
  })
  
  if(!is.list(colorDF)){
    colorDF<-as.list(colorDF[,1])
  }
  
  
  featureIds<-rownames(z)#as.matrix(EnsemblIDToGeneID[match(rownames(z),EnsemblIDToGeneID$ENSEMBL_ID),2])[,1]
  
  #Prepare Annotation
  ha_column = HeatmapAnnotation(df = sampleDGroup, 
                                col= colorDF)
  
  ht1 = Heatmap(z, name="z-score", top_annotation = ha_column, 
                column_title = Title,row_labels = featureIds,
                column_names_rot = 45,column_names_gp = gpar(fontsize = FontSColumn),row_names_gp = gpar(fontsize = FontSRow),
                cluster_columns = clColumns,show_row_names=TRUE,
                column_title_gp = gpar(fontsize = FontSTitle, fontface = "bold"),...)
  
  return(ht1)
}


# Function MDS
#'@description This function create a multidimensional scaling (MDS) plot that visualizing the level of (dis)similarity between individual cases of a dataset.
#'MDS is a distance-preserving algorithms defined by the pairwise distances of data points.
#'MDS takes a matrix D where Dij represents the dissimilarity between points i and j and produces a mapping on a lower dimension, preserving the dissimilarities as closely as possible. 
#'There are many different ways of calculating dissimilarity among samples, Euclidean distance the default here 
#'@param expression Numerical matrix with samples as columns (e.g gene expression)
#'@param group logical. If TRUE labels for the sample are shown 
#'@param point.size Name of the input table
#'@param LegendName_Color logical. If TRUE dots are resized based on the values of a selected input column
#'@param LegendName_Shape Number or name of the column in the input table to use for the size of the plot dots
#'@return It returns a heatmap of the expression of all features
#'@examples
#'@import tidyverse 
#'@import ggrepel
#'
plot_MDS<-function(expression,group,point.size,
                   LegendName_Color="group",
                   LegendName_Shape="shape",...){
  
  ## calculate distance for the sample
  data <- expression %>%
    t() %>%
    dist() %>%
    as.matrix()
  
  ## convert distance matrix to Classical multidimensional scaling(MDS)
  mdsData <- data.frame(cmdscale(data))
  mds <- cbind(mdsData, as.data.frame(data)) # combine with distance with mds
  
  ## plot in ggplot2
  plotmds <- ggplot(mds, aes(X1, X2,color=group)) +
    geom_point(size = point.size) +
    ggTheme(1) +
    labs(x = "Leading LogFC dim 1", y = "Leading LogFC dim 2", title = "MDS plot") +
    labs(shape=LegendName_Shape, col=LegendName_Color)+
    ggrepel::geom_text_repel(data = mds,aes(label = rownames(mds)),max.overlaps = Inf)
  return(plotmds)
  
  
}



# Sample distance heatmap
#'@description This function created to show sample (dis) similarity in a heatmap using distance or correlation methods. 
#'@param expression Numerical matrix with samples as columns (e.g gene expression)
#'@return It returns a heatmap of the expression of all features
#'@import tidyverse 
#'@import pheatmap 
#'@import RColorBrewer
#'@example 
#'
plot_SampleDistance<-function(expression,...){
  ## calculate distance(Dissimilarity Measure )
  sampleDis <- pheatmap::pheatmap(expression %>%
                                    t() %>%
                                    dist() %>%
                                    as.matrix(),
                                  col = rev(RColorBrewer::brewer.pal(n = 8, name = "RdBu")),
                                  main = "Sample dissmilarity",
                                  fontsize = 9.23,
                                  fontsize_col = 12,
                                  fontsize_row = 12)
  
  
  ## calculate sample correlation(Similarity Measure )
  sampleCor <- pheatmap::pheatmap(expression %>%
                                    cor() %>%
                                    as.matrix(),
                                  col = rev(RColorBrewer::brewer.pal(n = 8, name = "RdBu")),
                                  main = "Sample similarity",
                                  fontsize = 9.23,
                                  fontsize_col = 12,
                                  fontsize_row = 12)
  listplots<-list()
  listplots[[1]]<-sampleDis
  listplots[[2]]<-sampleCor
  return(listplots)
}


#' @keywords Internal
ggTheme<-function(theme=1,size_title=12){
  if(theme==1){
    ggTheme<-theme_minimal() +
      theme(plot.title = element_text(size = size_title, face = "bold",hjust = 0.5),
            legend.title = element_text(color = "Black", size = 12, face = "bold"),
            legend.text=element_text(color = "Black", size = 12, face = "bold"),
            legend.position = "bottom",
            axis.title=element_text(size = 12,face="bold", colour = "black"),
            axis.text = element_text(size = 12),
            axis.ticks = element_line(colour='black'),
            panel.background = element_blank(),
            panel.grid.major =  element_line(colour = "grey90", size = 0.2),
            panel.grid.minor =  element_line(colour = "grey98", size = 0.5),
            panel.border = element_rect(color='black',fill=NA))
    return(ggTheme)
  }
}



