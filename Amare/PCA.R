# The BioFig package contains a collection of R functions used by biologists with some programming knowledge. 
# The package is applicable to visualize the high dimension data to understand biological questions.
plot_2DPCA<-function(expression = c("matrix", "data.frame", "DGEList", "DESeqDataSet"),
                     group,colors=NULL,shape=NULL,
                     samplenames,title="",LegendName_Color="group",
                     LegendName_Shape="shape",LegendName="group",
                     ggrepelLab=TRUE,size_gglab=5,size_title=14,
                     point.size=4,scl=T,ntop=NULL,transform=NULL,MahalanobisEllips=F){
  
  #Add option to use only the top n genes that explain most of the variance
  #Add option transform data (could be packed in another function)
  if(!is.null(transform) && transform=="vst"){
    expression<-vst(expression,fitType = "local")
    scl<-F
  }
  #If the class of the object is dds
  ## edger approach
  if(class(expression)=="DESeqDataSet"){
    group<-colData(expression)[,group]
    samplenames<-colData(expression)[,samplenames]
    expression<-assay(expression)
  }
  #If the class of the object is DGEList
  if(class(expression)=="DGEList"){
    edgR_list_object<- edgeR::plotMDS.DGEList(expression)
    ## Then extract MDS matrix from edgR_list_object 
    ## Save the extracted MDS matrix as ed_list
    ed_list = as.data.frame(edgR_list_object$cmdscale.out) %>% dplyr::rename(DM1 = V1, DM2 = V2)
    
    ## Then pass the the data to ggplot2 to get the plot 
    pplot <- ggplot(ed_list, aes(x=DM1,y=DM2))+
    geom_point(size=3) + 
    theme(axis.title=element_text(size = 12,face="bold", colour = "black"),
                     axis.text = element_text(size = 12),
                     axis.ticks = element_line(colour='black'),
                     plot.title = element_text(hjust = 0.5,size=12,face="bold"),
                     legend.position = "bottom",
                     legend.title = element_text(color = "Black", size = 12, face = "bold"),
                     legend.text=element_text(color = "Black", size = 12, face = "bold"),
                     panel.background = element_blank(),
                     panel.grid.major =  element_line(colour = "grey90", size = 0.2),
                     panel.grid.minor =  element_line(colour = "grey98", size = 0.5),
                     panel.border = element_rect(color='black',fill=NA)) +
     labs(x = "Leading LogFC dim 1", y = "Leading LogFC dim 2", title = "MDS plot") +
      ggrepel::geom_text_repel(data = ed_list,aes(label = rownames(ed_list)))
    return(pplot)
     } 
     
  ## Then, it follows the code below for normalization and visualization.
  df_pca<-prcomp(t(expression),scale=scl)
  df_out <- as.data.frame(df_pca$x)
  df_out$group<-group
  df_out$sample_name<-samplenames
  
  percentage <- round(df_pca$sdev^2 / sum(df_pca$sdev^2)*100,2)
  percentage <- paste0(colnames(df_out)[grep("^PC",colnames(df_out))], " (", paste0(as.character(percentage), "% variance", ")") )
  print(percentage)
  
  if(is.null(shape)){
    p<-ggplot(df_out,aes(x=PC1,y=PC2,color=group))
    if(!is.null(colors)){p<-p+scale_color_manual(values=colors,name=LegendName)}
  }else{
    p<-ggplot(df_out,aes(x=PC1,y=PC2,color=group,shape=shape))
    if(!is.null(colors)){p<-p+scale_color_manual(values=colors,name=LegendName)}
  }
  
  p<-p+geom_point(size=point.size)+ xlab(percentage[1]) + ylab(percentage[2])+
    #geom_text(label=row.names(df_out_raw),show.legend = FALSE,hjust=0,vjust=0.2)
    theme(plot.title = element_text(size = size_title, face = "bold",hjust = 0.5),
          legend.title=element_text(size=14), legend.text=element_text(size=12,),
          axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),
          panel.background = element_blank(),
          panel.grid.major =  element_line(colour = "grey90", size = 0.2),
          panel.grid.minor =  element_line(colour = "grey98", size = 0.5),
          panel.border = element_rect(color='black',fill=NA))
  if(ggrepelLab){
    p<-p + geom_text_repel(aes(label=sample_name),show.legend = FALSE,
                           size=size_gglab,
                           force = 2,max.overlaps = Inf) +
      labs(shape=LegendName_Shape, col=LegendName_Color)+
      ggtitle(title)
  }
  
  if(MahalanobisEllips){
    p<-p+stat_ellipse()
  }
  return(p)
}


## Test Functions
## PC plot 
library(dplyr)
library(ggplot2)
expression <- read.csv("cts.csv", row.names = 1)
expression <- expression[apply(expression, 1, function(row) all(row !=0)), ]
expression <- log2(expression)

Sampledata <- read.csv("Sample.csv", row.names = 1) %>% dplyr::select(-2)
expression <- edgeR::DGEList(counts = expression, group = t(Sampledata))

debug(plot_2DPCA)
plot_2DPCA(expression)


# calculate sample distances
library(tidyverse)
plotmds <- function(data){
  
  ## calculate distance for the sample
  data <- expression %>%
    t() %>%
    dist() %>%
    as.matrix()
  
  ## convert distance matrix to Classical multidimensional scaling(MDS)
  mdsData <- data.frame(cmdscale(data))
  mds <- cbind(mdsData, as.data.frame(data)) # combine with distance with mds
  
  ## plot in ggplot2
  ggplot(mds, aes(X1, X2)) +
    geom_point(size = 3) +
    theme_minimal() +
    theme(axis.title=element_text(size = 12,face="bold", colour = "black"),
          axis.text = element_text(size = 12),
          axis.ticks = element_line(colour='black'),
          plot.title = element_text(hjust = 0.5,size=12,face="bold"),
          legend.position = "bottom",
          legend.title = element_text(color = "Black", size = 12, face = "bold"),
          legend.text=element_text(color = "Black", size = 12, face = "bold"),
          panel.background = element_blank(),
          panel.grid.major =  element_line(colour = "grey90", size = 0.2),
          panel.grid.minor =  element_line(colour = "grey98", size = 0.5),
          panel.border = element_rect(color='black',fill=NA)) +
    labs(x = "Leading LogFC dim 1", y = "Leading LogFC dim 2", title = "MDS plot") +
    ggrepel::geom_text_repel(data = mds,aes(label = rownames(mds)))
  
}
  





