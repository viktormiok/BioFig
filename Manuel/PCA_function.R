expression <- read.csv("Data_expression_normalized.csv", row.names = 1)
head(expression) 

plot_2DPCA<-function(expression,group,colors=NULL,shape=NULL,
                     samplenames,title="PCA",LegendName_Color="group",
                     LegendName_Shape="shape",LegendName="group",
                     ggrepelLab=TRUE,size_gglab=5,size_title=1,
                     point.size=4,scl=T,ntop=NULL,transform=NULL,
                     MahalanobisEllips=F, plottype = c("PCA", "MDS", "sample/geneDistance")){
  
  #Functions used
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
  
  # if (!is(scale, "logical")) {
  #   stop("Input (scale) is of wrong class.")
  # } 
  # if (!is(MahalanobisEllips, "logical")) {
  #   stop("Input (MahalanobisEllips) is of wrong class.")
  # } 
  
  # if (!is(plottype, "character")) {
  #   stop("Input (plottype) is of wrong class.")
  # }
  
  #Deseq2 option
  if(class(expression)=="DESeqDataSet"){
    group<-colData(expression)[,group]
    samplenames<-colData(expression)[,samplenames]
    expression<-assay(expression)
  }
  
  #EdgeR option
  if(class(expression)=="DGEList" && attr(attributes(expression)[[1]],"package")=="edgeR"){
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
  
  if(plotype == "PCA"){
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
  
  ## MDS
  if(plotype == "MDS"){
    
    ## calculate distance for the sample
    data <- expression %>%
      t() %>%
      dist() %>%
      as.matrix()
    
    ## convert distance matrix to Classical multidimensional scaling(MDS)
    mdsData <- data.frame(cmdscale(data))
    mds <- cbind(mdsData, as.data.frame(data)) # combine with distance with mds
    
    ## plot in ggplot2
    plotmds <- ggplot(mds, aes(X1, X2)) +
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
    return(plotmds)
    
  }
  
  ## Sample Heatmap
  if(plotype == "sample/geneDistance"){
    ## calculate distance(Dissimilarity Measure )
    sampleDis <- pheatmap::pheatmap(expression %>%
                                      t() %>%
                                      dist() %>%
                                    as.matrix(),
                                    col = rev(RColorBrewer::brewer.pal(n = 8, name = "RdBu")),
                                    main = "Sample dissmilarity",
                                    fontsize = 14)
    
    
    ## calculate sample correlation(Similarity Measure )
    sampleCor <- pheatmap::pheatmap(expression %>%
                                    cor() %>%
                                    as.matrix(),
                                    col = rev(RColorBrewer::brewer.pal(n = 8, name = "RdBu")),
                                    main = "Sample similarity",
                                    fontsize = 14)
    
    return(cowplot::plot_grid(sampleDis, sampleCor, labels = "AUTO"))
  }
}

