plot_2DPCA<-function(expression,group,colors=NULL,shape=NULL,
                     samplenames,title="PCA",LegendName_Color="group",
                     LegendName_Shape="shape",LegendName="group",
                     ggrepelLab=TRUE,size_gglab=5,size_title=14,
                     point.size=4,scl=T,ntop=NULL,
                     transform=c("no","vst","rlog"),scale=F,MahalanobisEllips=F){
  
  #Include check data type if not stop
  #(Viktorian)
  
  #Deseq2 option
  if(class(expression)=="DESeqDataSet"){
    group<-colData(expression)[,group]
    samplenames<-colData(expression)[,samplenames]
    expression<-assay(expression)
  }
  #EdgeR option
  if(class(expression)=="DGEList" && attr(attributes(expression)[[1]],"package")=="edgeR"){
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