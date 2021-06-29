plot_ExploratoryPlots<-function(expression,group,colors=NULL,shape=NULL,
                                samplenames,title="PCA",LegendName_Color="group",
                                LegendName_Shape="shape",LegendName="group",
                                ggrepelLab=TRUE,size_gglab=5,size_title=1,
                                point.size=4,scale=T,ntop=NULL,transform="no",
                                MahalanobisEllips=F, plottype = c("all","PCA", "MDS", "Distance","Expression"),...){
  
  
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
  
  
  
  switch(plottype,
         all={
           #PCA
           p_PCA<-plot_2DPCA(expression=expression,scl=scl,
                         group=group,shape=shape,samplenames=samplenames,
                         title = title,colors=colors,ggrepelLab = ggrepelLab,
                         size_gglab = size_gglab,size_title = size_title,
                         point.size = point.size,LegendName_Color = LegendName_Color,
                         LegendName_Shape = LegendName_Shape,LegendName = LegendName,
                         MahalanobisEllips = MahalanobisEllips,...)
           #MDS
           p_MDS<-plot_MDS(expression=expression,group=group)
           
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
                         title = title,colors=colors,ggrepelLab = ggrepelLab,
                         size_gglab = size_gglab,size_title = size_title,
                         point.size = point.size,LegendName_Color = LegendName_Color,
                         LegendName_Shape = LegendName_Shape,LegendName = LegendName,
                         MahalanobisEllips = MahalanobisEllips,...)
           plot(p)
           ggsave(file="PCA.svg", plot=p, width=10, height=8)
           ggsave(file="PCA.pdf", plot=p, width=10, height=8)
         },
         MDS={
           p<-plot_MDS(expression=expression,group=group)
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



# plots_arrange<-arrangeGrob(grobs=pList,ncol=2)
# ggsave(file="SampleDistances.svg", plot=plots_arrange, width=10, height=8)
# do.call("grid.arrange",c(pList,ncol=2))

