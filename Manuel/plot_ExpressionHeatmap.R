#' Expression Heatmap
#' @param expTable The expression table to use to create the heatmap
#' @param Title Title of the plot
#' @param sampleDGroup A data frame specifying groups/factors for every sample
#' @param clColumns logical. If TRUE columns will be clustered
#' @param FontSRow Font size for the row labels
#' @param FontSColumn Font size for the column labels
#' @param FontSTitle Font size for the title
#' @param setseed Set seed to obtain the same colors
#' @import ComplexHeatmap
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


