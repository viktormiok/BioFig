# Sample distance heatmap
#'@description This function created to show sample (dis) similarity in a heatmap using distance or correlation methods. 
#'@param expression Numerical matrix with samples as columns (e.g gene expression)
#'@import tidyverse 
#'@import pheatmap 
#'@import RColorBrewer
#'@example 


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




# generated with the help of rnorm()
set.seed(78)
m  <- matrix(rnorm(100) , nrow = 20)
rownames(m) <-  paste0(rep("gene",20)) 
colnames(m) <- c("S1","S2","S3", "S4","S5")
group <- colnames(m)
plot_SampleDistance(m)



