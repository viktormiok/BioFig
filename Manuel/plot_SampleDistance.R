plot_SampleDistance<-function(expression,...){
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
  listplots<-list()
  listplots[[1]]<-sampleDis
  listplots[[2]]<-sampleCor
  return(listplots)
}
