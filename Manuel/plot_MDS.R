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
#'@import tidyverse 
#'@import ggrepel
#'@example 

plot_MDS <-function(expression,group,point.size,
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
    #ggTheme(1) +
    labs(x = "Leading LogFC dim 1", y = "Leading LogFC dim 2", title = "MDS plot") +
    labs(shape=LegendName_Shape, col=LegendName_Color)+
    ggrepel::geom_text_repel(data = mds,aes(label = rownames(mds)),max.overlaps = Inf)
  return(plotmds)
    
}


# generated with the help of rnorm()
set.seed(78)
m  <- matrix(rnorm(100) , nrow = 20)
rownames(m) <-  paste0(rep("gene",20)) 
colnames(m) <- c("S1","S2","S3", "S4","S5")
group <- colnames(m)
plot_MDS(m, point.size = 5, group = group)






