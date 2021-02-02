#Test the functions
setwd("~/Documents/BioFig.git/Manuel/")
source("./R_Functions_Interest_JMMK.R")

data<-read.csv("./Data/Data_expression_normalized.csv",row.names = 1)

SampleInfo<-as.data.frame(cbind(file_name=c("Ctrl_1.CEL","Ctrl_2.CEL","Ctrl_3.CEL",
                                            "Dom_1.CEL","Dom_2.CEL","Dom_3.CEL",
                                            "N_D_3.CEL","N_D_4.CEL","N_D_5.CEL",
                                            "NEN_2.CEL","NEN_3.CEL","NEN_5.CEL"),
                                sample_name=c("Ctrl_1","Ctrl_2","Ctrl_3",
                                              "Dom_1","Dom_2","Dom_3",
                                              "N_D_3","N_D_4","N_D_5",
                                              "NEN_2","NEN_3","NEN_5"),
                                treatment = c(rep('Ctrl',3),rep('Dom',3),
                                              rep('N_D',3),rep('NEN',3)),
                                colors = c(rep('blue',3),rep('chartreuse4',3),
                                           rep('cyan',3),rep('purple',3))))

PCA_example<-plot_2DPCA(expression = data,
           group = SampleInfo$treatment,
           samplenames = SampleInfo$sample_name,
           title = "PCA",
           LegendName_Color = "Treatment",
           ggrepelLab = TRUE,size_title = 25)
ggsave(file=paste0("./Figures/PCA_example.svg"), plot=PCA_example, width=12, height=10)
ggsave(file=paste0("./Figures/PCA_example.pdf"), plot=PCA_example, width=12, height=10)


library(WGCNA)
library(gplots)
library(flashClust)
svg("./Figures/Dendrogram_example.svg",width = 8,height = 5)
Dendrogram_example<-sample_dendrogram(expression = data,
                                      colors = SampleInfo$colors,
                                      title = "Sample dendrogram",
                                      labels = c("Treatment"))
dev.off()


DomvsCtrl<-read.csv("./Data/DomperidonevsCtrlAllGenes.csv",row.names = 1)


