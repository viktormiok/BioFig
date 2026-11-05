################################################################################
#                                                                              #
#          Load the libraries and data sets                                    #    
#                                                                              #
################################################################################
rm(list = ls())

library(tidyverse)               
library(phyloseq)                    
library(DESeq2)                        
library(microbiome)               
library(vegan)                          
library(picante)                       
library(ALDEx2)                        
library(metagenomeSeq)          
library(HMP)                              
library(dendextend)                
library(selbal)                          
#library(rms)
library(breakaway)      

################################################################################
# Test reading data from raw container
#cont_curated <- storage_container(ad_endp_tok, "raw")
#list_storage_files(cont_curated, dir="rni/HM_STUDIES", recursive=TRUE)
# Load clinical data at 1 year
#physeq <- storage_read_csv(cont_curated, "rni/HM_STUDIES/count_abundance.csv")

physeq <- read.csv("C:/Users/User/Documents/Rugby/count_abundance.csv")

# remove the subjects with low quality 
drops <- c("LA1","LA2","AM1","AM2")
physeq <- physeq[ , !(names(physeq) %in% drops)]

################## test ########################################################
#cont_curated <- storage_container(ad_endp_tok, "curated")
#list_storage_files(cont_curated, dir="rni/HM_STUDIES", recursive=TRUE)
# Load clinical data at 1 year
#ps <- storage_load_rds(cont_curated, "rni/HM_STUDIES/ps_giloteaux_2016.rds")
################################################################################

################################################################################
#                                                                              #
#          Create a phyloseq object                                            #    
#                                                                              #
################################################################################

taxmat <- as.matrix(physeq[,2:10])
otumat <- as.matrix(physeq[,11:45])
rownames(taxmat) <-rownames(otumat) <- taxmat[,1]
taxmat <- taxmat[,-1]

p <- colnames(otumat)
sdmat <- as.data.frame(as.numeric(substr(p, nchar(p), nchar(p))),p)
rownames(sdmat) <- p
#sdmat <- as.matrix(sdmat)
colnames(sdmat) <- "timepoint"
sdmat$subject <- substr(p, 1, 2)
sdmat$timepoint <- as.factor(sdmat$timepoint)

sdmat <- sdmat %>%
  mutate(treatment = case_when(
    timepoint == "1" ~ "Pre",
    timepoint == "2" ~ "Post",
    TRUE ~ timepoint  # Keep other values unchanged
  ))

sdmat <- sdmat %>%
  mutate(subjectID = case_when(
    subject == "CC" ~ "S1",
    subject == "CS" ~ "S2",
    subject == "DB" ~ "S3",
    subject == "DP" ~ "S4",
    subject == "DS" ~ "S5",
    subject == "GR" ~ "S6",
    subject == "JS" ~ "S7",
    subject == "KA" ~ "S8",
    subject == "LM" ~ "S9",
    subject == "LN" ~ "S10",
    subject == "MD" ~ "S11",
    subject == "MI" ~ "S12",
    subject == "ON" ~ "S13",
    subject == "RM" ~ "S14",
    subject == "SV" ~ "S15",
    subject == "VN" ~ "S16",
    subject == "VP" ~ "S17",
    subject == "WR" ~ "S18",
    TRUE ~ subject  # Keep other values unchanged
  ))

#sdmat$names <- rownames(sdmat)
#rownames(sdmat) <- sdmat$subjectID 


OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
SD = sample_data(sdmat)


ps = phyloseq(OTU, SD, TAX)
ps

################################################################################
#                                                                              #
#          Perform visualisations                                              #    
#                                                                              #
################################################################################

#Get count of phyla
table(phyloseq::tax_table(ps)[, "Phylum"])

sample_data(ps)$treatment <- factor(
  sample_data(ps)$treatment, 
  levels = c("Pre", "Post")
)

#Plot
phyloseq::plot_bar(ps, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "subject", y = "Relative Abundance\n") +
  facet_wrap(~ treatment, scales = "free") +
  scale_x_discrete(labels = ps@sam_data$subjectID) # +
  #theme(panel.background = element_blank(),
  #      axis.text.x=element_blank(),
  #      axis.ticks.x=element_blank())

#Agglomerate to phylum-level and rename
ps_phylum <- phyloseq::tax_glom(ps, "Phylum")
phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
phyloseq::otu_table(ps_phylum)[1:5, 1:5]

#Melt and plot
phyloseq::psmelt(ps_phylum) %>%
  ggplot(data = ., aes(x = treatment, y = Abundance, fill = treatment)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")


#Extract OTU table and compute BC
ps_rel_otu <- data.frame(phyloseq::otu_table(ps))
ps_rel_otu <- t(ps_rel_otu)
bc_dist <- vegan::vegdist(ps_rel_otu, method = "bray")
as.matrix(bc_dist)[1:5, 1:5]


#Save as dendrogram
ward <- as.dendrogram(hclust(bc_dist, method = "ward.D2"))
#Provide color codes
meta <- data.frame(phyloseq::sample_data(ps))
colorCode <- c(Control = "red", `Chronic Fatigue` = "blue")
labels_colors(ward) <- colorCode[meta$Status][order.dendrogram(ward)]
#Plot
plot(ward)

################################################################################
#                                                                              #
#          Beta diversity                                                      #    
#                                                                              #
################################################################################

#CLR transform
(ps_clr <- microbiome::transform(ps, "clr"))   

phyloseq::otu_table(ps)[1:5, 1:5]

phyloseq::otu_table(ps_clr)[1:5, 1:5]

#PCA via phyloseq
ord_clr <- phyloseq::ordinate(ps_clr, "RDA")
#Plot scree plot
phyloseq::plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

#Examine eigenvalues and % prop. variance explained
head(ord_clr$CA$eig)  

sapply(ord_clr$CA$eig[1:5], function(x) x / sum(ord_clr$CA$eig)) 

#Scale axes and plot ordination
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
phyloseq::plot_ordination(ps, ord_clr, type="samples", color="treatment", label="subjectID") + 
  geom_point(size = 2) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = treatment), linetype = 2)

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(ps_clr, method = "jaccard") # 
#ADONIS test
vegan::adonis2(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$treatment)

#Dispersion test and plot
dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(ps_clr)$timepoint)
dispr


plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Jaccard Distance", sub = "")

permutest(dispr)

boxplot(dispr, main = "", xlab = "")
################################################################################
#                                                                              #
#          Differential abundance testing                                      #    
#                                                                              #
################################################################################


#Generate data.frame with OTUs and metadata
ps_wilcox <- data.frame(t(data.frame(phyloseq::otu_table(ps_clr))))
ps_wilcox$timepoint <- phyloseq::sample_data(ps_clr)$timepoint
#Define functions to pass to map
wilcox_model <- function(df){
  wilcox.test(abund ~ timepoint, data = df) #wilcox
}
wilcox_pval <- function(df){
  wilcox.test(abund ~ timepoint, data = df)$p.value
}
#Create nested data frames by OTU and loop over each using map 
wilcox_results <- ps_wilcox %>%
  gather(key = OTU, value = abund, -timepoint) %>%
  group_by(OTU) %>%
  nest() %>%
  mutate(wilcox_test = map(data, wilcox_model),
         p_value = map(data, wilcox_pval))                       
#Show results
head(wilcox_results)

wilcox_results$wilcox_test[[1]]

wilcox_results$p_value[[1]]


wilcox_results <- wilcox_results %>%
  dplyr::select(OTU, p_value) %>%
  unnest()

head(wilcox_results)

fin <- as.data.frame(wilcox_results)

fin <- fin[order(fin$p_value, decreasing = F), ]
head(fin)

write.csv(fin,"C:/Users/Miokvi/OneDrive - Danone/Documents/Expertise/Rugby/Differential_abundance_1.csv", row.names = TRUE)

#Adding taxonomic labels
taxa_info <- data.frame(tax_table(ps_clr))
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")
#Computing FDR corrected p-values
wilcox_res <- wilcox_results %>%
  full_join(taxa_info) %>%
  arrange(p_value) %>%
  mutate(BH_FDR = p.adjust(p_value, "fdr")) %>%
  filter(BH_FDR < 1) %>%
  dplyr::select(OTU, p_value, BH_FDR, everything())

#Printing results
print.data.frame(wilcox_res[,1:2]) 

boxplot(wilcox_res[,2], main = "P-values of all bacterias, red line: p-value=0.05", ylab = "p-values")
abline(h=0.05, col="red")


################################################################################
#   Alpha diversity
################################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)

setwd("C:/Users/User/Documents/Rugby")

adiv <- read.table(file = 'alpha.tsv', sep = '\t', header = TRUE)

adiv$Cohort <- as.character(adiv$Cohort)
adiv <- adiv %>%
  mutate(Treatment = case_when(
    Cohort == "1" ~ "Pre",
    Cohort == "2" ~ "Post",
    TRUE ~ Cohort  # Keep other values unchanged
  ))

adiv$Treatment <- factor(
  adiv$Treatment, 
  levels = c("Pre", "Post")
)

colnames(adiv) <- c("Name","Cohort","Chao1","Simpson","Shannon","Treatment")

adiv %>%
  gather(key = metric, value = value, c("Chao1", "Simpson", "Shannon")) %>%
  mutate(metric = factor(metric, levels = c("Chao1", "Simpson", "Shannon"))) %>%
  ggplot(aes(x = Treatment, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Treatment), height = 0, width = .2) +
  labs(x = "Treatment", y = "Diversity index values") +
  facet_wrap(~ metric, scales = "free") #+
  #theme(legend.position="none")

#Summarize
adiv %>%
  group_by(Cohort) %>%
  dplyr::summarise(median_observed = median(chao1),
                   median_shannon = median(simpson),
                   median_pd = median(shannon))

#Wilcoxon test of location
wilcox.test(chao1 ~ Cohort, data = adiv, exact = FALSE, conf.int = TRUE)
t.test(simpson ~ Cohort, data = adiv, exact = FALSE, conf.int = TRUE)
t.test(shannon ~ Cohort, data = adiv, exact = FALSE, conf.int = TRUE)
