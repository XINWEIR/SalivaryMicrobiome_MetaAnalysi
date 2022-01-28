library(qiime2R)
library(biomformat)
library(tidyverse)
library(reshape2)
library(vegan)
library(ggplot2)
######################
#Taxonomy
######################
#Colour
fill_col <- c("#498EAF","#E5BB4B","#D66C44", "#499360", "#631F16", "#E29E93","#E4DACE", "grey")
phylum_col =c("#D18237", "#9D1E31", "#4A746A", "#324856", "#CA4026", "#F5E68B", "#D1825B", "grey")

#Taxonomy-plot
library(readxl)
setwd("~/Documents")
Taxonomy_studyname <- read_excel("Taxonomy-studyname.xlsx")
Taxonomy_studyname$Taxonomy <-as.factor(Taxonomy_studyname$Taxonomy)
Taxonomy_studyname$study <-as.factor(Taxonomy_studyname$study)

taxa1 <-
Taxonomy_studyname %>%
  mutate(study= factor(study, levels= order2)) %>%
  ggplot(aes(x=value, y =study , fill = Taxonomy)) +
  geom_bar(stat = "identity",position="fill", width=0.85, col = "black")+ 
  scale_x_continuous(labels = scales::percent, expand = c(0,0.01)) +
  scale_fill_manual(values = phylum_col)+ 
  labs(x="Relative abundance(%)", y = NULL)+
  ggtitle("A)")+
  theme(text=element_text(family="sans", size=10, face = "bold"),
       legend.position ="none",
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid") #??????
        )
taxa1

#########################################################################
#Alpha Diversity("phyloseq" + "amplicon")
library("phyloseq")

otu_HOMD <- otu_HOMD[, intersect(rownames(metadata), colnames(otu_HOMD))]
#intersect(rownames(metadata), colnames(otu_HOMD))
#setdiff(rownames(meta_2267), colnames(otu_HOMD))
pq <- phyloseq(
  otu_table(otu_HOMD, taxa_are_rows = T),
  sample_data(metadata)
)
#Do not need to transform or filter before estmate_richness
#pq_rela <- transform_sample_counts(pq, function(x) {x/sum(x)})
View(otu_table(pq_rela))
pq_richness <- estimate_richness(pq, split = TRUE, measures = NULL)

###################################
#Group by ethnicity, sequence=study
###################################
pq_richness_meta <- left_join(add_rownames(pq_richness, var="SampleID"), add_rownames(metadata, var="SampleID"), by= "SampleID")
library("dplyr")
pq_richness_arrange <- arrange(pq_richness_meta, ethnicity,study)
order <- unique(pq_richness_arrange$study)
ap1 <-pq_richness_arrange %>%
  mutate(study= factor(study, levels= order)) %>%
  ggplot(aes(x=Shannon, y=study, fill=ethnicity))+
  geom_boxplot(size=0.5, col= "black")+
  scale_fill_manual(values = fill_col)+
  ggtitle("B)")+
  theme_classic()+
  theme(axis.text.y = element_blank(),axis.text.x = element_text(size = 10,face = "bold"))+
  theme(axis.text.x = element_text(hjust = 1))+
  theme(axis.title.x=element_text(size=11))+
  theme(axis.title.y=element_text(size=11))+
  theme(legend.position ="none") +
  labs(x="Shannon", y=NULL)+
  theme(plot.title=element_text(hjust=0.5), 
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        legend.title=element_blank())
ap1
ap2 <-pq_richness_arrange %>%
  mutate(study= factor(study, levels= order)) %>%
  ggplot(aes(x=Chao1, y=study, fill=ethnicity))+
  geom_boxplot(size=0.5, col= "black")+
  scale_fill_manual(values = fill_col)+
  ggtitle("C)")+
  theme_classic()+
  theme(axis.text.y = element_blank(),axis.text.x = element_text(size = 10,face = "bold"))+
  theme(axis.text.x = element_text(hjust = 1))+
  theme(axis.title.x=element_text(size=11))+
  theme(axis.title.y=element_text(size=11))+
  labs(x="Chao1", y = NULL)+
  theme(legend.position ="none") +
  theme(plot.title=element_text(hjust=0.5), 
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        legend.title=element_blank())
ap2
ap4 <-pq_richness_arrange %>%
  mutate(study= factor(study, levels= order)) %>%
  ggplot(aes(x=Simpson, y=study, fill=ethnicity))+
  geom_boxplot(size=0.5, col= "black")+
  scale_fill_manual(values = fill_col)+
  ggtitle("D)")+
  theme_classic()+
  theme(axis.text.y = element_blank(),axis.text.x = element_text(size = 10,face = "bold"))+
  theme(axis.text.x = element_text(hjust = 1))+
  theme(axis.title.x=element_text(size=11))+
  theme(axis.title.y=element_text(size=11))+
  labs(x="Simpson", y = NULL)+
  theme(legend.position ="none") +
  theme(plot.title=element_text(hjust=0.5), 
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        legend.title=element_blank())
ap4

#######################################################################################
#Make Plot
#######################################################################################
#get legends
legend_ap2 = get_legend(ap2)
legend_tax1 = get_legend(taxa1)
plot_grid(legend_tax1, legend_ap2, ncol =1, align = "v")

#Paste plots together
library(cowplot)
plot_3 <- plot_grid(taxa1, ap1, ap2, ap4, "",align = "h", ncol =5, rel_widths = c(3,1,1,1,1.5))
plot_3