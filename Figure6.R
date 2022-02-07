library(metacoder)
library(vegan)
library(taxa)
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(ape)
library(qiime2R)
#############################
#A:Data preparation
#############################
otu <- read_qza("HOMD-2206.qza")$data

#metadata
sample_data <- readRDS("meta_eth_2206.rds")
sample_data$ethnicity <- factor(sample_data$ethnicity, levels = c("China", "NorthAmerica", "Europe"), ordered = TRUE)
levels(sample_data$ethnicity)

#otu table
otu_data <- otu[, sample_data$SampleID]
dim(otu_data)
otu_data <- add_rownames(data.frame(otu_data), var = "OTU ID")

#taxonomy
setwd("~/Documents/2206-otu")
HOMD_taxonomy<-read_qza("ref-1522-taxonomy.qza")$data
HOMD_taxonomy_sep<-parse_taxonomy(HOMD_taxonomy)
colnames(HOMD_taxonomy)[1] <- "OTU ID"
HOMD_taxonomy_sep <- add_rownames(HOMD_taxonomy_sep, var= "OTU ID")
taxa_data <- left_join(HOMD_taxonomy, HOMD_taxonomy_sep, by=c('OTU ID'))

#Combine
taxa_data$`OTU ID` <- as.character(taxa_data$`OTU ID`) # Must be same type for join to work
otu_data$`OTU ID` <- as.character(otu_data$`OTU ID`) # Must be same type for join to work
otu_data <- left_join(otu_data, taxa_data,
                      by = c("OTU ID")) # identifies cols with shared IDs

#Making obj
obj_eth <- parse_tax_data(otu_data,
                          class_cols = "Taxon",
                          class_sep = ";",
                          class_regex = "^([a-z]{0,1})_{0,2}(.*)$",
                          class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))
###############################
#Metacoder process
###############################
print(obj_eth)
names(obj_eth$data) <- "otu_counts"
obj_eth$data$otu_props <- calc_obs_props(obj_eth, "otu_counts", other_cols = TRUE)
print(obj_eth$data$otu_counts)
print(obj_eth$data$otu_props)
obj_eth$data$tax_abund <- calc_taxon_abund(obj_eth, "otu_props")

obj_eth$data$diff_table <- compare_groups(obj_eth, data = "tax_abund",
                                          cols = sample_data$SampleID,
                                          groups = sample_data$ethnicity)
print(obj_eth$data$diff_table)
obj_eth <- mutate_obs(obj_eth, "diff_table",
                      wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))
obj_eth$data$diff_table$log2_median_ratio[obj_eth$data$diff_table$wilcox_p_value > 0.05] <- 0
set.seed(1)
obj_eth
obj_eth %>%
  metacoder::filter_taxa(taxon_ranks == "s", supertaxa = TRUE, reassign_obs = FALSE) %>%
  mutate_obs("cleaned_names", gsub(taxon_names, pattern = "\\[|\\]", replacement = "")) %>%
  metacoder::filter_taxa(grepl(cleaned_names, pattern = "^[a-zA-Z]+$"), reassign_obs = FALSE) %>%
  heat_tree_matrix(data = "diff_table",
                   node_label = cleaned_names,
                   node_size = n_obs, # number of OTUs
                   node_color = log2_median_ratio, # difference between groups
                   node_color_trans = "linear",
                   node_color_interval = c(-3, 3), # symmetric interval
                   edge_color_interval = c(-3, 3), # symmetric interval
                   node_color_range =  c("#BC5F6A", "gray", "#19B3B1"), # diverging colors
                   node_size_axis_label = "OTU count",
                   node_color_axis_label = "Log 2 ratio of median counts",
                   row_label_color = "#19B3B1",
                   col_label_color = "#BC5F6A",
                   layout = "da", initial_layout = "re",
                   key_size = 0.67,
                   seed = 2)

##############################
#BC-alpha diversity
##############################
meta_eth <- readRDS("meta_eth_2206.rds")
otu_eth <- data.frame(otu[,rownames(meta_eth)])
library(phyloseq)
pq_eth <- phyloseq(
  otu_table(otu_eth, taxa_are_rows = T),
  sample_data(meta_eth)
)
pq_richness <- estimate_richness(pq_eth, split = TRUE, measures = NULL)
library(dplyr)
pq_richness_meta <- left_join(add_rownames(pq_richness, var="SampleID"), meta_eth[,1:5], by= "SampleID")
pq_richness_arrange <- arrange(pq_richness_meta, ethnicity)
res <- wilcox.test(Shannon ~ ethnicity, data = pq_richness_arrange,
                   exact = FALSE)
res_C <- wilcox.test(Chao1 ~ ethnicity, data = pq_richness_arrange,
                     exact = FALSE)

library(ggplot2)
library(ggprism)
df_p_val <- data.frame(
  group1 = "Chinese",
  group2 = "Western",
  label = "p<0.001",
  y.position = 5.2
)
ad_s <-pq_richness_arrange %>%
  ggplot(aes(x=ethnicity, y=Shannon))+
  geom_boxplot(size=0.5, aes(fill = ethnicity))+
  scale_fill_manual(values = c("#57A6C5", "#D99787"))+
  labs(x="Geographic Location") +
  add_pvalue(df_p_val,
             xmin = "group1",
             xmax = "group2",
             label = "label",
             y.position = "y.position")+
  theme(
    #        text=element_text(family="sans", size=10, face = "bold"),
    axis.title = element_text(family="sans", size=10, face = "bold"),
    axis.text.y = element_text(family="sans", size=10, face = "bold"),
    axis.text.x = element_text(family="sans", size=10),
    #axis.line = element_line(colour = "black", size = 0.5),
    panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
    #        panel.grid.major.x = element_line(color = "grey", size = 0.2),
    panel.background = element_blank(),
    legend.title = element_text(size = 8, face = "bold"),
    legend.key = element_rect(color = "white"))
ad_s


df_p_val_c <- data.frame(
  group1 = "Chinese",
  group2 = "Western",
  label = "p<0.001",
  y.position = 600
)
ad_c <-pq_richness_arrange %>%
  ggplot(aes(x=ethnicity, y=Chao1))+
  geom_boxplot(size=0.5, aes(fill = ethnicity))+
  scale_fill_manual(values = c("#57A6C5", "#D99787"))+
  labs(x="Geographic Location") +
  add_pvalue(df_p_val_c,
             xmin = "group1",
             xmax = "group2",
             label = "label",
             y.position = "y.position")+
  theme(
    #        text=element_text(family="sans", size=10, face = "bold"),
    axis.title = element_text(family="sans", size=10, face = "bold"),
    axis.text.y = element_text(family="sans", size=10, face = "bold"),
    axis.text.x = element_text(family="sans", size=10),
    #axis.line = element_line(colour = "black", size = 0.5),
    #        panel.grid.major.x = element_line(color = "grey", size = 0.2),
    panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
    panel.background = element_blank(),
    legend.title = element_text(size = 8, face = "bold"),
    legend.key = element_rect(color = "white"))
ad_c
ad_cw <- plot_grid(ad_c+ theme(legend.position ="none"),
                   ad_s+ theme(legend.position ="none"),  
                   align = "h", nrow =1,  rel_widths = c(1:1))
ad_cw
ggsave(ad_cw, file = "Ethnicity_cvsw.pdf", width = 150, height = 75, unit = 'mm')




