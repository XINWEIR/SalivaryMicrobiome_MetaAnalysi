library(tidyr)
library(reshape2)
library(vegan)
library(qiime2R)
library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
########################
#Formating Genus Names
########################
genus.name.list<- strsplit(DA_eth_l6_intersect$taxa,";")
genus.name <- vector(mode="numeric",length=0)
for(i in 1:48){
  genus.name[i]<- genus.name.list[[i]][6]
}
genus.name <- substring(genus.name, 4)
genus.name.df <- data.frame(name=genus.name, taxa=DA_eth_l6_intersect$taxa)
write.table(genus.name.df, file = "DA_eth_genus_names.txt",row.names=FALSE,col.names=TRUE, sep="\t")

DA_eth_genus_names <- read_excel("DA_eth_genus_names.xlsx")
DA_eth_l6_intersect <- left_join(DA_eth_l6_intersect, DA_eth_genus_names[,2:3])
########################
#Plot-ANCOMBC
########################
include_ABC_eth_l6 <- DA_eth_l6_intersect[,6:9]
colnames(include_ABC_eth_l6)[4] <- "taxa"
#spotplot for ANCOMBC
include_ABC_eth_l6$q_value = round(include_ABC_eth_l6$q_value, 2)
include_ABC_eth_l6=include_ABC_eth_l6[order(include_ABC_eth_l6$`W Statistic`,decreasing = F),]
include_ABC_eth_l6$Significance <- ifelse(include_ABC_eth_l6$q_value <= 0.01, "p<0.01", 
                                          ifelse(include_ABC_eth_l6$q_value <=0.05, "p<0.05", "p>0.05"))

include_ABC_eth_l6$taxa = factor(include_ABC_eth_l6$taxa, levels = include_ABC_eth_l6$taxa)
levels(include_ABC_eth_l6$taxa)
include_ABC_eth_l6$Significance = factor(include_ABC_eth_l6$Significance, levels = c("p>0.05", "p<0.05", "p<0.01"))
levels(include_ABC_eth_l6$Significance)
#Making plot-1
ABC_W <- 
  ggplot(include_ABC_eth_l6, aes(y=`W Statistic`, x=taxa, fill= Significance))+
  geom_point(size = 4, shape = 21, col = "black", stroke = 0.8) +
  #scale_color_viridis(discrete=TRUE) +
  scale_fill_manual(values = c("#FEDCD2", "#D99787"))+
  labs(y='W Statistic', x=NULL) +
  coord_flip() + 
  geom_hline(yintercept=0, linetype="dashed") +
  scale_y_continuous(breaks = seq(-18,10,4))+
  theme(
    #        text=element_text(family="sans", size=10, face = "bold"),
    axis.title = element_text(family="sans", size=10, face = "bold"),
    axis.text.y = element_text(family="sans", size=8, face = "bold"),
    axis.text.x = element_text(family="sans", size=8),
    axis.line = element_line(colour = "black", size = 0.5),
    #        panel.grid.major.x = element_line(color = "grey", size = 0.2),
    panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
    panel.background = element_blank(),
    legend.title = element_text(size = 8, face = "bold"),
    legend.key = element_rect(color = "white")
  )
ABC_W

###########################
#Plot-relative abundance
###########################
rel_eth_l6 <-decostand(HOMD_l6_eth, method="total", MARGIN=2)
rel_eth_l6_DA <- rel_eth_l6[DA_eth_l6_intersect$taxa,]
rel_eth_l6_DA <- add_rownames(data.frame(rel_eth_l6_DA), var = "taxa")
test <- left_join(rel_eth_l6_DA, DA_eth_genus_names[,2:3])
test <- test[,-1]
rownames(test) <- test$name2
rel_eth_l6_DA <- as.matrix(test)[, -1990]
dim(rel_eth_l6_DA)
rel_eth_l6_DA <- add_rownames(data.frame(t(rel_eth_l6_DA)), var="SampleID")
meta_eth_1 <- add_rownames(data.frame(meta_eth), var="SampleID")
rel_eth_l6_DA <- left_join(rel_eth_l6_DA, meta_eth_1[, c(1,5)])

l6_eth_DA_box <- rel_eth_l6_DA %>%
  gather(taxa, relabun, -SampleID, -ethnicity)

l6_eth_DA_box$taxa <- factor(l6_eth_DA_box$taxa, levels = levels(include_ABC_eth_l6$taxa))
levels(l6_eth_DA_box$taxa)
l6_eth_DA_box$relabun <- as.numeric(l6_eth_DA_box$relabun)
#Making plot-2
require(scales)
legend_title_relabun <- "Geographic location"
plot_rel_eth_l6 <- 
ggplot(l6_eth_DA_box,aes(x = relabun, y = taxa)) +
  geom_boxplot(aes(fill = ethnicity), position = position_dodge(preserve = "single")) +
  scale_fill_manual(legend_title_relabun, values = c("#57A6C5", "#D99787"))+
  labs(x='Relative Abundance (log10 scale)', y=NULL) + 
  #ylab(NULL) + #remove y axis title
  theme(
    #    text=element_text(family="sans", size=10, face = "bold"),
    axis.title = element_text(family="sans", size=10, face = "bold"),
    axis.text.y = element_text(family="sans", size=8, face = "bold"),
    axis.text.x = element_text(family="sans", size=8),
    axis.line = element_line(colour = "black", size = 0.5),
    panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
    #panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), #??????
    panel.background = element_blank(),
    #       legend.position ="none",
    legend.title = element_text(size = 8, face = "bold"),
    legend.key = element_rect(color = "white"))+
  scale_x_continuous(trans='log10', 
                     limits = c(10^-5, 0.5),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))

plot_rel_eth_l6


###########################
#Plot-random forest
###########################
include_imp_eth_l6 <- DA_eth_l6_intersect[, c(2,3,4,5,9)]
colnames(include_imp_eth_l6)[5] <- "taxa"
include_imp_eth_l6$taxa <- factor(include_imp_eth_l6$taxa, levels = levels(include_ABC_eth_l6$taxa))
levels(include_imp_eth_l6$taxa)
#Making plot
RF_imp <-
  ggplot(include_imp_eth_l6, aes(x = MeanDecreaseAccuracy, y = taxa)) +
  geom_bar(stat = "identity", fill = "#57A6C5",width=0.85, col = "black") +
  #  scale_fill_manual("#57A6C5")+
  labs(x='Mean Decrease Accuracy', y=NULL) + 
  #ylab(NULL) + #remove y axis title
  theme(
    #text=element_text(family="sans", size=10, face = "bold"),
    axis.title = element_text(family="sans", size=10, face = "bold"),
    axis.text.y = element_text(family="sans", size=8, face = "bold"),
    axis.text.x = element_text(family="sans", size=8),
    axis.line = element_line(colour = "black", size = 0.5),
    panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
    panel.background = element_blank(),
    legend.title = element_text(size = 8, face = "bold")
  )
RF_imp

###########################
#Plot-combine
###########################
library(cowplot)
legendA <- get_legend(ABC_W)
legendB <- get_legend(plot_rel_eth_l6)
legendC <- get_legend(RF_imp)

legendABC <- plot_grid(legendA, legendB, legendC, nrow = 3, align = "vh", rel_heights =c(1,1,0.8))
legendABC
eth_l6_plot <- plot_grid(ABC_W + theme(legend.position="none"), 
                         plot_rel_eth_l6 + theme(legend.position="none", axis.text.y = element_blank()), 
                         RF_imp + theme(legend.position="none", axis.text.y = element_blank()), 
                         align = "h", 
                         ncol =3, rel_widths = c(5,4,3.5))
eth_l6_plot
setwd("~/Documents/Official/Figure7/D_DA_ethnicity")
saveRDS(eth_l6_plot, "eth_l6_plot_nolegend.rds")
ggsave(eth_l6_plot, file = "Ethnicity_DA_level6_sizeL.pdf", width = 230, height = 200, unit = 'mm')





