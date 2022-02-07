####################################
#Alluvial diagram
####################################
library(ggplot2)
library(ggalluvial)
library(viridis)
library(readxl)
connecttable <- read_excel("connecttable.xlsx")
connecttable=connecttable[order(connecttable$occ,decreasing = T),]
  
#White
alluvial_white <-ggplot(data = connecttable,
       aes(axis1 = Genus, axis2 = otu, 
           axis3 = zOTU,
           axis4 =zotuGenus,
           y = occ)) +
  scale_x_discrete(limits = c("Genus"
                              ,"otu"
                              , "zOTU"
                              , "zotuGenus"
                              ), expand = c(.1, .1)) +
  xlab("Demographic") +
  scale_fill_viridis(discrete=T)+
  geom_alluvium(aes(fill = zOTU)) +
  scale_y_continuous(expand=c(0.01,0.01))+
  geom_stratum(color = "grey") +
#  geom_flow(stat = "alluvium", lode.guidance = "rightleft", 
#            color = "darkgray") + 
  geom_text(stat = "stratum", size=3, aes(label = after_stat(stratum))) +
  theme(
    #text=element_text(family="sans", size=10, face = "bold"),
    axis.title = element_text(family="sans", size=10, face = "bold"),
    axis.text.y = element_text(family="sans", size=8, face = "bold"),
    axis.text.x = element_text(family="sans", size=8),
    #axis.line = element_line(colour = "black", size = 0.5),
    panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
    panel.background = element_blank(),
    legend.title = element_text(size = 8, face = "bold")
  )
alluvial_white

########################################
#Alpha diversity-zotu(same as HOMD-otu)
#######################################
library(ggplot2)
library(ggprism)
pq_rich_xav_z <- estimate_richness(pq_xav_z, split = TRUE, measures = NULL)
library(dplyr)
pq_rich_meta_xav_z <- left_join(add_rownames(pq_rich_xav_z, var="SampleID"), metadata.xav[,1:5], by= "SampleID")
#Shannon
pq_rich_arr_xav_z <- arrange(pq_rich_meta_xav_z, ethnicity)
res_z <- wilcox.test(Shannon ~ ethnicity, data = pq_rich_arr_xav_z,
                   exact = FALSE)
res_z
df_p_val_z <- data.frame(
  group1 = "Chinese",
  group2 = "Western",
  label = "p=0.073",
  y.position = 4.5
)
ad_s_z <-pq_rich_arr_xav_z %>%
  ggplot(aes(x=ethnicity, y=Shannon))+
  geom_boxplot(size=0.5, aes(fill = ethnicity))+
  scale_fill_manual(values = c("#57A6C5", "#D99787"))+
  labs(x="Geographic Location") +
  add_pvalue(df_p_val_z,
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
ad_s_z
#Chao1
res_C_z <- wilcox.test(Chao1 ~ ethnicity, data = pq_rich_arr_xav_z,
                     exact = FALSE)
res_C_z
df_p_val_c_z <- data.frame(
  group1 = "Chinese",
  group2 = "Western",
  label = "p<0.001",
  y.position = 410
)
ad_c_z <-pq_rich_arr_xav_z %>%
  ggplot(aes(x=ethnicity, y=Chao1))+
  geom_boxplot(size=0.5, aes(fill = ethnicity))+
  scale_fill_manual(values = c("#57A6C5", "#D99787"))+
  labs(x="Geographic Location") +
  add_pvalue(df_p_val_c_z,
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
ad_c_z
res_sp_z <- wilcox.test(Simpson ~ ethnicity, data = pq_rich_arr_xav_z,
                      exact = FALSE)
res_sp_z
library(cowplot)
library(ggplot2)
ad_z <- plot_grid(ad_s_z+ theme(legend.position="none"), 
                         ad_c_z + theme(legend.position="none"))
ad_z