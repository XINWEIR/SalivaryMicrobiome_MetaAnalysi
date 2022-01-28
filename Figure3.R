setwd("~/Documents/ADONIS/ADONIS-2206")
library(ggplot2)
library(viridis)
clr_euc <- readRDS("adonis_df_e_CLR.rds")
rar_bray <- readRDS("adonis_df_rar_b.rds")
rar_euc <- readRDS("adonis_df_rar_e.rds")
rar_wUF <- readRDS("adonis_df_rar_u.rds")
rra_bray <- readRDS("adonis_df_rra_b.rds")
rra_euc <- readRDS("adonis_df_rra_e.rds")
rra_wUF <- readRDS("adonis_df_tss_u.rds")
tss_bray <- readRDS("adonis_df_tss_b.rds")
tss_euc <- readRDS("adonis_df_tss_e.rds")
tss_wUF <- readRDS("adonis_df_tss_u.rds")
         




#OTU
norm_matrix <- rbind(tss_wUF, rra_wUF, tss_bray[43:49,], rra_bray[43:49,], tss_euc[43:49,], rra_euc[43:49,], clr_euc[43:49,])
norm_matrix$R2 <- round(as.numeric(norm_matrix$R2)*100, digits =3)
norm_matrix$Fmodel <- as.numeric(norm_matrix$Fmodel)

P_ADONIS_OTU <- 
ggplot(data=norm_matrix, aes(x=R2, y=Category, fill=Index)) + 
  theme_bw() +
  scale_fill_manual(values = c("#01A08A","#01CBAE","#FF2500", "#F98400","#F2AD00", "#5BBCD6", "#A1D8E7"))+
  geom_point(color ="black", size=4, shape =24, stroke =0.5) +
  labs(title="D) OTU level ", y= "")+
  theme(plot.title=element_text(size=15, face="bold"),
        axis.text.y=element_text(size=10, face = "bold"),
        legend.position = "top") +
  scale_x_continuous(breaks=seq(0, 12, 1))
P_ADONIS_OTU

#ADONIS L5
norm_matrix_l5 <- rbind(tss_bray[22:28,], rra_bray[22:28,], tss_euc[22:28,], rra_euc[22:28,], clr_euc[22:28,])
norm_matrix_l5$R2 <- round(as.numeric(norm_matrix_l5$R2)*100, digits =3)
norm_matrix_l5$Fmodel <- as.numeric(norm_matrix_l5$Fmodel)
write.table(norm_matrix_l5, file= "adonis-l5.txt", row.names = F, col.names = T, sep = "\t")

P_ADONIS_l5 <- 
  ggplot(data=norm_matrix_l5, aes(x=R2, y=Category, fill=Index)) + 
  theme_bw() +
  scale_fill_manual(values = c("#01A08A","#01CBAE","#FF2500", "#F98400","#F2AD00", "#5BBCD6"))+
  geom_point(color ="black", size=4, shape =24, stroke =0.5) +
  labs(title="A) Family Level ", y = "")+
  theme(plot.title=element_text(size=15, face="bold"),
        axis.text.y=element_text(size=10, face = "bold")) +
  scale_x_continuous(breaks=seq(0, 10, 0.5))

P_ADONIS_l5

#ADONIS L6
norm_matrix_l6 <- rbind(tss_bray[29:35,], rra_bray[29:35,], tss_euc[29:35,], rra_euc[29:35,], clr_euc[29:35,])
norm_matrix_l6$R2 <- round(as.numeric(norm_matrix_l6$R2)*100, digits =3)
norm_matrix_l6$Fmodel <- as.numeric(norm_matrix_l6$Fmodel)
write.table(norm_matrix_l6, file= "adonis-l6.txt", row.names = F, col.names = T, sep = "\t")

P_ADONIS_l6 <- 
    ggplot(data=norm_matrix_l6, aes(x=R2, y=Category, fill=Index)) + 
    theme_bw() +
    scale_fill_manual(values = c("#01A08A","#01CBAE","#FF2500", "#F98400","#F2AD00", "#5BBCD6"))+
    geom_point(color ="black", size=4, shape =24, stroke =0.5) +
    labs(title="B) Genus Level ", y = "")+
    theme(plot.title=element_text(size=15, face="bold"),
          axis.text.y=element_text(size=10, face = "bold")) +
    scale_x_continuous(breaks=seq(0, 10, 0.5))
P_ADONIS_l6

#ADONIS L7
norm_matrix_l7 <- rbind(tss_bray[36:42,], rra_bray[36:42,], tss_euc[36:42,], rra_euc[36:42,], clr_euc[36:42,])
norm_matrix_l7$R2 <- round(as.numeric(norm_matrix_l7$R2)*100, digits =3)
norm_matrix_l7$Fmodel <- as.numeric(norm_matrix_l7$Fmodel)
write.table(norm_matrix_l7, file= "adonis-l7.txt", row.names = F, col.names = T, sep = "\t")

P_ADONIS_l7 <- 
  ggplot(data=norm_matrix_l7, aes(x=R2, y=Category, fill=Index)) + 
  theme_bw() +
  scale_fill_manual(values = c("#01A08A","#01CBAE","#FF2500", "#F98400","#F2AD00", "#5BBCD6"))+
  geom_point(color ="black", size=4, shape =24, stroke =0.5) +
  labs(title="C) Species Level ", y= "")+
  theme(plot.title=element_text(size=15, face="bold"),
        axis.text.y=element_text(size=10, face = "bold")) +
  scale_x_continuous(breaks=seq(0, 10, 0.5))

P_ADONIS_l7

#######################################################################################
#Make Plot
#######################################################################################
library(cowplot)
plot_adonis_56 <- plot_grid(P_ADONIS_l5 + theme(legend.position ="none"), 
                         P_ADONIS_l6 + theme(legend.position ="none"), 
                         align = "h", nrow =1,  rel_widths = c(4:3))
plot_adonis_56
plot_adonis_7O <- plot_grid(P_ADONIS_l7 + theme(legend.position ="none"), 
                            P_ADONIS_OTU + theme(legend.position ="none"), 
                            align = "h", nrow =1,  rel_widths = c(4:3))
plot_adonis_7O
plot_adonis <-plot_grid(plot_adonis_56,plot_adonis_7O,align = "v", nrow =2)
plot_adonis
legend_p <-get_legend(P_ADONIS_OTU)
plot_adonis_legend <-plot_grid("",legend_p, "",ncol = 1, rel_heights = c(1, 1, 1))
plot_adonis_legend
plot_adonis_all <- plot_grid(plot_adonis, plot_adonis_legend, ncol = 2, rel_widths = c(5:1))
plot_adonis_all

#Save "Square"
ggsave(plot_adonis_all, file = "adonis-plot-2206_v1.pdf", width = 450, height = 220, unit = 'mm')
#Save version 2
ggsave(plot_adonis_all, file = "adonis-plot-2206_v2.pdf", width = 450, height = 160, unit = 'mm')


#Add PCoA Plot
setwd("~/Documents/ADONIS/ADONIS-2206/PCA-rds")
p_tss_bray_eth <- readRDS("p_tss_bray_eth.rds")
p_tss_bray_eth
p_tss_wUF_eth <- readRDS("p_tss_wUF_eth.rds")
p_tss_wUF_eth
p_PCA <- plot_grid(p_tss_bray_eth, p_tss_wUF_eth, ncol =1, rel_heights = c(1,1))
p_PCA
p_AF <- plot_grid(plot_adonis, p_PCA, ncol =2, rel_widths = c(1.8,1))
p_AF

