########################################################
#Differential pathways-Geographic location-KEGG_Level3
########################################################
library(dplyr)
HOMD_l6_eth <- readRDS("HOMD_l6_eth.rds") #genus table
meta_eth <- readRDS("meta_eth.rds")  #metadata
KEGG_l3 <- readRDS("KEGG_l3.rds")  #KEGG pathway table
rownames(L3_eth) <- KEGG_l3$Pathway

#ANCOMBC
library(phyloseq)
pq_L3_eth<-phyloseq(
  otu_table(L3_eth, taxa_are_rows = T),
  sample_data(meta_eth_5000)
)
library(ANCOMBC)
l3_eth_out = ancombc(phyloseq = pq_L3_eth,
                    formula = "ethnicity + hypervariable",
                    p_adj_method = "holm",
                    zero_cut = 0.9,
                    lib_cut = 0,
                    group = "ethnicity",
                    struc_zero = TRUE,
                    neg_lb = TRUE,
                    tol = 1e-05,
                    max_iter = 100,
                    conserve = TRUE,
                    alpha = 0.05,
                    global = FALSE
)
res_eth_l3 = l3_eth_out$res
diff_ethl3_groups = res_eth_l3$diff_abn
diff_ethl3_groups <-add_rownames(diff_ethl3_groups)
colnames(diff_ethl3_groups)[1] <- "CvsW"
diff_ethl3_w = res_eth_l3$W
diff_ethl3_w <-add_rownames(diff_ethl3_w)
colnames(diff_ethl3_w)[2] <- "W Statistic"
colnames(diff_ethl3_w)[1] <- "CvsW"
diff_ethl3_groups <- left_join(diff_ethl3_groups[1:2], diff_ethl3_w[,1:2])
#Differential pathways-ANCOMBC
Diff_eth_l3 <- subset(diff_ethl3_groups, ethnicityWestern==TRUE)
rownames(Diff_eth_l3) <- Diff_eth_l3$CvsW
colnames(Diff_eth_l3)[1] <- "Pathway"
Eth_l3_full <- left_join(Diff_eth_l3, KO_l1_l3)
ETH_l3_full <- subset(Eth_l3_full, PathwayL1 == "Metabolism")


#Random Forest
pq_L3_eth  = transform_sample_counts(pq_L3_eth, function(x) x / sum(x))
View(otu_table(pq_L3_eth))
response_eth_l3 <- as.factor(sample_data(pq_L3_eth)$ethnicity)
predictors_eth_l3 <- t(otu_table(pq_L3_eth))
#predictors <- data.frame(predictors)
dim(predictors_eth_l3)

library(randomForest)
library(dplyr)
set.seed(123)
rf_model_eth_l3 <- randomForest(predictors_eth_l3, response_eth_l3, ntree=501, importance=T, confusion=T, err.rate=T)
print(rf_model_eth_l3)
write.table(rf_model_eth_l6$confusion, file = "eth_l3_confusion.txt", sep = "\t", quote = F, row.names = T, col.names = T)
#Importance
varImpPlot(rf_model_eth_l3)
imp_eth_l3 = as.data.frame(round(importance(rf_model_eth_l3), 2))
imp_eth_l3=imp_eth_l3[order(imp_eth_l3$MeanDecreaseAccuracy,decreasing = F),]
#Cross-validation
result_eth_l3= replicate(5, rfcv(predictors_eth_l3, response_eth_l3, cv.fold = 5, scale = "log", step = 0.9), simplify=FALSE)
error.cv.eth.l3= sapply(result_eth_l3, "[[", "error.cv")
colnames(error.cv.eth.l3) = paste('err',1:5,sep='.')
error.cv.eth.l3 <-data.frame(error.cv.eth.l3, err.mean = apply(error.cv.eth.l3,1,mean))
allerr_eth_l3 <- add_rownames(data.frame(error.cv.eth.l3), var = "num")
allerr_eth_l3$num <- as.numeric(allerr_eth_l3$num)
#Quick visulisation
matplot(result_eth_l3[[1]]$n.var, cbind(rowMeans(error.cv.eth.l3), error.cv.eth.l3), type="l",
        lwd=c(2, rep(1, ncol(error.cv.eth.l3))), col=1, lty=1, log="x",
        xlab="Number of variables", ylab="CV Error")
#Define the optimal number of important feature
optimal= 100
library(ggplot2)
n_imp_eth_l3 = ggplot() + 
  #  main_theme +
  theme_classic() + 
  geom_line(aes(x = allerr_eth_l3$num, y = allerr_eth_l3$err.1), colour = 'grey') + 
  geom_line(aes(x = allerr_eth_l3$num, y = allerr_eth_l3$err.2), colour = 'grey') + 
  geom_line(aes(x = allerr_eth_l3$num, y = allerr_eth_l3$err.3), colour = 'grey') + 
  geom_line(aes(x = allerr_eth_l3$num, y = allerr_eth_l3$err.4), colour = 'grey') + 
  geom_line(aes(x = allerr_eth_l3$num, y = allerr_eth_l3$err.5), colour = 'grey') + 
  geom_line(aes(x = allerr_eth_l3$num, y = allerr_eth_l3$err.mean), colour = 'black') + 
  geom_vline(xintercept = optimal, colour='black', lwd=0.36, linetype="dashed") + 
  #  geom_hline(yintercept = min(allerr$err.mean), colour='black', lwd=0.36, linetype="dashed") + 
  coord_trans(x = "log2") +
  scale_x_continuous(breaks = c(1, 2, 5, 10, 20, 30, 50, 100, 200)) + # , max(allerr$num)
  labs(title=paste('Geographic Location (n = ', dim(predictors_eth_l3)[1],')', sep = ''), 
       x='Number of OTUs ', y='Cross-validation error rate') + 
  annotate("text", x = optimal+80, y = max(allerr_eth_l3$err.mean), label=paste("optimal = ", optimal, sep="")) 
n_imp_eth_l3
setwd("~/Documents/PICRUSt2/KEGG0923")
ggsave(n_imp_eth_l3, file = "Ethnicity_rfcv_plot_pathway3_new.pdf", width = 150, height = 75, unit = 'mm')
#Differential pathways-random forest
imp_eth_l3 = tail(imp_eth_l3, n = optimal)
imp_eth_l3 <- add_rownames(imp_eth_l3, var = "Pathway")
Imp_eth_l3 <- left_join(imp_eth_l3, KO_l1_l3)
#Select pathways related to Metabolism
IMP_l3_full <- subset(Imp_eth_l3, PathwayL1 == "Metabolism")

############################################
#Differential pathways-Western vs Chinese
############################################
DA_eth_l3 <- full_join(ETH_l3_full,IMP_l3_full)
DP_eth_l3 <- na.omit(DA_eth_l3)
DP_eth_l3 <- DP_eth_l3[-8,]
DP_eth_l3 = DP_eth_l3[order(DP_eth_l3$`W Statistic`,decreasing = F),]
DP_eth_l3[1:29,] = DP_eth_l3[1:29,][order(DP_eth_l3$`PathwayL2`[1:29],decreasing = F),]
#rownames(L3_eth)[5] <- "Valine, leucine and isoleucine degradation [PATH:ko00280]"
library(vegan)
L3_eth_rel <-decostand(L3_eth, method="total", MARGIN=2)
L3_eth_rel <- add_rownames(L3_eth_rel, var = "Pathway")
eth_info_l3 <-left_join(DP_eth_l3, L3_eth_rel)
eth_l3_rel <- eth_info_l3[, -c(1:9)]
rownames(eth_l3_rel) <- eth_info_l3$Pathway
rownames(eth_l3_rel) <- gsub("\\[.*\\]", "", rownames(eth_l3_rel))


############################################
#Differential genus-Western vs Chinese
############################################
setwd("~/Documents/PICRUSt2/KEGG")
DA_eth_l6_intersect <-readRDS("DA_eth_l6_intersect.rds")
DA_eth_l6_intersect = DA_eth_l6_intersect [order(DA_eth_l6_intersect$`W Statistic`,decreasing = T),]
library(vegan)
HOMD_l6_eth_rel <-decostand(HOMD_l6_eth, method="total", MARGIN=2)
DA_l6_eth_rel <- HOMD_l6_eth_rel[DA_eth_l6_intersect$taxa,]
DA_l6_eth_rel <- DA_l6_eth_rel[, intersect(colnames(KEGG_L3), colnames(HOMD_l6_eth))]
l6_eth_rel <- DA_l6_eth_rel


#######################
#Spearman correlation
######################
library(corrplot)
matrix_rel_eth <- cor(t(l6_eth_rel),t(eth_l3_rel),method="spearman",use="complete.obs")

library(psych)
metrix_p <-corr.test(t(l6_eth_rel),t(eth_l3_rel),method="spearman",adjust="fdr")
matrix_q <- data.frame(metrix_p["p.adj"])
colnames(matrix_q) <- colnames(matrix_rel_eth)
matrix_q <- as.matrix(matrix_q)
str(matrix_rel_eth)
str(matrix_q)
corrplot(matrix_rel_eth,
         #title = "Correlation Plot",
         method = "circle",
         outline = T, 
         addgrid.col = "grey", 
         p.mat= matrix_q,
         sig.level = 0.01,
         insig ="blank",
         #mar = c(0,0,1,0), 
         addrect = 4, 
         #rect.col = "grey", 
         #rect.lwd = 1,
         cl.pos = "b", 
#        tl.pos = "ld",
         tl.col = "black", #colour of text labels 
         tl.cex = 0.5, #Size of text labels
         cl.cex = 0.5)
write.table(DA_eth_l3, file = "DA_eth_l3.txt", sep = "\t", quote = F, row.names = F, col.names = T)




