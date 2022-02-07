############################
#RandomForest-important OTUs
############################
#This script were used four times for four Categories: Ethnicity/ Hypervariable Regions/ Sample Type/ Smoking
#Category:Ethnicity were used as example here
library(qiime2R)
library(phyloseq)
library(nloptr)
library(dplyr)

#Ethnicity-pq
setwd("~/Documents/OTUTABLE")
HOMD_eth <- read_qza("HOMD-ethnicity-table.qza")$data
setwd("~/Documents/metaRDS")
meta_eth <- readRDS("meta_eth.rds")

pq_eth<-phyloseq(
  otu_table(HOMD_eth, taxa_are_rows = T),
  sample_data(meta_eth)
)
pq_eth

#TSS transformation
pq_eth  = transform_sample_counts(pq_eth, function(x) x / sum(x) )
View(otu_table(pq_eth))

response_eth <- as.factor(sample_data(pq_eth)$ethnicity)
predictors_eth <- t(otu_table(pq_eth))
#predictors <- data.frame(predictors)
dim(predictors_eth)
########################
#train the model
#######################
library(randomForest)
library(dplyr)
set.seed(42)
rf_model_eth <- randomForest(predictors_eth, response_eth, ntree=501, importance=T, confusion=T, err.rate=T)
print(rf_model_eth)
write.table(rf_model_eth$confusion, file = "eth_otu_confusion.txt", sep = "\t", quote = F, row.names = T, col.names = T)
#Importance
varImpPlot(rf_model_eth)
imp_eth = as.data.frame(round(importance(rf_model_eth), 2))
imp_eth=imp_eth[order(imp_eth$MeanDecreaseAccuracy,decreasing = F),]
###################################################################
#Cross-validation
result_eth= replicate(5, rfcv(predictors_eth, response_eth, cv.fold = 5, scale = "log", step = 0.9), simplify=FALSE)
error.cv.eth= sapply(result_eth, "[[", "error.cv")
colnames(error.cv.eth) = paste('err',1:5,sep='.')
error.cv.eth <-data.frame(error.cv.eth, err.mean = apply(error.cv.eth,1,mean))
allerr_eth <- add_rownames(data.frame(error.cv.eth), var = "num")
allerr_eth$num <- as.numeric(allerr_eth$num)
#Quick visulisation
matplot(result_eth[[1]]$n.var, cbind(rowMeans(error.cv.eth), error.cv.eth), type="l",
        lwd=c(2, rep(1, ncol(error.cv.eth))), col=1, lty=1, log="x",
        xlab="Number of variables", ylab="CV Error")
#Define the optimal number of important feature (The first point that stop dropping)
optimal= 109
#ggplot2 visualise result of cross validation
library(ggplot2)
n_imp_eth = ggplot() + 
  #  main_theme +
  theme_classic() + 
  geom_line(aes(x = allerr_eth$num, y = allerr_eth$err.1), colour = 'grey') + 
  geom_line(aes(x = allerr_eth$num, y = allerr_eth$err.2), colour = 'grey') + 
  geom_line(aes(x = allerr_eth$num, y = allerr_eth$err.3), colour = 'grey') + 
  geom_line(aes(x = allerr_eth$num, y = allerr_eth$err.4), colour = 'grey') + 
  geom_line(aes(x = allerr_eth$num, y = allerr_eth$err.5), colour = 'grey') + 
  geom_line(aes(x = allerr_eth$num, y = allerr_eth$err.mean), colour = 'black') + 
  geom_vline(xintercept = optimal, colour='black', lwd=0.36, linetype="dashed") + 
  #  geom_hline(yintercept = min(allerr$err.mean), colour='black', lwd=0.36, linetype="dashed") + 
  coord_trans(x = "log2") +
  scale_x_continuous(breaks = c(1, 2, 5, 10, 20, 30, 50, 100, 200)) + # , max(allerr$num)
  labs(title=paste('Geographic Location (n = ', dim(predictors_eth)[1],')', sep = ''), 
       x='Number of OTUs ', y='Cross-validation error rate') + 
  annotate("text", x = optimal+60, y = max(allerr_eth$err.mean), label=paste("optimal = ", optimal, sep="")) 
n_imp_eth


imp_eth_tax <- left_join(imp_eth, HOMD.tax)
imp_eth_tax <- left_join(imp_eth_tax, Core.MRA)
colourfill <- data.frame(Core = c("core", "NA"), bdcol = c("black", "white"))
imp_eth_tax[is.na(imp_eth_tax)] <- "NA"
imp_eth_tax <- left_join(imp_eth_tax, colourfill)
imp_eth_tax$taxa = factor(imp_eth_tax$taxa, levels = imp_eth_tax$taxa)

#ggplot2-FigureS8
phylum_col =c("#7FB1D0", "#7CBEAE", "#F4AE6B", "#E1586C", "#F5E68B", "#06948E", "grey")
ggplot(imp_eth_tax, aes(x = MeanDecreaseAccuracy, y = taxa, fill = Phylum)) +
  geom_bar(stat = "identity", width=0.85, col = imp_eth_tax$bdcol) +
  scale_fill_manual(values = phylum_col)+
  labs(title=paste('Geographic Location - Top ', optimal,' important features', sep = ''), 
       x='Mean Decrease Accuracy ', y='OTU') + 
  #ylab(NULL) + #remove y axis title
  theme(text=element_text(family="sans", size=10, face = "bold"),
        axis.line = element_line(colour = "black", size = 0.5),
        #panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), #??????
        panel.background = element_blank(),
        legend.position ="none")




