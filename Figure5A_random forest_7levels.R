library(reshape2)
library(dplyr)
library(phyloseq)
library(vegan)
library(ggplot2)
library(microbiome)
library(randomForest)
######################################################
#1.Random Forest with OOB-relative abundance
#Transform to relative abundance = RELA
l_rela_all<-lapply(l_phylo_all, function(x) transform_sample_counts(x, function(x){x / sum(x)}))
l_rela_rf_all <-l_rela_all
for (i in 1:length(l_rela_all)){ 
  predictors <- t(otu_table(l_rela_all[[i]]))
  response <- as.factor(sample_data(l_rela_all[[i]])[[l_category[i], exact = FALSE]])
  set.seed(123)
  l_rela_rf_all[[i]] <- randomForest(predictors, response, importance=TRUE, proximity=TRUE, ntree = 501, confusion = T, err.rate =T)
}

l_rela_rf_all[[37]]

saveRDS(l_rela_rf_all, "l_rela_rf_all.rds")
######################################################
#2.Random Forest with OOB-rar(rarefied)
l_rar_all<-lapply(l_phylo_all, function(x) rarefy_even_depth(x, 
                                                             sample.size = 5000, 
                                                             replace = FALSE, 
                                                             trimOTUs = TRUE, 
                                                             rngseed = 56987))
l_rar_rf_all <-l_rar_all
for (i in 1:length(l_rar_all)){ 
  predictors <- t(otu_table(l_rar_all[[i]]))
  response <- as.factor(sample_data(l_rar_all[[i]])[[l_category[i], exact = FALSE]])
  set.seed(123)
  l_rar_rf_all[[i]] <- randomForest(predictors, response, importance=TRUE, proximity=TRUE, ntree = 501, confusion = T, err.rate =T)
}
l_rar_rf_all[[42]]
saveRDS(l_rar_rf_all, "l_rar_rf_all_new.rds")
######################################################
#3.Random Forest with OOB-rra(rarefied+relative abundance =rra)
l_rra_all<-lapply(l_rar_all, function(x) transform_sample_counts(x, function(x){x / sum(x)})) 
l_rra_rf_all <-l_rra_all
for (i in 1:length(l_rra_all)){ 
  predictors <- t(otu_table(l_rra_all[[i]]))
  response <- as.factor(sample_data(l_rra_all[[i]])[[l_category[i], exact = FALSE]])
  set.seed(123)
  l_rra_rf_all[[i]] <- randomForest(predictors, response, importance=TRUE, proximity=TRUE, ntree = 501, confusion = T, err.rate =T)
}
l_rra_rf_all[[49]]
saveRDS(l_rra_rf_all, "l_rra_rf_all_new.rds")
#######################################################
#4.Random Forest with OOB-rar(CLR)
library(zCompositions)
library(CoDaSeq)
l_phylo_all <-lapply(l_phylo_all, function(x) prune_samples(sample_sums(x) > 0, x))
l_otu_all<-lapply(l_phylo_all, function(x) as.data.frame(otu_table(x)))
l_n0_all<-lapply(l_otu_all, function(x) cmultRepl(t(x), method="CZM", label=0))
l_clr_all<-lapply(l_n0_all, function(x) codaSeq.clr(x))
l_phylo_clr_all <- l_phylo_all
for(i in 1:length(l_phylo_clr_all)){
  otu_table(l_phylo_clr_all[[i]])<-otu_table(as.data.frame(l_clr_all[[i]]), taxa_are_rows = F)
}

# Set prunescale 
prunescale = 0.0001
# Prune out rare OTUs by mean relative abundance set by prunescale
l_phylo_clr.prune_all <-l_phylo_clr_all
for (i in 1:length(l_phylo_clr_all)){ 
  taxmean<- taxa_sums(l_rela_all[[i]])/nsamples(l_rela_all[[i]])
  l_phylo_clr.prune_all[[i]] <- prune_taxa(taxmean > prunescale, l_phylo_clr_all[[i]])
}
l_phylo_clr.prune_all[46]

for(i in 1:length(l_phylo_clr.prune_all)){
  if(is.null(sample_data(l_phylo_clr.prune_all[[i]])[[l_category[i], exact = FALSE]])){
    stop("Category is not match!")
  }else{print(i)}
}

l_clr_rf_all <-l_phylo_clr.prune_all
for (i in 1:length(l_phylo_clr.prune_all)){ 
  predictors <- otu_table(l_phylo_clr.prune_all[[i]])
  response <- as.factor(sample_data(l_phylo_clr.prune_all[[i]])[[l_category[i], exact = FALSE]])
  set.seed(123)
  l_clr_rf_all[[i]] <- randomForest(predictors, response, importance=TRUE, proximity=TRUE, ntree = 501, confusion = T, err.rate =T)
}
l_clr_rf_all[[37]]
saveRDS(l_clr_rf_all, "l_clr_rf_all_new.rds")

#Test
l_clr_rf_all[[41]]

