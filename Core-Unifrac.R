#Modified from Shade & Stopnisek(https://github.com/ShadeLab/PAPER_Shade_CurrOpinMicro/blob/master/script/sncm.fit.R)
#          and unifrac() from package:ampvis2 (https://rdrr.io/github/MadsAlbertsen/ampvis2/src/R/internals.R)
library(foreach)
library(tidyverse)
library(reshape2)
library(vegan)
library(ape)
library(qiime2R)
library(dbplyr)
library(tidyr)
library(dplyr)
library(phyloseq)

nReads=5000
wUFaddition <- NULL
tree<-read_qza("HOMD_V15.1_mod2_rooted_tree_0222.qza")$data
otu <- readRDS("OTU__2206_rar.rds")
###########################
#Rank
###########################
head(otu[,1:10])
otu_PA <- 1*((otu>0)==1)
head(otu_PA[,1:10])
metadata<-read_q2metadata("sample-metadata.txt")
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)
otu_rel <- apply(decostand(otu, method="total", MARGIN=2),1, mean) 
occ_abun <- tibble::rownames_to_column(as.data.frame(cbind(otu_occ, otu_rel)),'otu')
head(occ_abun[,1:3])
PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>% 
  gather(SampleID, abun, -otu) %>%
  left_join(metadata, by = 'SampleID') %>%
  group_by(otu, study) %>%
  summarise(genotype_freq=sum(abun>0)/length(abun),
            coreGenotype=ifelse(genotype_freq == 1, 1, 0)) %>% # 1 only if occupancy 1 with specific genotype, 0 if not
  group_by(otu) %>%
  summarise(sumF=sum(genotype_freq),
            sumG=sum(coreGenotype),
            nS=length(study)*2,
            Index=(sumF+sumG)/nS)
view(PresenceSum)
head(PresenceSum[,1:5])
otu_ranked <- occ_abun %>%
  left_join(PresenceSum, by='otu') %>%
  transmute(otu=otu,
            rank=Index) %>%
  arrange(desc(rank))
head(otu_ranked[,1:2]) #The end of ranking

###################################
#beta diversity calculation=Start(2)
###################################
#Setup otu table and tree
wUFaddition <- NULL
otu_start=otu_ranked$otu[1:2]                  #set out_start as the first row of otu_ranked
start_matrix <- as.matrix(otu[otu_start,])     #the row of otu_start in table(otu) (abundance in all samples)
tips_keep <- otu_start
otu_start_tree<-drop.tip(tree,tree$tip.label[-match(tips_keep, tree$tip.label)])
str(otu_start_tree)

#Manipulate the tree
# Create a matrix that maps each internal node to its 2 descendants
node.desc <- matrix(otu_start_tree$edge[order(otu_start_tree$edge[, 1]), ][, 2], byrow = TRUE, ncol = 2)
if( !all(rownames(start_matrix) == taxa_names(otu_start_tree)) ){
  start_matrix <- start_matrix[taxa_names(otu_start_tree), ]
}
ntip <- length(otu_start_tree$tip.label)
#create a matrix (OTUX(tip(abundance)+node(abundance))
edge_array <- matrix(0,
                     nrow = ntip + otu_start_tree$Nnode,
                     ncol = ncol(start_matrix),
                     dimnames = list(NULL, sample_names = colnames(start_matrix))
)
edge_array[1:ntip, ] <- start_matrix
ord.node <- order(ape::node.depth(otu_start_tree))[(ntip + 1):(ntip + otu_start_tree$Nnode)]
for (i in ord.node) {
  edge_array[i, ] <- colSums(edge_array[node.desc[i - ntip, ], , drop = FALSE], na.rm = TRUE)
}
edge_array <- edge_array[otu_start_tree$edge[, 2], ]


samplesums <- colSums(start_matrix)
spn <- combn(colnames(start_matrix), 2, simplify = FALSE)
distlist <- foreach::foreach(i = spn) %dopar% {
  A <- i[1]
  B <- i[2]
  AT <- samplesums[A]
  BT <- samplesums[B]
  wUF_branchweight <- abs(edge_array[, A] - edge_array[, B])/nReads
  numerator <- sum({otu_start_tree$edge.length * wUF_branchweight},na.rm = TRUE)
  return(numerator)
}
matIndices <- data.frame(do.call(rbind, spn)[, 2:1])
s_names <- paste(matIndices$X1, matIndices$X2, sep =' - ')
dis_start <-data.frame(cbind(s_names,unlist(distlist)))
names(dis_start)[2] <- 2 
wUFaddition <- rbind(wUFaddition,dis_start)
Process <- c(1:500)

##########################
#Loop(3-500)
##########################
for(i in 3:500){
  otu_add=otu_ranked$otu[i] 
  otu_start <- c(otu_start, otu_add)
  add_matrix <- as.matrix(otu[otu_add,])
  start_matrix <- rbind(start_matrix, add_matrix)
  tips_keep <- otu_start
  otu_loop_tree<-drop.tip(tree,tree$tip.label[-match(tips_keep, tree$tip.label)])
  node.desc <- matrix(otu_loop_tree$edge[order(otu_loop_tree$edge[, 1]), ][, 2], byrow = TRUE, ncol = 2)
  if( !all(rownames(start_matrix) == taxa_names(otu_loop_tree)) ){
    start_matrix <- start_matrix[taxa_names(otu_loop_tree), ]
  }
  ntip <- length(otu_loop_tree$tip.label)
  #create a matrix (OTUX(tip(abundance)+node(abundance))
  edge_array <- matrix(0,
                       nrow = ntip + otu_loop_tree$Nnode,
                       ncol = ncol(start_matrix),
                       dimnames = list(NULL, sample_names = colnames(start_matrix))
  )
  edge_array[1:ntip, ] <- start_matrix
  ord.node <- order(ape::node.depth(otu_loop_tree))[(ntip + 1):(ntip + otu_loop_tree$Nnode)]
  for (j in ord.node) {
    edge_array[j, ] <- colSums(edge_array[node.desc[j - ntip, ], , drop = FALSE], na.rm = TRUE)
  }
  edge_array <- edge_array[otu_loop_tree$edge[, 2], ]
  
  #Calculate weighted unifrac
  samplesums <- colSums(start_matrix)
  spn <- combn(colnames(start_matrix), 2, simplify = FALSE)
  distlist <- foreach::foreach(k = spn) %dopar% {
    A <- k[1]
    B <- k[2]
    AT <- samplesums[A]
    BT <- samplesums[B]
    wUF_branchweight <- abs(edge_array[, A] - edge_array[, B])/nReads
    numerator <- sum({otu_loop_tree$edge.length * wUF_branchweight},na.rm = TRUE)
    return(numerator)
  }
  matIndices <- data.frame(do.call(rbind, spn)[, 2:1])
  s_names <- paste(matIndices$X1, matIndices$X2, sep =' - ')
  dis_add <-data.frame(cbind(s_names,unlist(distlist)))
  names(dis_add)[2] <- i 
  wUFaddition <- left_join(wUFaddition, dis_add, by=c('s_names'))
  print(Process[i])
}
saveRDS(wUFaddition, "wUFaddtion.rds")

##########################
#Full
##########################
tips_keep <- row.names(otu)
trial.tree<-drop.tip(tree,tree$tip.label[-match(tips_keep, tree$tip.label)])
str(trial.tree)
# Create a matrix that maps each internal node to its 2 descendants
node.desc <- matrix(trial.tree$edge[order(trial.tree$edge[, 1]), ][, 2], byrow = TRUE, ncol = 2)
OTU <- as.matrix(otu)
if( !all(rownames(OTU) == taxa_names(trial.tree)) ){
  OTU <- OTU[taxa_names(trial.tree), ]
}
ntip <- length(trial.tree$tip.label)
#create a matrix (OTUX(tip(abundance)+node(abundance))
edge_array <- matrix(0,
                     nrow = ntip + trial.tree$Nnode,
                     ncol = ncol(OTU),
                     dimnames = list(NULL, sample_names = colnames(OTU))
)
edge_array[1:ntip, ] <- OTU
ord.node <- order(ape::node.depth(trial.tree))[(ntip + 1):(ntip + trial.tree$Nnode)]
for (i in ord.node) {
  edge_array[i, ] <- colSums(edge_array[node.desc[i - ntip, ], , drop = FALSE], na.rm = TRUE)
}
edge_array <- edge_array[trial.tree$edge[, 2], ]
samplesums <- colSums(OTU)

spn <- combn(colnames(OTU), 2, simplify = FALSE)
distlist <- foreach::foreach(i = spn) %dopar% {
  A <- i[1]
  B <- i[2]
  AT <- samplesums[A]
  BT <- samplesums[B]
  wUF_branchweight <- abs(edge_array[, A] - edge_array[, B])/nReads
  numerator <- sum({trial.tree$edge.length * wUF_branchweight},na.rm = TRUE)
  return(numerator)
}
matIndices <- data.frame(do.call(rbind, spn)[, 2:1])
s_names <- paste(matIndices$X1, matIndices$X2, sep =' - ')
dis_full <-data.frame(cbind(s_names,unlist(distlist)))
names(dis_full)[2] <- length(rownames(otu))
dis_full <-readRDS("dis_full.rds")
dis_full[,2] <-as.numeric(dis_full[,2])

rownames(dis_full) <- dis_full$x_names
temp_wUF <- dis_full
temp_wUF$x_names <- NULL
temp_wUF_matrix <- as.matrix(temp_wUF)

wUF_ranked <- data.frame(rank = as.factor(row.names(t(temp_wUF_matrix))),t(temp_wUF_matrix)) %>% 
  gather(comparison, wUF, -rank) %>%
  group_by(rank) %>%
  summarise(MeanwUF=mean(wUF)) %>%            # mean Bray-Curtis dissimilarity
  arrange(desc(-MeanwUF)) %>%
  mutate(proportionwUF=MeanwUF/max(MeanwUF))   # proportion of the dissimilarity explained by the n number of ranked OTUs
Increase=wUF_ranked$MeanwUF[-1]/wUF_ranked$MeanwUF[-length(wUF_ranked$MeanwUF)]
increasewUF <- data.frame(IncreasewUF=c(0,(Increase)), rank=factor(c(1:(length(Increase)+1))))
wUF_ranked <- left_join(wUF_ranked, increaseWUF)
wUF_ranked <- wUF_ranked[-nrow(wUF_ranked),]


# Final increase in wUF similarity of equal or greater then 2% 
wUF_Increase<- subset(wUF_ranked,IncreasewUF>=1.02)
LastCall_wUF <- max(as.numeric(as.character(wUF_Increase$rank)))
wUF_ranked$color <- ifelse(as.numeric(as.character(wUF_ranked$rank))<=LastCall_wUF, 'palegreen4')

#Creating plot of wUF distance
ggplot(wUF_ranked[1:150,], aes(x=factor(wUF_ranked$rank[1:150], levels=wUF_ranked$rank[1:150]))) +
  geom_point(aes(y=proportionwUF), color=wUF_ranked$color[1:150]) +
  theme_classic() + theme(strip.background = element_blank(),axis.text.x = element_blank(), axis.title.y = element_text(face="bold"), axis.title.x = element_text(face="bold")) +
  theme(axis.ticks.x = element_blank(),axis.text.y =element_text(face="bold")) +
  theme(axis.line = element_line(colour = "black", size = 1.0)) +
  geom_vline(xintercept=elbow_wUF, lty=3, col='#ffae19', cex=1.2) +
  geom_vline(xintercept=LastCall_wUF, lty=3, col='palegreen4', cex=1.2) +
  labs(x='ranked OTUs',y='weighted UniFrac similarity') +
  annotate(geom="text", x=elbow_wUF+40, y=.1, label=paste("elbow_wUF method"," (",elbow_wUF,")", sep=''), color="#ffae19",cex=4.5, face="bold")+    
  annotate(geom="text", x=LastCall_wUF+50, y=.48, label=paste("Last 2% increase (",LastCall_wUF,")",sep=''), color="palegreen4",cex=4.5, face="bold")
