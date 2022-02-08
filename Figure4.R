library("ComplexHeatmap")
library("circlize")
library(viridis)

Core.MRA.2206 <- readRDS("Core_MRA_OTU.rds")
Core_MRA_log <- readRDS("Core_MRA_subgroups.rds")

#Taxonomy
setwd("~/Documents")
library(qiime2R)
taxonomy_1522 <- read_qza("ref-1522-taxonomy.qza")$data
taxonomy_1522_sep <- parse_taxonomy(taxonomy_1522)
HOMD_1522 <- add_rownames(taxonomy_1522_sep, var= "otu")
HOMD_taxonomy_core_new <- left_join(HOMD_1522, Core.MRA.2206)


#Column=Phylum
otu_phylum <- data.frame(Phylum = HOMD_taxonomy_core_new$Phylum)
rownames(otu_phylum) <- HOMD_taxonomy_core_new$otu
otu_phylum <- as.matrix(otu_phylum)
phylum.name <- "Saccharibacteria"
otu_phylum[24,1]<- phylum.name
  #Order of hist 1 & hist 2 == match
otu_phylum <- as.data.frame(otu_phylum[rownames(Core_MRA_log), ])
colnames(otu_phylum)[1]<-"Phylum"

#Column=MRA
otu_MRA_core <- HOMD_taxonomy_core_new[, c(1,8)]
otu_MRA_core$MRA <- round(as.numeric(otu_MRA_core$MRA)*100,2)
otu_MRA_core_matrix <- as.matrix(otu_MRA_core[,-1])
rownames(otu_MRA_core_matrix) <- otu_MRA_core$otu
otu_MRA_core <- as.data.frame(otu_MRA_core_matrix[rownames(Core_MRA_log), ])
  #Order of hist 1 & hist 3 == match
rownames(Core_MRA_log) == rownames(otu_MRA_core)
colnames(otu_MRA_core)[1]<-"MRA"

#Taxonomy annotation=l6&l7
tax_ano <- data.frame(otu=rownames(Core_MRA_log),Core=1)
Core_g_s <- apply(HOMD_taxonomy_core_new_new[,7:8],1, function(x)paste0(paste(names(x[!is.na(x)]),x[!is.na(x)], sep = "="),collapse = ";"))
Core_gs_df <- data.frame(otu =HOMD_taxonomy_core_new_new$otu, Core_g_s)
tax_ano <-left_join(tax_ano, Core_gs_df)
rownames(Core_MRA_log) == tax_ano$otu
Core_g_s <- tax_ano[,3]

#Top annotation=Figure4D
group_test <- readRDS("group_relabun.rds")
ha_df <- data.frame(Category = c(rep("Ethnicity", 3), rep("AlcoholUse", 2), rep("TobaccoUsage", 2), rep("Gender", 2), rep("SampleType", 3), rep("HypervariableRegion", 2), rep("AgeRange", 3)))
ha_col<- list(Category = c(AgeRange="#E4DACE", AlcoholUse = "#E5BB4B", Ethnicity ="#498EAF", SampleType = "#631F16", Gender = "#499360", HypervariableRegion = "#E29E93", TobaccoUsage = "#D66C44"))
bx_col <- c(rep("#498EAF",3),rep("#E5BB4B", 2), rep("#D66C44",2), rep("#499360",2), rep("#631F16",3), rep("#E29E93",2), rep("#E4DACE",3))
ha_column = HeatmapAnnotation(
  prop = anno_boxplot(group_test, height = unit(4, "cm"), gp = gpar(fill = bx_col), outline = FALSE),
  df = ha_df,
  col = ha_col,
  annotation_name_side = "left"
)

#heatmap1=relative abundance
Category = c(rep("Ethnicity", 3), rep("AlcoholUse", 2), rep("TobaccoUsage", 2), rep("Gender", 2), rep("SampleType", 3), rep("HypervariableRegion", 2), rep("AgeRange", 3))
ht1 = ComplexHeatmap::Heatmap(Core_MRA_log,
                        col = viridis(50),
                        name = "Log_MRA",
                        show_column_dend = FALSE,
                        row_dend_width = unit(4, "cm"),
                        column_order = colnames(Core_MRA),
                        column_names_rot = 45,
                        row_names_gp = gpar(fontsize = 5.6),
                        row_km = 4, show_parent_dend_line = FALSE,
                        top_annotation = ha_column, #column annotation
                        column_split = Category,
                        rect_gp = gpar(col = "black")
) 
ht1

#heatmap2=Phylum
col_phylum = c(Actinobacteria = "#324856", Bacteroidetes = "#4A746A", Firmicutes ="#D18237", Fusobacteria = "#CA4026", Proteobacteria = "#9D1E31", Saccharibacteria = "#F5E68B")
ht2 = Heatmap(otu_phylum, 
              name = "Phylum", 
              col = col_phylum, 
              row_names_gp = gpar(fontsize = 5.6),
              rect_gp = gpar(col = "black"))
ht2

#heatmap3=MRA
col_runif = colorRamp2(c(0, max(otu_MRA_core)), c("#E8C7B6", "#D5683D"))
ha_row = rowAnnotation(foo = anno_text(Core_g_s, location = 0 , just = "left", gp = gpar(fontsize = 6)))
ht3 = Heatmap(otu_MRA_core, 
              name = "MRA", 
              col = col_runif, 
              right_annotation = ha_row,
              row_names_gp = gpar(fontsize = 5.6),
              rect_gp = gpar(col = "black"),
              width = 2,
              cell_fun = function(j, i, x, y, width, height, fill) {
               grid.text(sprintf("%.1f", otu_MRA_core[i, j]), x, y, gp = gpar(fontsize = 5.6))}
               )
ht3

#Combine heatmaps
ht1 + ht2 + ht3

