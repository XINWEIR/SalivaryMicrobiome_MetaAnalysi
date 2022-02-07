library(readr)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)

setwd("~/Documents/Official/Figure6/Figure 6B")
err_rf_2206_1120 <- read_excel("err_rf_2206_1120.xlsx")

err_rf <- err_rf_2206_1120[, -1]
err_rf <- err_rf[, -29]
err_rf <- err_rf*100
Category = c(rep("Phyla",4), rep( "Classes", 4), rep("Orders",4), rep("Families",4), rep("Genera",4), rep("Species",4), rep("OTU",4))
rownames(err_rf) <- c("Hypervariable Region", "Geographic Location","Tobacco Usage", "Sample Type","Alcohol Use","Study","Gender","Age Range")


#box plots
rf_df <- data.frame(Category = c("Hypervariable Region", "Geographic Location", "Tobacco Usage", "Sample Type", "Alcohol Use", "Study", "Gender", "Age Range"))
rf_col<- list(Category = c('Age Range'="#E4DACE", 'Alcohol Use' = "#E5BB4B", 'Geographic Location' ="#498EAF", 'Sample Type' = "#631F16", Gender = "#499360", 'Hypervariable Region' = "#E29E93", 'Tobacco Usage' = "#D66C44", Study = "#A9B7C0"))
rf_bx_col <- c("#E29E93", "#498EAF", "#D66C44", "#631F16", "#E5BB4B", "#A9B7C0","#499360", "#E4DACE")
err_rf_row = data.frame(t(err_rf))
rf_ha_column = rowAnnotation(
  err_rate = anno_boxplot(err_rf_row, height = unit(4, "cm"), gp = gpar(fill = rf_bx_col), outline = FALSE),
  df = rf_df,
  col = rf_col,
  width = unit(6, "cm")
)
#Heatmap
library(viridis)
rf1 <-
  Heatmap(err_rf, name = "Err rate(%)",
        col = viridis(50),
        column_labels = c(rep(c("CLR", "RAR","RRA", "TSS"), 7)),
        column_names_gp = gpar(fontsize = 8),
        cluster_columns = FALSE,
        show_row_dend = FALSE,
        row_order = rownames(err_rf), 
        row_names_gp = gpar(fontsize = 8),
        column_order = order(as.numeric(gsub("column", "", colnames(err_rf)))),
        column_split = factor(Category, levels =c("Phyla", "Classes", "Orders", "Families", "Genera", "Species", "OTU")),
        right_annotation = rf_ha_column,
#       cluster_column_slices = FALSE,
        rect_gp = gpar(col = "black"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", err_rf[i, j]), x, y, gp = gpar(fontsize = 8.0, fontface = "bold", col = "white"))})

ggsave(rf1, file = "rf-2206.pdf", width = 300, height = 100, unit = 'mm')

