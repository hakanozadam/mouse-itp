## Code to merge with Hakan
library(plotly)
library(ribor)
library(Seurat)
library(reshape2)
library(edgeR)
library(EnhancedVolcano)
library(pheatmap)
library(stringr)

setwd("~/projects/mice_cell_development/mouse-itp/figures/proteomics")

correlation_file = "./TE_Prot_Cors_Corrected.csv"

this_cor_matrix = read.csv(correlation_file, row.names = 1)

row.names(this_cor_matrix) = str_replace(row.names(this_cor_matrix), c("_") , c(" "))
row.names(this_cor_matrix) = str_replace(row.names(this_cor_matrix), c("X") , c(""))
row.names(this_cor_matrix) = str_replace(row.names(this_cor_matrix), c("\\.") , c(""))

colnames(this_cor_matrix) = row.names(this_cor_matrix)

paletteLength = 100


input_heatmap = this_cor_matrix[1:6, 7:11]
input_heatmap

myBreaks <- c(seq(-0.5, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(input_heatmap)/paletteLength, max(input_heatmap), length.out=floor(paletteLength/2)))


myColor       = colorRampPalette(c("#E2E60A", "white", "#228B22"  ))(paletteLength)


te_and_prot_correlation_plot = 
pheatmap (
  input_heatmap,
  #color        =  colorRampPalette(brewer.pal( 9 ,"Greens"))(100),
  #color        =  colorRampPalette(brewer.pal( 9 ,"PRGn"))(100),
  color = myColor,
  breaks = myBreaks,
  #color       =  colorRampPalette( rev(grDevices::heat.colors(100) )  )(length(breaks_manual)),
  legend_breaks = seq(-1, 1.5 , 0.5),
  legend_labels = seq(-1, 1.5 , 0.5),
  cluster_rows = F,
  cluster_cols = F,
  display_numbers = T,
  main   = "Cor. of proteomics and TE",
  fontsize          = 7,
  fontsize_col      = 7,
  fontsize_row      = 7)

te_and_prot_correlation_plot


get_output_file_path = function(file_name, output_folder = "pdf"){
  this_path = paste( output_folder, file_name, sep = "/"  )
  return(this_path)
}

save_plot_pdf = function(filename, this_plot, width = NA, height = NA){
  this_file = get_output_file_path(filename)
  print(this_file)
  ggsave(this_file, 
         plot   = this_plot, 
         #device = cairo_pdf, 
         width  = width,
         height = height,
         dpi    = PDF_resolution )
}



save_plot_pdf("te_prot_cors_corrected.pdf", te_and_prot_correlation_plot, width = 2.6, height = 2.4)
