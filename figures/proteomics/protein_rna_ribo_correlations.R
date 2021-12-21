## Code to merge with Hakan
library(plotly)
library(ribor)
library(Seurat)
library(reshape2)
library(edgeR)
library(EnhancedVolcano)
library(pheatmap)
library(cowplot)

source('./Ribo_Summary_Function.R')
source("../qc/Ribo_Summary_Function.R")

mm = Ribo('../../../mouse-itp_v5.ribo', rename = rename_default)
mouse_rna_cds = read.csv('../../../cds_counts.csv.gz')
prots = read.csv ('../../../Protein_MouseEmbryo_PMID_29281840.csv')

min_len = 29
max_len = 35

ribo_rc <- get_region_counts(mm,
                             range.lower = min_len,
                             range.upper = max_len,
                             length      = TRUE,
                             transcript  = FALSE,
                             tidy = F,
                             alias       = TRUE,
                             region      = c("CDS"), 
                             compact = F)

rcw = dcast(ribo_rc, transcript ~ experiment)  

colnames(mouse_rna_cds)  = gsub(colnames(mouse_rna_cds), pattern = ".", fixed = T, replacement = "-")

for (id in 1:length(mouse_rna_cds$transcript)) { 
  mouse_rna_cds$transcript[id] = rename_default(mouse_rna_cds$transcript[id])
}
## Remove 20210607.RNAseq.4cell.cross.B
mouse_rna_cds = mouse_rna_cds[,-11]

# Mouse_RNA_CDS as before and doesn't include "X20210607-RNAseq-4cell-cross-B"
all_counts = merge(mouse_rna_cds, rcw, by = "transcript")

# Convert ribo and rna-seq into density 
cds_len = get_region_lengths(mm, alias = T)[,c(1,4)]
# Density corresponds to Reads per 1kb 

allcounts_density = all_counts
for (col_id in 2:48) { 
  allcounts_density[,col_id] = 1000* allcounts_density[,col_id] / 
    cds_len$CDS[match( allcounts_density$transcript, cds_len$transcript )]  
}

# CLR normalize density and average across replicates
# Note that we switched the margin
set.seed(3)
allcounts_cds_clr_density = NormalizeData(allcounts_density[,-1], 
                                          scale.factor = 10000, margin = 2, 
                                          normalization.method= "CLR",verbose = TRUE)
rownames(allcounts_cds_clr_density) = all_counts[,1]
colSums(allcounts_cds_clr_density)

stage = unlist(lapply ( strsplit(colnames(all_counts)[-1], split = "-" ) , "[[", 3 )  ) 
method = unlist(lapply ( strsplit(colnames(all_counts)[-1], split = "-" ) , "[[", 2 )  ) 
groups = as.factor(paste(method, stage,  sep= "_") ) 

mean_clr_density = t(apply(allcounts_cds_clr_density, 1, FUN = function (x){ tapply (as.numeric(x), groups, mean)}) ) 


## Add protein expression
prots_rna_ribo = mean_clr_density[strip_extension(rownames(mean_clr_density) ) %in% prots$Gene.Name, ] 
prots_rna_ribo = cbind( prots_rna_ribo, 
                        prots[match(strip_extension(row.names(prots_rna_ribo)) , prots$Gene.Name), ]  ) 

## Remove annotations/blastocyst from the proteomics data
prots_rna_ribo = prots_rna_ribo[, -c(19:22)]

# Note that 451 ids in proteomics data did not match ids in ribo/rna
dim(prots_rna_ribo)

# Subset to non-zero expression; 29 transcripts with no detectable rna or ribo 
prots_rna_ribo = prots_rna_ribo[apply(prots_rna_ribo[,1:12], 1, sum) > 0, ]
dim(prots_rna_ribo)


## We should be using Spearman's correction before reporting these. 
## Estimation of the reliability of the assays based on replicate correlation: 
rna_reliability = 0.79
ribo_reliability = 0.71
prot_reliability = 0.8


## Note that there are some negative correlation that will get replaced with zero in the following figure
## If using this version might either change color scheme or label as less than zero
## Will split this into RNA - Prot and Ribo-Prot 
rna_prot = prots_rna_ribo[ , c(13, 11, 12, 7:10, 14:18)]
ribo_prot = prots_rna_ribo[, c(13, 5, 6, 1:4,14:18 )]

rna_prot_cors = cor(rna_prot[,-1], method = "spearman") / sqrt(prot_reliability * rna_reliability)
ribo_prot_cors = cor(ribo_prot[,-1], method = "spearman") / sqrt(prot_reliability * ribo_reliability)

rna_prot_cors = rna_prot_cors[1:6, 7:11]
ribo_prot_cors = ribo_prot_cors[1:6, 7:11]

## RNA - Prot and Ribo - Prot Correlations for Heatmap
rna_prot_cors
ribo_prot_cors

################################################################################

labels_heatmap = c("GV", "MII", "1-cell", "2-cell", "4-cell", "8-cell", "Morula")

rna_prot_na_inserted_df = cbind( rep(NA, 6),rep(NA, 6), rna_prot_cors    )
rna_prot_na_inserted_df = rbind( rna_prot_na_inserted_df, rep(NA, 7 ))
rownames(rna_prot_na_inserted_df) = labels_heatmap
colnames(rna_prot_na_inserted_df) = labels_heatmap

## 4-8 cell color list
link_list$color = c("#D1B6B1", "#0A0A55", "#D1B6B1", "#090934", 
                    "#CEA7A0", "#3C3C9F",  "#D0AFA9", "#1E1E95")

heatmap_breaks    = c(-1,   0,    0.31, 0.35,  0.39, 
                      0.43, 0.47, 0.51, 0.55,  0.59)
heatmap_colors    = c("#D1B6B1" , "#D3D3D3", "#B4B4C8" ,"#9696BE", "#7878B4", 
                      "#5A5AA9" ,"#3C3C9F", "#2d2d9a" ,"#1E1E95", "#0f0f8f", 
                      "#00008B")

combined_values = sort( c( rna_prot_cors[rna_prot_cors   > 0.34],
                           ribo_prot_cors[ribo_prot_cors > 0.34]) )


rna_breaks_pre = seq(from = combined_values[1], to = combined_values[length(combined_values)], length.out = 7 )

rna_heatmap_colors = colorRampPalette(
  colors = c("#B4B4C8", "darkblue"))(length(rna_breaks_pre) )

rna_heatmap_colors = c("#D1B6B1", "#D3D3D3",  rna_heatmap_colors)

rna_breaks = c(-1, 0 ,0.32, rna_breaks_pre)

heatmap_rna_prot_cor = pheatmap (
  rna_prot_na_inserted_df, 
  color             =  rna_heatmap_colors,
  breaks            = rna_breaks, 
#  main              = "Protein & RNA Correlation",
  cluster_rows      = F, 
  cluster_cols      = F,
  legend            = F,
  na_col            = "white",
  display_numbers   = T,
  number_color      = "#FA8128",
  fontsize          = 7,
  fontsize_col      = 7,
  fontsize_row      = 7)
print(heatmap_rna_prot_cor)

################################################################################

## Ribo Heatmap

ribo_prot_cors

ribo_prot_na_iserted = cbind( rep(NA, 6), rep(NA, 6) , ribo_prot_cors )
ribo_prot_na_iserted = rbind(ribo_prot_na_iserted, rep(NA, 7))

rownames(ribo_prot_na_iserted) = labels_heatmap
colnames(ribo_prot_na_iserted) = labels_heatmap

# NOTE: We use RNA heatmap breaks and colors to make sure that 
# both heatmaps use the same color scheme.

#ribo_breaks_pre = sort(ribo_prot_na_iserted[ribo_prot_na_iserted > 0]  )

#ribo_heatmap_colors = colorRampPalette(
#  colors = c("lightgray", "darkblue"))(length(ribo_breaks_pre) )

#ribo_heatmap_colors = c("#D1B6B1",  ribo_heatmap_colors)

#ribo_breaks = c(-1, 0 , ribo_breaks_pre)

heatmap_ribo_prot_cor = pheatmap (
  ribo_prot_na_iserted, 
  color             = rna_heatmap_colors,
  breaks            = rna_breaks, 
  #main              = "Protein & Ribo Correlation",
  cluster_rows      = F, 
  cluster_cols      = F,
  legend            = F,
  width             = 0.3,
  height            = 0.3,
  na_col            = "white",
  display_numbers   = T,
  number_color      = "#FA8128",
  fontsize          = 7,
  fontsize_col      = 7,
  fontsize_row      = 7)
print(heatmap_ribo_prot_cor)

################################################################################
PDF_resolution = 600

get_output_file_path = function(file_name, output_folder = "pdf"){
  this_path = paste( output_folder, file_name, sep = "/"  )
  return(this_path)
}

save_plot_pdf = function(filename, this_plot, width = NA, height = NA){
  this_file = get_output_file_path(filename)
  print(this_file)
  ggsave(this_file, 
         plot   = this_plot, 
         device = cairo_pdf, 
         width  = width,
         height = height,
         dpi    = PDF_resolution )
  
}

rna_title = ggdraw() + 
  draw_label(
    "Protein & RNA correlation",
    fontface = 'bold',
    fontfamily = 'helvetica',
    hjust = 0.5
  )

save_plot_pdf("rna_protein_cor_heatmap.pdf", heatmap_rna_prot_cor, width = 2.4, height = 2.4)

save_plot_pdf("ribo_protein_cor_heatmap.pdf", heatmap_ribo_prot_cor, width = 2.4, height = 2.4)

combined_plot = plot_grid( heatmap_rna_prot_cor$gtable, heatmap_ribo_prot_cor$gtable, nrow = 1 )

save_plot_pdf("ribo_protein_rna_cor_heatmap.pdf", combined_plot, width = 2.4 * 2, height = 2.4)
