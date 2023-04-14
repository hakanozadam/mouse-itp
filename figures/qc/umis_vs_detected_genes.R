library(ribor)
library(Seurat)
library(reshape2)
library(smatr)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(cowplot)
library(dplyr)

################################################################################

# If you encounter problems with this script, run correlation.R or correlation_v2.R
# first. Then run this script.

setwd("~/projects/ribo-itp/repos/mouse-itp/figures/qc/")
source('./Ribo_Summary_Function.R')
source('./rename_experiments.R')

human_ribo_file = '../../../../data/human-itp_v4.ribo'
mouse_ribo_file = '../../../../data/mouse-itp_v5.ribo'

mouse_rnaseq_count_file = '../../raw_mouse_rnaseq_cds_counts.csv.gz'

external_file = "~/projects/ribo-itp/repos/mouse-itp/external_data/GSE162060_EE_RPFcounts.csv.gz"

human = Ribo(human_ribo_file, rename_default)
mm    = Ribo(mouse_ribo_file, rename = rename_default)



################################################################################


BURNT_ORANGE = "#bf5700"
UT_BLUE      = "#005f86"

MOUSE_MIN_LENGTH = 29
MOUSE_MAX_LENGTH = 35

FONT_LABEL_SIZE = 8
FONT_TITLE_SIZE = 9

PDF_resolution = 600
FIGURE_FONT    = "sans"


ribo_orange = rgb(228,88,10 , maxColorValue = 255)
rna_blue    = rgb(55,135,192, maxColorValue = 255)

AXIS_THICKNESS  = 0.35

basic_text_element = element_text(family = FIGURE_FONT, 
                                  face   = "plain", 
                                  size   = FONT_LABEL_SIZE)

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

################################################################################

min_len = MOUSE_MIN_LENGTH
max_len = MOUSE_MAX_LENGTH

mouse_exp_name_mapper = make_name_mapper( mm@experiments  )


experiment_name_mapper = mouse_exp_name_mapper

rename_experiments = function(experiment_name){
  return(experiment_name_mapper[[experiment_name]])
}

rename_experiments = Vectorize(rename_experiments, USE.NAMES = FALSE)


experiments_to_use = experiments(mm)


ribo_rc <- get_region_counts(mm,
                             experiment = experiments_to_use,
                             range.lower = min_len,
                             range.upper = max_len,
                             length      = TRUE,
                             transcript  = FALSE,
                             tidy = F,
                             alias       = TRUE,
                             region      = c("CDS"), 
                             compact = F)

ribo_rc$experiment = rename_experiments(ribo_rc$experiment)

rcw = dcast(ribo_rc, transcript ~ experiment) 

# GV:    3, 4, 5
# MII:   2, 3, 4
# 1cell: 1, 2, 5
selected_samples = c(2,3,6, 23, 24, 25, 21, 20, 19)
selected_rcw     = rcw[, selected_samples]

set.seed(3)
sub_40k = subsampleMatrix( selected_rcw, 40000)
sub_30k = subsampleMatrix( selected_rcw, 30000)
sub_20k = subsampleMatrix( selected_rcw, 20000)
sub_10k = subsampleMatrix( selected_rcw, 10000)
sub_5k  = subsampleMatrix( selected_rcw, 5000)
sub_2_5k  = subsampleMatrix( selected_rcw, 2500)
sub_1_250k  = subsampleMatrix( selected_rcw, 1250)
sub_0   = subsampleMatrix( selected_rcw, 0)

# # Plot number of genes detected as a function of UMI count
# # We should plot a line graph that shows cell to cell relationship as well the sum across cells. 
# old_genes_detected_df = data.frame( 
#   sample_depth = c( colSums(selected_rcw ), colSums(sub_40k), colSums(sub_30k), 
#                     colSums(sub_20k), colSums(sub_10k), colSums(sub_5k), colSums(sub_0)), 
#   genes_detected = c(   apply (selected_rcw, 2, function(x){sum(x > 0)}), 
#                         apply (sub_40k, 2, function(x){sum(x > 0)}), 
#                         apply (sub_30k, 2, function(x){sum(x > 0)}), 
#                         apply (sub_20k, 2, function(x){sum(x > 0)}), 
#                         apply (sub_10k, 2, function(x){sum(x > 0)}),
#                         apply (sub_5k, 2, function(x){sum(x > 0)}), 
#                         apply (sub_0, 2, function(x){sum(x > 0)})
#   ), 
#   sample_id =  rep(colnames(selected_rcw), 7)
#   
# )
# 
# 
# 
# umis_vs_genes_detected_plot = 
#   ggplot(
#     data = genes_detected_df, 
#     aes(x = sample_depth, y = genes_detected, group=sample_id) ) +
#   geom_line(linetype = "dashed", aes(col = sample_id) )  +  ylim(c( 0, 7000) ) +
#   geom_point(size = 0.7, aes(color = sample_id) ) + 
#   xlab("Number of CDS Mapping UMIs") + 
#   ylab("Genes Detected") + 
#   theme_bw() + 
#   theme(legend.title      = element_blank(),
#         axis.text.y       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
#         axis.text.x       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
#         axis.title.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
#         axis.title.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
#         legend.text       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE)) + 
#   #scale_x_continuous( expand = c(0.05, 0.05) )
#   scale_x_continuous( expand = c(0.05, 0.05), minor_breaks = NULL, labels = scales::label_number_si()) +
#   scale_y_continuous( expand = c(0.05, 0.05), minor_breaks = NULL, labels = scales::label_number_si(), limits = c(0, 7000))
# 
# 
# save_plot_pdf("umis_vs_genes_detected.pdf", umis_vs_genes_detected_plot, 
#               width = unit(3.45, "in"), height = unit(2.7, "in") )
################################################################################
### Detection upto 20k.
# Why 20k? Because the deepest experiments of GSE162060 have depth around 20k.



external_counts = read.csv(external_file, row.names = 1 )

depths_df = sort(apply( external_counts, 2, sum) )

# We get the 3 experiments having median depth.


len_external_experiments = length( depths_df  )
experiments_with_max_depth = depths_df[c(len_external_experiments-2, len_external_experiments - 1, len_external_experiments)]

external_deepest_experiment_counts = external_counts[, names(experiments_with_max_depth) ] 
colnames(external_deepest_experiment_counts)

set.seed(3)
external_sub_20k = subsampleMatrix( external_deepest_experiment_counts, 20000)
external_sub_10k = subsampleMatrix( external_deepest_experiment_counts, 10000)
external_sub_5k  = subsampleMatrix( external_deepest_experiment_counts, 5000)
external_sub_2_5k  = subsampleMatrix( external_deepest_experiment_counts, 2500)
external_sub_1_250k  = subsampleMatrix( external_deepest_experiment_counts, 1250)
external_sub_0   = subsampleMatrix( external_deepest_experiment_counts, 0)


external_genes_detected_upto_20k_df = data.frame( 
  sample_depth = c( 
                    colSums(external_sub_20k), colSums(external_sub_10k), colSums(external_sub_5k), 
                    colSums(external_sub_2_5k), colSums(external_sub_1_250k), colSums(external_sub_0)   ), 
  genes_detected = c(   
                        apply (external_sub_20k, 2, function(x){sum(x > 0)}), 
                        apply (external_sub_10k, 2, function(x){sum(x > 0)}), 
                        apply (external_sub_5k, 2, function(x){sum(x > 0)}), 
                        apply (external_sub_2_5k, 2, function(x){sum(x > 0)}),
                        apply (external_sub_1_250k, 2, function(x){sum(x > 0)}), 
                        apply (external_sub_0, 2, function(x){sum(x > 0)})
  ), 
  sample_id =  rep(colnames(external_deepest_experiment_counts), 6)
  
)


our_genes_detected_upto_20k_df = data.frame( 
  sample_depth = c( 
                    colSums(sub_20k), colSums(sub_10k), colSums(sub_5k), 
                    colSums(sub_2_5k), colSums(sub_1_250k), colSums(sub_0)), 
  genes_detected = c( 
                        apply (sub_20k, 2, function(x){sum(x > 0)}), 
                        apply (sub_10k, 2, function(x){sum(x > 0)}), 
                        apply (sub_5k, 2, function(x){sum(x > 0)}), 
                        apply (sub_2_5k, 2, function(x){sum(x > 0)}),
                        apply (sub_1_250k, 2, function(x){sum(x > 0)}), 
                        apply (sub_0, 2, function(x){sum(x > 0)})
  ), 
  sample_id =  rep(colnames(selected_rcw), 6)
  
)

genes_detected_upto_20k_df = rbind(external_genes_detected_upto_20k_df, our_genes_detected_upto_20k_df)

umis_vs_genes_detected_upto_20k_plot = 
  ggplot(
    data = genes_detected_upto_20k_df, 
    aes(x = sample_depth, y = genes_detected, group=sample_id) ) +
  geom_line(linetype = "dashed", aes(col = sample_id) )  +  ylim(c( 0, 7000) ) +
  geom_point(size = 0.7, aes(color = sample_id) ) + 
  xlab("Number of CDS Mapping UMIs") + 
  ylab("Genes Detected") + 
  theme_bw() + 
  theme(legend.title      = element_blank(),
        axis.text.y       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.text.x       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.title.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.title.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        legend.text       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE)) + 
  #scale_x_continuous( expand = c(0.05, 0.05) )
  scale_x_continuous( expand = c(0.05, 0.05), minor_breaks = NULL, labels = scales::label_number_si()) +
  scale_y_continuous( expand = c(0.05, 0.05), minor_breaks = NULL, labels = scales::label_number_si(), limits = c(0, 7000))

save_plot_pdf("umis_vs_genes_detected_upto_20k.pdf", umis_vs_genes_detected_upto_20k_plot, 
              width = unit(4.45, "in"), height = unit(2.7, "in") )

################################################################################

our_genes_detected_df = data.frame( 
  sample_depth = c( colSums(selected_rcw ), colSums(sub_40k), colSums(sub_30k), 
                    colSums(sub_20k), colSums(sub_10k), colSums(sub_5k), colSums(sub_2_5k),colSums(sub_1_250k), colSums(sub_0)), 
  genes_detected = c(   apply (selected_rcw, 2, function(x){sum(x > 0)}), 
                        apply (sub_40k, 2, function(x){sum(x > 0)}), 
                        apply (sub_30k, 2, function(x){sum(x > 0)}),
                        apply (sub_20k, 2, function(x){sum(x > 0)}), 
                        apply (sub_10k, 2, function(x){sum(x > 0)}),
                        apply (sub_5k, 2, function(x){sum(x > 0)}), 
                        apply (sub_2_5k, 2, function(x){sum(x > 0)}),
                        apply (sub_1_250k, 2, function(x){sum(x > 0)}), 
                        apply (sub_0, 2, function(x){sum(x > 0)})
  ), 
  sample_id =  rep(colnames(selected_rcw), 9),
  group = rep("This study" , length(colnames(selected_rcw)) * 9)
  
)


external_genes_detected_df = data.frame( 
  sample_depth = c( colSums(external_deepest_experiment_counts),
    colSums(external_sub_20k), colSums(external_sub_10k), colSums(external_sub_5k), 
    colSums(external_sub_2_5k), colSums(external_sub_1_250k), colSums(external_sub_0)   ), 
  genes_detected = c(   
    apply (external_deepest_experiment_counts, 2, function(x){sum(x > 0)}), 
    apply (external_sub_20k, 2, function(x){sum(x > 0)}), 
    apply (external_sub_10k, 2, function(x){sum(x > 0)}), 
    apply (external_sub_5k, 2, function(x){sum(x > 0)}), 
    apply (external_sub_2_5k, 2, function(x){sum(x > 0)}),
    apply (external_sub_1_250k, 2, function(x){sum(x > 0)}), 
    apply (external_sub_0, 2, function(x){sum(x > 0)})
  ), 
  sample_id =  rep(colnames(external_deepest_experiment_counts), 7)
  
)

genes_detected_df = rbind(external_genes_detected_df, our_genes_detected_df)

umis_vs_genes_detected_plot = 
  ggplot(
    data = genes_detected_df, 
    aes(x = sample_depth, y = genes_detected, group=sample_id) ) +
  geom_line(linetype = "dashed", aes(col = sample_id) )  +  ylim(c( 0, 7000) ) +
  geom_point(size = 0.7, aes(color = sample_id) ) + 
  xlab("Number of CDS Mapping UMIs") + 
  ylab("Genes Detected") + 
  theme_bw() + 
  theme(legend.title      = element_blank(),
        axis.text.y       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.text.x       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.title.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.title.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        legend.text       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE)) + 
  #scale_x_continuous( expand = c(0.05, 0.05) )
  scale_x_continuous( expand = c(0.05, 0.05), minor_breaks = NULL, labels = scales::label_number_si()) +
  scale_y_continuous( expand = c(0.05, 0.05), minor_breaks = NULL, labels = scales::label_number_si(), limits = c(0, 7000)) 


umis_vs_genes_detected_plot

save_plot_pdf("umis_vs_genes_detected.pdf", umis_vs_genes_detected_plot, 
              width = unit(4.45, "in"), height = unit(2.7, "in") )


################################################################################

median_index = as.integer(length( depths_df ) / 2 )
experiments_with_median_depth = depths_df[c(median_index-1, median_index, median_index + 1)]

external_median_experiment_counts = external_counts[, names(experiments_with_median_depth) ] 

set.seed(3)
external_median_sub_5k  = subsampleMatrix( external_median_experiment_counts, 5000)
external_median_sub_2_5k  = subsampleMatrix( external_median_experiment_counts, 2500)
external_median_sub_1_250k  = subsampleMatrix( external_median_experiment_counts, 1250)
external_median_sub_0 = subsampleMatrix( external_median_experiment_counts, 0)

our_genes_detected_upto_5k_df = data.frame( 
  sample_depth = c(  colSums(sub_5k), colSums(sub_2_5k),colSums(sub_1_250k), colSums(sub_0)), 
  genes_detected = c(   
                        apply (sub_5k, 2, function(x){sum(x > 0)}), 
                        apply (sub_2_5k, 2, function(x){sum(x > 0)}),
                        apply (sub_1_250k, 2, function(x){sum(x > 0)}), 
                        apply (sub_0, 2, function(x){sum(x > 0)})
  ), 
  sample_id =  rep(colnames(selected_rcw), 4),
  group = rep("This study" , length(colnames(selected_rcw)) * 4)
)

external_genes_detected_median_upto_5k_df = data.frame( 
  sample_depth = c(  colSums(external_median_sub_5k), colSums(external_median_sub_2_5k),
                     colSums(external_median_sub_1_250k), colSums(external_median_sub_0)), 
  genes_detected = c(   
    apply (external_median_sub_5k, 2, function(x){sum(x > 0)}), 
    apply (external_median_sub_2_5k, 2, function(x){sum(x > 0)}),
    apply (external_median_sub_1_250k, 2, function(x){sum(x > 0)}), 
    apply (external_median_sub_0, 2, function(x){sum(x > 0)})
  ), 
  sample_id =  rep(colnames(external_median_experiment_counts), 4),
  group = rep("GSE162060" , length(colnames(external_median_experiment_counts)) * 4)
  
)

genes_detected_upto_5k_median_df = rbind(external_genes_detected_median_upto_5k_df, our_genes_detected_upto_5k_df)

genes_detected_upto_5k_median_df$sample_id = factor(genes_detected_upto_5k_median_df$sample_id, 
                                                    levels = c(colnames(external_median_experiment_counts), colnames(selected_rcw)) )

umis_vs_genes_detected_upto_5k_median_plot = 
  ggplot(
    data = genes_detected_upto_5k_median_df, 
    aes(x = sample_depth, y = genes_detected, group=sample_id) ) +
  geom_line( aes(col = sample_id, linetype = group) )  +  ylim(c( 0, 2500) ) +
  geom_point(size = 0.7, aes(color = sample_id) ) + 
  xlab("Number of CDS Mapping UMIs") + 
  ylab("Genes Detected") + 
  theme_bw() + 
  theme(legend.title      = element_blank(),
        axis.text.y       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.text.x       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.title.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.title.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        legend.text       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE)) + 
  #scale_x_continuous( expand = c(0.05, 0.05) )
  scale_x_continuous( expand = c(0.05, 0.05), minor_breaks = NULL, labels = scales::label_number_si()) +
  scale_y_continuous( expand = c(0.05, 0.05), minor_breaks = NULL, labels = scales::label_number_si(), limits = c(0, 2500)) +
  scale_linetype_manual(values=c("dotted", "twodash"))
  


umis_vs_genes_detected_upto_5k_median_plot

save_plot_pdf("umis_vs_genes_detected_upto_5k_median.pdf", umis_vs_genes_detected_upto_5k_median_plot, 
              width = unit(5.5, "in"), height = unit(4, "in") )

#######

external_genes_detected_median_df = data.frame( 
  sample_depth = c( colSums(external_median_experiment_counts) , colSums(external_median_sub_5k), colSums(external_median_sub_2_5k),
                     colSums(external_median_sub_1_250k), colSums(external_median_sub_0)), 
  genes_detected = c(   
    apply (external_median_experiment_counts, 2, function(x){sum(x > 0)}),
    apply (external_median_sub_5k, 2, function(x){sum(x > 0)}), 
    apply (external_median_sub_2_5k, 2, function(x){sum(x > 0)}),
    apply (external_median_sub_1_250k, 2, function(x){sum(x > 0)}), 
    apply (external_median_sub_0, 2, function(x){sum(x > 0)})
  ), 
  sample_id =  rep(colnames(external_median_experiment_counts), 5),
  group = rep("GSE162060" , length(colnames(external_median_experiment_counts)) * 5)
  
)


genes_detected_median_df = rbind(external_genes_detected_median_df, our_genes_detected_df)

genes_detected_median_df$sample_id = factor(genes_detected_median_df$sample_id, 
                                                    levels = c(colnames(external_median_experiment_counts), colnames(selected_rcw)) )

umis_vs_genes_detected_median_plot = 
  ggplot(
    data = genes_detected_median_df, 
    aes(x = sample_depth, y = genes_detected, group=sample_id) ) +
  geom_line( aes(col = sample_id, linetype = group) )  +  ylim(c( 0, 7000) ) +
  geom_point(size = 0.7, aes(color = sample_id) ) + 
  xlab("Number of CDS Mapping UMIs") + 
  ylab("Genes Detected") + 
  theme_bw() + 
  theme(legend.title      = element_blank(),
        axis.text.y       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.text.x       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.title.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.title.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        legend.text       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE)) + 
  #scale_x_continuous( expand = c(0.05, 0.05) )
  scale_x_continuous( expand = c(0.05, 0.05), minor_breaks = NULL, labels = scales::label_number_si()) +
  scale_y_continuous( expand = c(0.05, 0.05), minor_breaks = NULL, labels = scales::label_number_si(), limits = c(0, 7000)) +
  scale_linetype_manual(values=c("dotted", "twodash"))

umis_vs_genes_detected_median_plot


save_plot_pdf("umis_vs_genes_detected_median.pdf", umis_vs_genes_detected_median_plot, 
              width = unit(5.5, "in"), height = unit(4, "in") )


inset_plot = umis_vs_genes_detected_upto_5k_median_plot + theme(legend.position = "none", axis.title.x = element_blank(),
                                                                axis.title.y = element_blank())

umis_vs_genes_detected_with_inset_median_plot = umis_vs_genes_detected_median_plot + annotation_custom( grob=ggplotGrob(inset_plot),
                                                        ymin = 0, ymax = 3400, xmin = 22000, xmax = 80000 )


save_plot_pdf("umis_vs_genes_detected_with_inset_median.pdf", umis_vs_genes_detected_with_inset_median_plot, 
              width = unit(6, "in"), height = unit(4.5, "in") )
