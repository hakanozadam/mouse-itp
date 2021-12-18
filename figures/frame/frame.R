library(ribor)
library(reshape2)
library(edgeR)
library(smatr)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(cowplot)
library(dplyr)

###############################################################################

mouse_raw_frame_file      = "./mouse_raw_frame_percentages.csv"
mouse_adjusted_frame_file = "./mouse_adjusted_frame_percentages.csv"

human_raw_frame_file       = "./human_raw_frame_percentages.csv"
human_adjustted_frame_file = "./human_adjusted_frame_percentages.csv" 

###############################################################################

BURNT_ORANGE = "#bf5700"
UT_BLUE      = "#005f86"

PERCENTAGE_BARPLOT_PATERNAL_COLOR =  "#f0756e"
PERCENTAGE_BARPLOT_MATERNAL_COLOR = "#7ccba2"
PERCENTAGE_BARPLOT_OTHER_COLOR    = "#045175"

FOURCELL_BACKGROUND_COLOR = "#d9d9d9"

RATIO_RIBO_COLOUR = BURNT_ORANGE
RATIO_RNA_COLOUR  = "navy"

PERCENTAGE_DASHED_COLOR = "#088f99"

CDS_GREEN  = rgb(124, 203, 162, maxColorValue = 255)
UTR5_BLUE  = rgb(4,   82,  117, maxColorValue = 255)
UTR3_GREEN = rgb(240, 116, 110,maxColorValue = 255)

AXIS_THICKNESS = 0.35

FONT_LABEL_SIZE = 7
FONT_TITLE_SIZE = 8

PDF_resolution = 600
FIGURE_FONT = "helvetica"


################################################################################

# Mouse Data
mouse_raw_df           = read.csv(mouse_raw_frame_file, row.names = 1)
colnames(mouse_raw_df) = c(0, 1, 2)

mouse_adjusted_df           = read.csv(mouse_adjusted_frame_file, row.names = 1)
colnames(mouse_adjusted_df) = c(0, 1, 2)

# Huamn Data
human_raw_df           = read.csv(human_raw_frame_file, row.names = 1)
colnames(human_raw_df) = c(0, 1, 2)

human_adjusted_df           = read.csv(human_adjustted_frame_file, row.names = 1)
colnames(human_adjusted_df) = c(0, 1, 2)

################################################################################

mouse_single_cell_idx       = c(1:9, 20:25)
mouse_single_cell_raw       = mouse_raw_df[mouse_single_cell_idx,] 
mouse_single_cell_adjusted  = mouse_adjusted_df[mouse_single_cell_idx,] 

human_monosome_idx    = 1:3
human_100_cell_idx    = 4:6

human_monosome_adjusted   = human_adjusted_df[human_monosome_idx,]
human_100cell_adjusted = human_adjusted_df[human_100_cell_idx,]

human_monosome_raw   = human_raw_df[human_monosome_idx,]
human_100cell_raw = human_raw_df[human_100_cell_idx,]

################################################################################

compute_se = function(x){
  return( sd(x)/sqrt(length(x))   )
}

################################################################################

mouse_single_cell_raw_means      = apply( mouse_single_cell_raw, MARGIN=2, FUN = mean)
mouse_single_cell_adjusted_means = apply( mouse_single_cell_adjusted, MARGIN=2, FUN = mean)

mouse_single_cell_raw_se      = apply( mouse_single_cell_raw, MARGIN=2, FUN = compute_se)
mouse_single_cell_adjusted_se = apply( mouse_single_cell_adjusted, MARGIN=2, FUN = compute_se)

mouse_raw_df      = data.frame(means  = mouse_single_cell_raw_means, 
                               se     = mouse_single_cell_raw_se,
                               frames = 0:2)
mouse_adjusted_df = data.frame( means  = mouse_single_cell_adjusted_means, 
                                se     = mouse_single_cell_adjusted_se,
                                frames = 0:2)

################################################################################

human_100cell_adjusted_means  = apply( human_100cell_adjusted, MARGIN = 2, FUN = mean )
human_monosome_adjusted_means = apply( human_monosome_adjusted, MARGIN = 2, FUN = mean )

human_100cell_adjusted_se  = apply( human_100cell_adjusted, MARGIN = 2, FUN = compute_se )
human_monosome_adjusted_se = apply( human_monosome_adjusted, MARGIN = 2, FUN = compute_se )

human_100cell_adjusted_df = data.frame(means  = human_100cell_adjusted_means,
                                       se     = human_100cell_adjusted_se,
                                       frames = 0:2)

human_monosome_adjusted_df = data.frame(means  = human_monosome_adjusted_means,
                                        se     = human_monosome_adjusted_se,
                                        frames = 0:2)

################################################################################

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



frame_barplot = function(df, plot_title = "Frame Percentage", ymax = 50){
  
  this_plot = ggplot(df, 
                     aes(x = frames, y=means)) +
    geom_bar(stat="identity", fill = BURNT_ORANGE, width = 0.8) +
    geom_errorbar( aes(  x        = frames, 
                         ymin     = means - se, 
                         ymax     = means + se), 
                   width    = 0.3, 
                   alpha    = 0.9, 
                   size     = 0.5) +
    theme_bw() + 
    theme(plot.title   = element_text(hjust = 0.5, family = FIGURE_FONT, face = "plain", size = FONT_TITLE_SIZE),
          panel.border = element_blank(),
          panel.grid   = element_blank(),
          axis.text.y       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
          axis.text.x       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
          axis.title.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
          axis.title.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
          axis.line         = element_line(colour = "black", size = AXIS_THICKNESS),
          legend.title      = element_blank() ) + 
    labs(title = plot_title, 
         fill  = "Region", 
         x     = "Frame", 
         y     = "Percentage") + 
    scale_y_continuous( expand = c(0,0), limits = c(0, ymax))+  
    scale_x_continuous( breaks = 0:2, labels = c("0", "1", "2") )
  
  return(this_plot)
}

################################################################################

## Mouse Plots

mouse_raw_barplot = frame_barplot(mouse_raw_df, plot_title = "Mouse Single Cell Frames", ymax = 40)

save_plot_pdf("mouse_raw_barplot.pdf", mouse_raw_barplot, width = 1.9, height = 2.3)

mouse_adjusted_barplot = frame_barplot(mouse_adjusted_df, plot_title = "Mouse Single Cell Frames", ymax = 50)

save_plot_pdf("mouse_adjusted_barplot.pdf", mouse_adjusted_barplot, width = 1.9, height = 2.3)

################################################################################

# Human Plots

human_100cell_adjusted_plot = frame_barplot(human_100cell_adjusted_df, plot_title = "Human 100-Cell Frames", ymax = 50)

human_monosome_adjusted_plot = frame_barplot(human_monosome_adjusted_df, plot_title = "Human 10M-Cell Frames", ymax = 50)

save_plot_pdf("human_100cell_adjusted.pdf", human_100cell_adjusted_plot, width = 1.9, height = 2.3)

save_plot_pdf("human_monosome_adjusted.pdf", human_monosome_adjusted_plot, width = 1.9, height = 2.3)

################################################################################

# Notes:
#   Mouse: Adjusted mean percentage, for frame 0 is 46.99% with standard error 0.18
#   For the mouse data, we used all single cell experiments from GV, MII and 1-cell embryo stages.
#   
#   Human 100-cell:  Adjusted mean percentage, for frame 0 is 48.34% with standard error 0.09
#   Human Monosome: Adjusted percentage, for frame 0 is 47.58% with standard error 0.13
#   
# Frame Adjustment:
#   For a read, we define its frame to be the difference between module 3 of the 
#   difference between the 5' end of the read and the translation start site.
#   
#   For each RPF length, we considered the last two nucleotides at the 3' end of the read
#   and one nucleotide downstream of the 3' end.
#   For each of these triplets of these nucleotides, we counted the number of reads
#   mapping to the frames 0,1,2.
#   
#   Then we applied cyclic shifts to each triplet so that the maximum value is at frame 0.
#   After aggregating all triplets for all RPF lengths,
#   we obtain the total adjusted counts of the frames 0, 1, 2.
#   
  
