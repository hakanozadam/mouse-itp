################################################################################
########                    L I B R A R I E S                          #########

library(ribor)

library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
library(cowplot)
library(RColorBrewer)

library(Cairo)

setwd("~/projects/ribo-itp/repos/mouse-itp/figures/qc/")
################################################################################

################################################################################

MOUSE_MIN_LENGTH = 29
MOUSE_MAX_LENGTH = 35

################################################################################

################################################################################
#########                  C O L O R I N G                             #########

BURNT_ORANGE = "#bf5700"
UT_BLUE      = "#005f86"

PERCENTAGE_BARPLOT_PATERNAL_COLOR =  "#f0756e"
PERCENTAGE_BARPLOT_MATERNAL_COLOR = "#7ccba2"
PERCENTAGE_BARPLOT_OTHER_COLOR    = "#045175"

#PERCENTAGE_BARPLOT_PATERNAL_COLOR =  "#999999"
#PERCENTAGE_BARPLOT_MATERNAL_COLOR = "#E69F00"
#PERCENTAGE_BARPLOT_OTHER_COLOR    = "#56B4E9"

PERCENTAGE_BARPLOT_COLORS = c(PERCENTAGE_BARPLOT_OTHER_COLOR,
                              PERCENTAGE_BARPLOT_MATERNAL_COLOR,
                              PERCENTAGE_BARPLOT_PATERNAL_COLOR
)




#MAIN_PERCENTAGE_RNASEQ_FILL_COLOR     = "skyblue"
#MAIN_PERCENTAGE_RIBOSEQ_FILL_COLOR    = "#abe39d"
#MAIN_PERCENTAGE_ERRORBAR_COLOR = "#413bb8"

FOURCELL_BACKGROUND_COLOR = "#d9d9d9"

COUNT_NORMALIZATION_FACTOR = 10000

#RATIO_RIBO_COLOUR = "orange"
#RATIO_RNA_COLOUR  = "blue"

RATIO_RIBO_COLOUR = BURNT_ORANGE
RATIO_RNA_COLOUR  = "navy"

PERCENTAGE_DASHED_COLOR = "#088f99"

#CDS_GREEN  = "#00BA38"
#UTR5_BLUE  = "#619CFF"
#UTR3_GREEN = "619CFF"

CDS_GREEN  = rgb(124, 203, 162, maxColorValue = 255)
UTR5_BLUE  = rgb(4,   82,  117, maxColorValue = 255)
UTR3_GREEN = rgb(240, 116, 110,maxColorValue = 255)

################################################################################
#########                 F O N T   S I Z E S                          #########

FONT_LABEL_SIZE = 8
FONT_TITLE_SIZE = 9

PDF_resolution = 600
FIGURE_FONT = "sans"

################################################################################

# These values come from RiboFlow statistics
# Percentages of reads filtered out
ribo_itp_ribosomal_rna_percentages = c(60.45, 51.71, 53.48)
conventional_ribosomal_rna_percentages = c(90.07, 82.84, 86.66)

df_riboitp = data.frame( percentage = ribo_itp_ribosomal_rna_percentages, 
                         type = "Ribo-ITP" , 
                         s_e = sd(ribo_itp_ribosomal_rna_percentages),
                         average_p = mean(ribo_itp_ribosomal_rna_percentages))

df_conventional = data.frame( percentage = conventional_ribosomal_rna_percentages,
                              type = "Conventional",
                              s_e = sd(conventional_ribosomal_rna_percentages),
                              average_p = mean(conventional_ribosomal_rna_percentages) )

df_percentages = rbind(df_riboitp, df_conventional)

df_percentages$type = factor(df_percentages$type, levels = c("Ribo-ITP", "Conventional"))

  
  this_plot = 
    ggplot(data=df_percentages, aes(x= type, y=average_p, fill = type )  )  + 
    geom_bar(position = "dodge", stat="identity", alpha = 1 ) + 
    geom_errorbar( aes(  x        = type, 
                         ymin     = average_p - s_e, 
                         ymax     = average_p + s_e), 
                   width    = 0.4, 
                   alpha    = 0.4, 
                   size     = 0.4,
                   position = position_dodge(width = 0.9)) +
    geom_point( aes(  y = percentage) , 
                stat = "identity", 
                alpha = 0.6,
                position = position_dodge(width = 0.9),
                shape = 19,
                size = 1) + 
    theme(plot.title       = element_text(hjust = 0.5, family = FIGURE_FONT, face = "plain", size = FONT_TITLE_SIZE),
          panel.border     = element_blank(),
          panel.grid       = element_blank(),
          panel.background = element_blank(),
          axis.text.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
          axis.title.y     = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
          axis.text.x      = element_blank(),
          axis.title.x     = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
          legend.title     = element_blank(),
          axis.ticks.x     = element_blank(),    
          legend.text      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
          axis.line        = element_line(colour = "black", size = 0.35)
          #legend.key.size  = unit(0.12, 'inches')
          ) + 
    labs(title= "", x="Method", y="Ribosomal RNA %") + 
    scale_y_continuous(limits = c(0, 100), expand = c(0,0)) + 
    scale_fill_manual("legend", 
                      values = c( "Ribo-ITP" = CDS_GREEN, "Conventional" = UTR3_GREEN), 
                      breaks = c(0, 25,50,75,100)) 
  


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


save_plot_pdf("ribosomal_rna_comparison.pdf", 
              this_plot, 
              width = 1, height = 2)
