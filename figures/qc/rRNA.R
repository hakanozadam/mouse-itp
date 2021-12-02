
library(reshape2)
library(tidyr)
library(dplyr)
library(ggplot2)

###############################################################
# We look at the ribosomal RNA abundance in the samples
###############################################################


monosome_stats  = "../../stats/monosome.csv"
hundred_1_stats = "../../stats/20201104-ITP-100-5mM-50-1.csv"
hundred_2_stats = "../../stats/20201209-ITP-100-5mM-6.csv"
hundred_3_stats = "../../stats/20210131-ITP-100-5mM-50_diluted-1.csv"

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

mm_stat_20210301_file = "../../stats/mm_20210301.csv"
mm_stat_20210318_file = "../../stats/mm_20210318.csv"
mm_stat_20210513_file = "../../stats/mm_20210513.csv"
mm_stat_20210614_file = "../../stats/mm_20210614.csv"


###############################################################

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

AXIS_THICKNESS = 0.35

################################################################################
#########                 F O N T   S I Z E S                          #########

FONT_LABEL_SIZE = 6
FONT_TITLE_SIZE = 8

PDF_resolution = 600
FIGURE_FONT = "helvetica"


stats_column_picker_vector <- c(1,2,3,4,5,6,7,8,9,16,17,18,19)
get_stats_summary_table <- function(csv_file){
  mytable <- read.csv(csv_file, header = TRUE, check.names = FALSE)
  return( mytable[stats_column_picker_vector,] )
}

monosome_raw_stats = get_stats_summary_table(monosome_stats)
hundred_1_raw_stats = get_stats_summary_table(hundred_1_stats)
hundred_2_raw_stats = get_stats_summary_table(hundred_2_stats)
hundred_3_raw_stats = get_stats_summary_table(hundred_3_stats)

monosome_rrna_percentages = monosome_raw_stats[5,-1]
colnames(monosome_rrna_percentages) = c("Monosome-1", "Monosome-2", "Monosome-3")
hundred_rrna_percentages = c( hundred_1_raw_stats[5,-1], hundred_2_raw_stats[5,-1], hundred_3_raw_stats[5,-1] )

percentage_df_pre = data.frame( monosome = as.vector(unlist(monosome_rrna_percentages)), hundred = hundred_rrna_percentages )

percentage_df = 
  melt(percentage_df_pre, value.name = "percentage", variable.name = "experiment") %>%
  group_by(experiment) %>%
  summarise(average = mean(percentage), se = sqrt(var(percentage) / 3))


rrna_percentage_plot = 
  ggplot(percentage_df, 
         aes(x=experiment, y=average)) +
  geom_bar(stat='identity', fill = BURNT_ORANGE) +
  geom_errorbar( aes(ymin=average-se, ymax= average + se), width = 0.4) + 
  theme_bw() +
  theme(plot.title   = element_text(hjust = 0.5, family = FIGURE_FONT, face = "plain", size = FONT_TITLE_SIZE),
        panel.border = element_blank(),
        panel.grid   = element_blank(),
        axis.text.y       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.text.x       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.title.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.title.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        legend.title      = element_blank(),
        axis.line.y       = element_line(colour = "black", size = AXIS_THICKNESS))  + 
  labs(title = "Ribosomal RNA", 
       x     = "Experiment", 
       y     = "Average %") + 
  scale_y_continuous(limits = c(0, 100), expand = c(0,0)) + 
  scale_x_discrete(labels = c("10M", "100"))


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

save_plot_pdf("rrna_percentage.pdf", rrna_percentage_plot, width=1.2, height=2)

# LEGEND
# Percentage of Ribosomal RNAs.
# Prior to mapping reads to the transcriptome, we first aligned them against ribosomal RNAS.
# The average of the percentages, from 3 replicates, in monosome perified bulk and 100 cell Ribo-ITP
# experiments are shown with standard error.

## Mouse RNA Mapped Reads

monosome_clipped_reads       = monosome_raw_stats[2,-1]
monosome_unique_mapped_reads = monosome_raw_stats[8,-1]
monosome_mapped_percentage   = (monosome_unique_mapped_reads / monosome_clipped_reads) * 100

hundred_rrna_clipped_reads     = c( hundred_1_raw_stats[2,-1], hundred_2_raw_stats[2,-1], hundred_3_raw_stats[2,-1] )
hundred_rrna_mapped_reads      = c( hundred_1_raw_stats[8,-1], hundred_2_raw_stats[8,-1], hundred_3_raw_stats[8,-1] )
hundred_rrna_mapped_percentage = (hundred_rrna_mapped_reads / hundred_rrna_clipped_reads) * 100

mapped_percentage_df_pre = data.frame( monosome = as.vector(unlist(monosome_mapped_percentage) ), 
                                       hundred = as.vector(unlist(hundred_rrna_mapped_percentage)) )

mapped_percentage_df = 
  melt(mapped_percentage_df_pre, value.name = "percentage", variable.name = "experiment") %>%
  group_by(experiment) %>%
  summarise(average = mean(percentage), se = sqrt(var(percentage) / 3))


mapped_percentage_plot = 
  ggplot(mapped_percentage_df, 
         aes(x=experiment, y=average)) +
  geom_bar(stat='identity', fill = BURNT_ORANGE) +
  geom_errorbar( aes(ymin=average-se, ymax= average + se), width = 0.4) + 
  theme_bw() +
  theme(plot.title   = element_text(hjust = 0.5, family = FIGURE_FONT, face = "plain", size = FONT_TITLE_SIZE),
        panel.border = element_blank(),
        panel.grid   = element_blank(),
        axis.text.y       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.text.x       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.title.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.title.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        legend.title      = element_blank(),
        axis.line.y       = element_line(colour = "black", size = AXIS_THICKNESS))  + 
  labs(title = "Mapped Reads", 
       x     = "Experiment", 
       y     = "Average %") + 
  scale_y_continuous(limits = c(0, 8), expand = c(0,0)) + 
  scale_x_discrete(labels = c("10M", "100"))

save_plot_pdf("mapped_read_percentage.pdf", mapped_percentage_plot, width=1.2, height=2)

# LEGEND
# Reads mapping to the transcriptome. Average percentage of the reads with error bars is shown.
# For averaging, we considered the percentage of reads, uniquely mapping to the transcriptome,
# in the set of reads having 3' adapters.

################################################################################

### Mouse rRNA Percentages


mm_stat_20210301 = get_stats_summary_table(mm_stat_20210301_file)
mm_stat_20210318 = get_stats_summary_table(mm_stat_20210318_file)
mm_stat_20210513 = get_stats_summary_table(mm_stat_20210513_file)
mm_stat_20210614 = get_stats_summary_table(mm_stat_20210614_file)

# Use inner_join so that the order of the rows is kept
mm_stats_all_raw = inner_join(
  inner_join(mm_stat_20210301, mm_stat_20210318),
  inner_join(mm_stat_20210513, mm_stat_20210614) )

###
experiment_list = 
  c("20210301-ITP-MII-25-B",
    "20210301-ITP-MII-50-A",
    "20210301-ITP-MII-50-B",
    "20210318-ITP-MII-50-B",
    "20210513-ITP-1cell-cross-50-A",
    "20210513-ITP-1cell-cross-50-B",
    "20210513-ITP-1cell-cross-50-C",
    "20210513-ITP-1cell-cross-50-D",
    "20210513-ITP-1cell-cross-50-E",
    "20210513-ITP-2cell-cross-50-B",
    "20210513-ITP-2cell-cross-50-C",
    "20210513-ITP-2cell-cross-50-F",
    "20210513-ITP-4cell-cross-50-B",
    "20210513-ITP-4cell-cross-50-C",
    "20210513-ITP-4cell-cross-50-D",
    "20210513-ITP-8cell-cross-50-A",
    "20210513-ITP-8cell-cross-50-B",
    "20210513-ITP-8cell-cross-50-C",
    "20210513-ITP-8cell-cross-50-D",
    "20210614-ITP-GV-50-A",
    "20210614-ITP-GV-50-B",
    "20210614-ITP-GV-50-C",
    "20210614-ITP-GV-50-E",
    "20210614-ITP-GV-50-F",
    "20210614-ITP-MII-50-D")

colnames(mm_stats_all_raw) = c("type", experiment_list)

mm_stats_all = mm_stats_all_raw[ c(colnames(mm_stats_all_raw)[1], experiment_list)  ]

mm_rrna_percentages = as.vector( unlist(mm_stats_all[5, -1] ) )

mm_percentage_df = data.frame(percentage = mm_rrna_percentages)

mouse_rrna_percentage_plot = 
ggplot(mm_percentage_df, aes(x=percentage)) + 
  geom_histogram(color = "white", fill = BURNT_ORANGE) + 
  theme_bw() +
  theme(plot.title   = element_text(hjust = 0.5, family = FIGURE_FONT, face = "plain", size = FONT_TITLE_SIZE),
        panel.border = element_blank(),
        panel.grid   = element_blank(),
        axis.text.y       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.text.x       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.title.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.title.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        legend.title      = element_blank(),
        axis.line.y       = element_line(colour = "black", size = AXIS_THICKNESS),
        axis.line.x       = element_line(colour = "black", size = AXIS_THICKNESS))  + 
  labs(title = "Mouse rRNA Content", 
       x     = "rRNA %", 
       y     = "# Samples") + 
  scale_y_continuous(labels = 0:4, limits = c(0, 4), expand = c(0,0))

save_plot_pdf("mouse_rrna_percentage.pdf", mouse_rrna_percentage_plot, width=2, height=2)

################################################################################

mouse_clipped_reads     = mm_stats_all[2,-1]
mouse_mapped_reads      = mm_stats_all[8,-1]
mouse_mapped_percentage = as.vector(unlist((mouse_mapped_reads / mouse_clipped_reads) * 100 ))

mouse_mapped_percentage_df = data.frame(percentage = mouse_mapped_percentage)

mouse_mapped_percentage_plot = 
  ggplot(mouse_mapped_percentage_df, aes(x=percentage)) + 
  geom_histogram(color = "white", fill = BURNT_ORANGE) + 
  theme_bw() +
  theme(plot.title   = element_text(hjust = 0.5, family = FIGURE_FONT, face = "plain", size = FONT_TITLE_SIZE),
        panel.border = element_blank(),
        panel.grid   = element_blank(),
        axis.text.y       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.text.x       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.title.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.title.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        legend.title      = element_blank(),
        axis.line.y       = element_line(colour = "black", size = AXIS_THICKNESS),
        axis.line.x       = element_line(colour = "black", size = AXIS_THICKNESS))  + 
  labs(title = "Mouse trans. mapping reads", 
       x     = "Mapped reads %", 
       y     = "# Samples") +
  scale_y_continuous( expand = c(0,0))

save_plot_pdf("mouse_mapped_percentage.pdf", mouse_mapped_percentage_plot, width=2, height=2)