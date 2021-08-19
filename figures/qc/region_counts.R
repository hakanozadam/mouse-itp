#REGION COUNTS

mouse_ribo_file = "../../../mouse-itp_v5.ribo"
human_ribo_file = "../../../../itp/human-itp_v4.ribo"


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

#Heatmap related packages
library(pheatmap)

# We don't need the following two libraries.
#library(seriation)
#library(dendextend)

################################################################################

mouse_ribo = Ribo(mouse_ribo_file, rename = rename_default)
human_ribo = Ribo(human_ribo_file, rename = rename_default)

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
FONT_TITLE_SIZE = 10

PDF_resolution = 600
FIGURE_FONT = "helvetica"


################################################################################
#######                   E X P E R I M E N T S                      ###########


make_name_mapper = function(experiment_names){
  exp_names        = sort( experiment_names )
  unique_raw_names = unique(exp_names )
  group_names      = sort(unique(unlist( lapply(strsplit( unique_raw_names, split = "-" ) , "[[", 3 ) ) ) )
  
  original_exp_names     = c()
  new_experiment_names   = c()
  list_counter           = 1
  experiment_name_mapper = list()
  
  for(g in group_names){
    elements_of_group = unique_raw_names[grep(  g, unique_raw_names )]
    size_of_group     = length( elements_of_group  )
    
    for( i in seq(1:size_of_group) ){
      this_basename = strsplit( elements_of_group[i], split = "-"   )[[1]][3]
      new_name = paste(  this_basename, i , sep = "-" )
      print( c( elements_of_group[i], new_name )  )
      
      original_exp_names   = c(original_exp_names, elements_of_group[i])
      new_experiment_names = c(new_experiment_names, new_name) 
      
      experiment_name_mapper[[list_counter]] = new_name
      list_counter = list_counter + 1
    }
    
  }
  
  names(  experiment_name_mapper ) = original_exp_names
  
  return(experiment_name_mapper)
}




mouse_exp_name_mapper = make_name_mapper( mouse_ribo@experiments  )


experiment_name_mapper = mouse_exp_name_mapper

rename_experiments = function(experiment_name){
  return(experiment_name_mapper[[experiment_name]])
}

rename_experiments = Vectorize(rename_experiments, USE.NAMES = FALSE)


human_exp_name_mapper = make_name_mapper( human_ribo@experiments  )

human_rename_experiments = function(experiment_name){
  return(human_exp_name_mapper[[experiment_name]])
}

human_rename_experiments = Vectorize(human_rename_experiments, USE.NAMES = FALSE)
################################################################################
################################################################################
################################################################################
################################################################################


################################################################################
###########               R E G I O N     C O U N T S           ################

mouse_region_counts = get_region_counts(ribo.object = mouse_ribo,
                                        region      = c("UTR5", "CDS", "UTR3"),
                                        range.lower = MOUSE_MIN_LENGTH,
                                        range.upper = MOUSE_MAX_LENGTH,
                                        compact     = FALSE)

human_region_counts = get_region_counts(ribo.object = human_ribo,
                                        region      = c("UTR5", "CDS", "UTR3"),
                                        range.lower = MOUSE_MIN_LENGTH,
                                        range.upper = MOUSE_MAX_LENGTH,
                                        compact     = FALSE)

mouse_region_counts$experiment = rename_experiments(mouse_region_counts$experiment)

human_region_counts$experiment = human_rename_experiments(human_region_counts$experiment)

plot_region_counts(x           = human_ribo,
                   range.lower = MOUSE_MIN_LENGTH,
                   range.upper = MOUSE_MAX_LENGTH)


plot_region_counts(x           = mouse_ribo,
                   range.lower = MOUSE_MIN_LENGTH,
                   range.upper = MOUSE_MAX_LENGTH)

###############################################################################

compute_region_percentages = function(df){
  this_df =
    df %>% 
    group_by(experiment) %>%
    mutate(experiment_total = sum(count))
  
  this_df =
    this_df %>%
    mutate( percentage = round(100 * (count / experiment_total), 1 )  ) %>%
    mutate(region = factor(region, levels = c("UTR3", "CDS", "UTR5") ))
    
  
  
  return(this_df)
}

################################################################################

customized_plot_region_counts = function(df){
  
  percentages <- replace(df$percentage, 
                         df$region != "CDS", "")
  
  this_df = df
  
  this_plot = 
    ggplot(this_df, 
           aes(x=experiment, y=percentage, fill=region)) +
    geom_col() +
    coord_flip() + 
    theme_bw() +
    theme(plot.title   = element_text(hjust = 0.5),
          panel.border = element_blank(),
          panel.grid   = element_blank(),
          axis.text.y       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
          axis.text.x       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
          axis.title.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
          axis.title.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE)) +
    scale_fill_discrete(breaks = c("UTR5", "CDS", "UTR3"), labels = c("5' UTR", "CDS", "3' UTR")) +
    geom_text(aes(x=.data$experiment, y=50, label=percentages), size=3) +
    scale_x_discrete(limits = rev) + 
    labs(title="Ribosome Occupancy", fill="Region", x="Experiment", y="Percentage")
  
  return(this_plot)
}

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

################################################################################
#####           S U P P L E M E N T A R Y     F I G U R E S             ########    

mouse_region_percentages = compute_region_percentages(mouse_region_counts)

human_region_percentages = compute_region_percentages(human_region_counts)

mouse_supplementary_plot = 
  customized_plot_region_counts(mouse_region_percentages)

human_supplementary_plot = 
  customized_plot_region_counts(human_region_percentages)

save_plot_pdf("mouse_region_counts_supp.pdf", mouse_supplementary_plot, width = 5, height = 10)

save_plot_pdf("human_region_counts_supp.pdf", human_supplementary_plot, width = 5, height = 3)

################################################################################


################################################################################
########           A G G R E G A T E D     F I G U R E S             ########### 

mouse_region_percentages = 
  mouse_region_percentages %>%
    mutate( group = unlist(lapply ( strsplit(  as.vector(experiment), split = "-" ), "[[", 1) )  ) %>%
    group_by(group) %>%
    mutate( replicate_count = length( unique( experiment )  ) ) 

mouse_region_percentages$group = factor(mouse_region_percentages$group , levels = c("GV", "MII", "1cell", "2cell", "4cell", "8cell")) 

mouse_region_percentages =
  mouse_region_percentages %>%
  group_by(group, region) %>%
  mutate(average_percentage = mean(percentage) ) %>%
  mutate(standard_error = sd(percentage) / sqrt(replicate_count) )
  

human_region_percentages = 
  human_region_percentages %>%
  mutate( group = unlist(lapply ( strsplit(  as.vector(experiment), split = "-" ), "[[", 1) )  ) %>%
  group_by(group) %>%
  mutate( replicate_count = length( unique( experiment )  ) ) 

human_region_percentages =
  human_region_percentages %>% 
  group_by(group, region) %>%
  mutate(average_percentage = mean(percentage) ) %>%
  mutate(standard_error = sd(percentage) / sqrt(replicate_count) )

plot_bar_with_error_bars = function(df){
  this_df = df
  this_df$region = factor(this_df$region, levels = c("UTR5", "CDS", "UTR3"))
  
  this_plot = 
    ggplot(data=this_df, aes(x= group, y=average_percentage, fill = region )  )  + 
      geom_bar(position = "dodge", stat="identity", alpha = 1 ) + 
      geom_errorbar( aes(  x        = group, 
                           ymin     = average_percentage - standard_error, 
                           ymax     = average_percentage + standard_error), 
                           width    = 0.4, 
                           alpha    = 0.4, 
                           size     = 0.7,
                           position = position_dodge(width = 0.9)) +
      theme(plot.title       = element_text(hjust = 0.5, family = FIGURE_FONT, face = "plain", size = FONT_TITLE_SIZE),
            panel.border     = element_blank(),
            panel.grid       = element_blank(),
            panel.background = element_blank(),
            axis.text.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
            axis.title.y     = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
            axis.text.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
            axis.title.x     = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
            legend.title     = element_blank(),
            legend.text      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE)) + 
      labs(title="Distribution of RPFs by transcript regions", fill="Region", x="Stage", y="Average percentage") + 
      scale_y_continuous(limits = c(0, 100), expand = c(0,0)) + 
      scale_fill_manual("legend", 
                        values = c("UTR5" = UTR5_BLUE, "UTR3" = UTR3_GREEN, "CDS" = CDS_GREEN), 
                        breaks = c("UTR5", "CDS", "UTR3"), 
                        labels = c("5' UTR", "CDS", "3' UTR")) 
  
  return(this_plot)
}

mouse_region_counts_with_error_bars = 
    plot_bar_with_error_bars(mouse_region_percentages)


plot_bar_with_error_bars(human_region_percentages)
