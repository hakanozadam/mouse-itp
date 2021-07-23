#SNP Figure

riboseq_overall_percentage_file = "../snp/snp_dataframes/riboseq_snp_percentages.csv"
rnaseq_overall_percentage_file  = "../snp/snp_dataframes/rnaseq_snp_percentages.csv"

overall_ribo_to_rna_ratios_file = "../snp/snp_dataframes/refined_ribo_to_rna_ratios.csv"
overall_ribo_to_rna_binary_file = "../snp/snp_dataframes/refined_ribo_to_rna_binary_significance.csv"

output_folder           = "./pdfs"

################################################################################
########                    L I B R A R I E S                          #########

library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
library(pheatmap)

################################################################################
#########                  C O L O R I N G                             #########

PERCENTAGE_BARPLOT_PATERNAL_COLOR =  "#999999"
PERCENTAGE_BARPLOT_MATERNAL_COLOR = "#E69F00"
PERCENTAGE_BARPLOT_OTHER_COLOR    = "#56B4E9"

PERCENTAGE_BARPLOT_COLORS = c(PERCENTAGE_BARPLOT_PATERNAL_COLOR, 
                              PERCENTAGE_BARPLOT_MATERNAL_COLOR,
                              PERCENTAGE_BARPLOT_OTHER_COLOR)

MAIN_PERCENTAGE_RNASEQ_FILL_COLOR     = "skyblue"
MAIN_PERCENTAGE_RIBOSEQ_FILL_COLOR    = "#abe39d"
MAIN_PERCENTAGE_ERRORBAR_COLOR = "#413bb8"

################################################################################
##### Experiments



rename_experiments = function(raw_names){
  plain_names = unlist(lapply ( strsplit(  raw_names, split = "-" ), "[[", 3) )
  
  uniq_names = unique( plain_names )
  
  enumerator = c()
  
  for(n in uniq_names){
    this_count = length( grep( n, plain_names )  )
    enumerator = c(enumerator, 1:this_count)
  }
  
  actual_names = paste( plain_names, enumerator, sep = "-" ) 
  
  return(actual_names)
}
################################################################################

##### PERCENTAGES SUPPLEMENTARY

riboseq_overall_percentages_pre = read.csv( riboseq_overall_percentage_file, 
                                        row.names = 1 )

rnaseq_overall_percentages_pre  = read.csv( rnaseq_overall_percentage_file, 
                                        row.names = 1 )

riboseq_raw_experiment_names               = row.names(riboseq_overall_percentages_pre)
riboseq_overall_percentages_pre$Experiment = rename_experiments(riboseq_raw_experiment_names)

rnaseq_raw_experiment_names               = row.names(rnaseq_overall_percentages_pre)
rnaseq_overall_percentages_pre$Experiment = rename_experiments(rnaseq_raw_experiment_names)

# Prep the Data for ggplot
riboseq_overall_percentages = melt(riboseq_overall_percentages_pre, value.name = "Percentage")
rnaseq_overall_percentages  = melt(rnaseq_overall_percentages_pre, value.name = "Percentage")

# Reorder the columns so that 1-cell appears on top
riboseq_overall_percentages$Experiment = 
    factor(riboseq_overall_percentages$Experiment, 
           levels = unique(rev( riboseq_overall_percentages$Experiment )) )

rnaseq_overall_percentages$Experiment = 
  factor(rnaseq_overall_percentages$Experiment, 
         levels = unique(rev( rnaseq_overall_percentages$Experiment )) )

################################################################################
## SUPPLEMENTARY FIGURE FOR DETAILED ALLELE PERCENTAGES ###############################

# Ribosome Profiling

ggplot(data=riboseq_overall_percentages, aes(x=Experiment, y=Percentage, fill = factor(variable, levels = c("Other", "Maternal", "Paternal")  ) ) )  +
  geom_bar(stat="identity" ) +
  coord_flip() + 
  guides(fill=guide_legend(title="Type")) + 
  theme(plot.title       = element_text(hjust = 0.5),
        panel.border     = element_blank(),
        panel.grid       = element_blank(),
        panel.background = element_blank())  + 
  scale_fill_manual(values = PERCENTAGE_BARPLOT_COLORS ) + 
  geom_hline(yintercept = 50, linetype="dotted", 
             color = "brown", size=0.6) + 
  labs(title = "Allele Percentages")


# RNA-Seq

ggplot(data=rnaseq_overall_percentages, aes(x=Experiment, y=Percentage, fill = factor(variable, levels = c("Other", "Maternal", "Paternal")  ) ) )  +
  geom_bar(stat="identity" ) +
  coord_flip() + 
  guides(fill=guide_legend(title="Type")) + 
  theme(plot.title       = element_text(hjust = 0.5),
        panel.border     = element_blank(),
        panel.grid       = element_blank(),
        panel.background = element_blank())  + 
  scale_fill_manual(values = PERCENTAGE_BARPLOT_COLORS ) + 
  geom_hline(yintercept = 50, linetype="dotted", 
             color = "brown", size=0.6) + 
  labs(title = "Allele Percentages")


################################################################################
#####       M A I N    F I G U R E   P E R C E N T A G E S      ################

get_percentage_averages = function( input_df, label ){
  result = 
    input_df %>%
    filter(variable == "Paternal") %>%
    select(Experiment, variable, Percentage) %>%
    group_by( stage = unlist(lapply ( strsplit(  as.vector(Experiment), split = "-" ), "[[", 1) )  ) %>%
    mutate( average_percentage = mean( (Percentage)  )  ) %>%
    filter(stage != "GV") %>%
    select(stage, average_percentage, Percentage) %>%
    summarise(average_percentage = mean( (Percentage)), sd_percentage = sd(Percentage) ) %>%
    mutate(experiment_type = label)
  
  return(result)
}

rnaseq_percentage_averages  = get_percentage_averages(rnaseq_overall_percentages, "RNA")
riboseq_percentage_averages = get_percentage_averages(riboseq_overall_percentages, "RIBO")    

percentage_averages = bind_rows( riboseq_percentage_averages, rnaseq_percentage_averages  )

percentage_averages$stage = factor( percentage_averages$stage, levels = c("MII", "1cell", "2cell", "4cell", "8cell") )

ggplot(data=percentage_averages, aes(x=stage, y=average_percentage, fill = experiment_type )  )  + 
  geom_bar(position = "dodge", stat="identity", alpha = 0.7 ) +
  geom_errorbar( aes(  x    = stage, 
                       ymin = average_percentage - sd_percentage, 
                       ymax = average_percentage + sd_percentage), 
                 width    = 0.4, 
                 colour   = MAIN_PERCENTAGE_ERRORBAR_COLOR, 
                 alpha    = 0.9, 
                 size     = 1.3,
                 position = position_dodge(width = 0.8)) + 
  
  theme(plot.title       = element_text(hjust = 0.5),
        panel.border     = element_blank(),
        panel.grid       = element_blank(),
        panel.background = element_blank()) + 
  
  scale_fill_manual( values = c(MAIN_PERCENTAGE_RIBOSEQ_FILL_COLOR,
                                MAIN_PERCENTAGE_RNASEQ_FILL_COLOR) ) + 
  
  guides(fill=guide_legend(title="")) + 
  
  labs(title = "Percentages of paternal alleles", 
       x     = "Percentage", 
       y     = "Stage")


################################################################################


################################################################################
#####################      H E A T M A P      ##################################

heatmap_df_raw = read.csv(overall_ribo_to_rna_ratios_file, row.names = 1)

# Adapted from https://stackoverflow.com/questions/31677923/set-0-point-for-pheatmap-in-r
paletteLength <- 100
myColor <- colorRampPalette(c("navy", "white", "red"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(heatmap_df_raw, na.rm = TRUE), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(heatmap_df_raw, na.rm = TRUE)/paletteLength, max(heatmap_df_raw, na.rm = TRUE), length.out=floor(paletteLength/2)))

myBreaks <- c(seq( -1 , 0, length.out=ceiling(paletteLength/2) + 1), 
              seq( 1 /paletteLength, 1, length.out=floor(paletteLength/2)))


pheatmap( heatmap_df_raw[-1,], cluster_rows = FALSE, 
          color= myColor,
          breaks = myBreaks,
          na_col = "white",
          labels_row = c("1cell", "2cell", "4cell", "8cell"))

### Trimmed down
### For the main figure
pheatmap( t(heatmap_df_raw[-1,]) , 
          #clustering_distance_rows = "correlation",
          #clustering_method = "centroid",
          #clustering_method = "complete",
          #clustering_method = "ward.D",
          show_rownames = FALSE,
          clustering_method = "median",
          cellwidth = 40,
          treeheight_row = 0,
          treeheight_col = 0,
          cluster_cols = FALSE, 
          color= myColor,
          breaks = myBreaks,
          na_col = "white",
          angle_col = 0,
          labels_col = c("1cell", "2cell", "4cell", "8cell"),
          main = "Paternal Ratio Difference",
          )


#### Detailed with row names
#### For the supplementary Figure
pheatmap( t(heatmap_df_raw[-1,]) , 
          #clustering_distance_rows = "correlation",
          #clustering_method = "centroid",
          #clustering_method = "complete",
          #clustering_method = "ward.D",
          legend_breaks = c(-1,  -0.5, 0, 0.5, 1),
          legend = T,
          show_rownames = TRUE,
          clustering_method = "median",
          cellwidth = 40,
          treeheight_row = 60,
          cluster_cols = FALSE, 
          color= myColor,
          breaks = myBreaks,
          na_col = "white",
          angle_col = 0,
          labels_col = c("1cell", "2cell", "4cell", "8cell"),
          main = "Paternal Ratio Difference",
          )

################################################################################



