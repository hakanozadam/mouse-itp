#SNP Figure

riboseq_overall_percentage_file = "../snp/snp_dataframes/riboseq_snp_percentages.csv"
rnaseq_overall_percentage_file  = "../snp/snp_dataframes/rnaseq_snp_percentages.csv"

output_folder           = "./pdfs"

################################################################################
########                    L I B R A R I E S                          #########

library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)


################################################################################
#########                  C O L O R I N G                             #########

PERCENTAGE_BARPLOT_PATERNAL_COLOR =  "#999999"
PERCENTAGE_BARPLOT_MATERNAL_COLOR = "#E69F00"
PERCENTAGE_BARPLOT_OTHER_COLOR    = "#56B4E9"

PERCENTAGE_BARPLOT_COLORS = c(PERCENTAGE_BARPLOT_PATERNAL_COLOR, 
                              PERCENTAGE_BARPLOT_MATERNAL_COLOR,
                              PERCENTAGE_BARPLOT_OTHER_COLOR)

MAIN_PERCENTAGE_FILL_COLOR     = "skyblue"
MAIN_PERCENTAGE_ERRORBAR_COLOR = "orange"

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

percentage_averages = 
  rnaseq_overall_percentages %>%
    filter(variable == "Paternal") %>%
    select(Experiment, variable, Percentage) %>%
    group_by( stage = unlist(lapply ( strsplit(  as.vector(Experiment), split = "-" ), "[[", 1) )  ) %>%
    mutate( average_percentage = mean( (Percentage)  )  ) %>%
    filter(stage != "GV") %>%
    select(stage, average_percentage, Percentage) %>%
    summarise(average_percentage = mean( (Percentage)), sd_percentage = sd(Percentage) )
    

ggplot(data=percentage_averages, aes(x=stage, y=average_percentage )  )  + 
  geom_bar(stat="identity", fill = MAIN_PERCENTAGE_FILL_COLOR, alpha = 0.7 ) +
  geom_errorbar( aes(x    = stage, 
                     ymin = average_percentage - sd_percentage, 
                     ymax = average_percentage + sd_percentage), 
                 width  = 0.4, 
                 colour = MAIN_PERCENTAGE_ERRORBAR_COLOR, 
                 alpha  = 0.9, 
                 size   = 1.3) + 
  theme(plot.title       = element_text(hjust = 0.5),
        panel.border     = element_blank(),
        panel.grid       = element_blank(),
        panel.background = element_blank()) + 
  labs(title = "Percentages of paternal alleles", 
       x     = "Percentage", 
       y     = "Stage")


################################################################################

