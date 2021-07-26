#SNP Figure

riboseq_overall_percentage_file = "../snp/snp_dataframes/riboseq_snp_percentages.csv"
rnaseq_overall_percentage_file  = "../snp/snp_dataframes/rnaseq_snp_percentages.csv"

overall_ribo_to_rna_ratios_file = "../snp/snp_dataframes/refined_ribo_to_rna_ratios.csv"
overall_ribo_to_rna_binary_file = "../snp/snp_dataframes/refined_ribo_to_rna_binary_significance.csv"



detailed_riboseq_count_file             = "../snp/snp_dataframes/riboseq_detailed_snps.csv.gz"
detailed_rnaseq_count_file              = "../snp/snp_dataframes/rnaseq_detailed_snps.csv.gz"

output_folder           = "./pdfs"

################################################################################
########                    L I B R A R I E S                          #########

library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
library(pheatmap)
library(cowplot)

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



# rename_experiments = function(raw_names){
#   plain_names = unlist(lapply ( strsplit(  raw_names, split = "-" ), "[[", 3) )
#   
#   uniq_names = unique( plain_names )
#   
#   enumerator = c()
#   
#   for(n in uniq_names){
#     this_count = length( grep( n, plain_names )  )
#     enumerator = c(enumerator, 1:this_count)
#   }
#   
#   actual_names = paste( plain_names, enumerator, sep = "-" ) 
#   
#   return(actual_names)
# }




# rename_experiments_2 = function(raw_names){
#   plain_names = unlist(lapply ( strsplit(  raw_names, split = "-" ), "[[", 3) )
#   
#   uniq_names = unique( plain_names )
#   
#   enumerator = c()
#   
#   for(n in uniq_names){
#     this_count = length( grep( n, uniq_names )  )
#     enumerator = c(enumerator, 1:this_count)
#   }
#   
#   actual_names = paste( plain_names, enumerator, sep = "-" ) 
#   
#   return(actual_names)
# }


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

temp_riboseq_df = read.csv(riboseq_overall_percentage_file, row.names = 1)
temp_rnaseq_df  = read.csv(rnaseq_overall_percentage_file, row.names = 1)


ribo_exp_name_mapper = make_name_mapper( row.names(temp_riboseq_df)  )
rna_exp_name_mapper  = make_name_mapper( row.names(temp_rnaseq_df) )

experiment_name_mapper = c(ribo_exp_name_mapper, rna_exp_name_mapper)

rename_experiments = function(experiment_name){
  return(experiment_name_mapper[[experiment_name]])
}

rename_experiments = Vectorize(rename_experiments, USE.NAMES = FALSE)




################################################################################

rename_genes = function(gene_names, split = "."){
  new_names = unlist(lapply ( strsplit(  gene_names, split = split, fixed = TRUE ), "[[", 1) )
  return(new_names)
}


# rename_genes = function(gene_names){
#   new_names = unlist(lapply ( strsplit(  gene_names, split = ".", fixed = TRUE ), "[[", 1) )
#   return(new_names)
# }


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

heatmap_df_raw             = read.csv(overall_ribo_to_rna_ratios_file, row.names = 1)
colnames( heatmap_df_raw ) = rename_genes(colnames( heatmap_df_raw ) )
heatmap_df                 = t(heatmap_df_raw[-1,])

heatmap_binary_raw             = read.csv(overall_ribo_to_rna_binary_file, row.names = 1)
colnames( heatmap_binary_raw ) = rename_genes(colnames( heatmap_binary_raw ) )
heatmap_binary                 = t(heatmap_binary_raw[-1,])
  
# Adapted from https://stackoverflow.com/questions/31677923/set-0-point-for-pheatmap-in-r
paletteLength <- 100
myColor <- colorRampPalette(c("navy", "white", "red"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
#myBreaks <- c(seq(min(heatmap_df_raw, na.rm = TRUE), 0, length.out=ceiling(paletteLength/2) + 1), 
#              seq(max(heatmap_df_raw, na.rm = TRUE)/paletteLength, max(heatmap_df_raw, na.rm = TRUE), length.out=floor(paletteLength/2)))

myBreaks <- c(seq( -1 , 0, length.out=ceiling(paletteLength/2) + 1), 
              seq( 1 /paletteLength, 1, length.out=floor(paletteLength/2)))


### Trimmed down
### For the main figure
pheatmap( heatmap_df , 
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
pheatmap( heatmap_df , 
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
          display_numbers = heatmap_binary,
          fontsize_number = 11
          )

################################################################################

#### GENE SPECIFIC PLOTS

detailed_riboseq_table_raw = read.csv(detailed_riboseq_count_file)
detailed_rnaseq_table_raw  = read.csv(detailed_rnaseq_count_file)

detailed_riboseq_table_raw$experiment = rename_experiments( detailed_riboseq_table_raw$experiment  )
detailed_riboseq_table_raw$transcript = rename_genes( detailed_riboseq_table_raw$transcript , split = "-" )

detailed_rnaseq_table_raw$experiment = rename_experiments( detailed_rnaseq_table_raw$experiment  )
detailed_rnaseq_table_raw$transcript = rename_genes( detailed_rnaseq_table_raw$transcript , split = "-" )

detailed_riboseq_table_raw$group = unlist(  lapply( strsplit(detailed_riboseq_table_raw$experiment, "-")  , "[[", 1 )  )
detailed_rnaseq_table_raw$group  = unlist(  lapply( strsplit(detailed_rnaseq_table_raw$experiment, "-")  , "[[", 1 )  )

experiment_groups = unique( detailed_riboseq_table_raw$group )
  
experiment_name_mapper

this_gene = "Wtap"



this_group = "8cell"

gene_riboseq_group_data = gene_riboseq_raw_data %>% filter( group == this_group )


gene_rnaseq_group_data = gene_rnaseq_raw_data %>% filter( group == this_group )


##############################################################

get_gene_percentages = function(this_df, this_gene){
  ## Add percentage column to the table
  
  gene_data = this_df %>% 
                          filter(transcript == this_gene) %>%
                          group_by(  experiment, position ) %>%
                          mutate( maternal_total = sum(maternal), paternal_total = sum(paternal) ) %>%
                          mutate( paternal_percentage = (paternal / (paternal_total + maternal_total )) * 100  )

  return(gene_data)
}

##############################################################

gene_percentage_ribo_data = get_gene_percentages( detailed_riboseq_table_raw, "Wtap" ) 


barplot_snp_counts_of_gene_helper = function( gene, df, exp_group, allele_type, ymax  ){
  gene_data = df %>% filter(transcript == gene & group == exp_group ) %>% arrange( position  )
  
  barplot_colors = colorRampPalette(c("navy", "yellow", "red"))( length( unique(gene_data$position  ) ) )
  
  if(allele_type == "paternal") { 
    p =  ggplot(data=gene_data %>% filter(group==exp_group), 
           aes(x    = experiment, 
               y    = paternal, 
               fill = factor(position, levels = sort(unique( gene_data$position  ) )  )  )  )
  }
  else{
    p =  ggplot(data=gene_data %>% filter(group==exp_group), 
                aes(x    = experiment, 
                    y    = maternal, 
                    fill = factor(position, levels = sort(unique( gene_data$position  ) )  )  )  )
  }

  
  this_plot = p +
    geom_bar(stat="identity" ) +
    theme(
      panel.border     = element_blank(),
      panel.grid       = element_blank(),
      panel.background = element_blank(),
      axis.title.x     = element_blank(),
      legend.position  = "none",
      plot.margin      = margin(0, 0, 0, 0, "cm"),
      # axis.title.y     = element_blank(),
      # axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x=element_blank(),
      # axis.ticks.y=element_blank()
    ) +
    #scale_fill_brewer(palette = "Blues")
    scale_fill_manual(values = barplot_colors) + 
    ylim(c(0, ymax))
  
   
  return(this_plot)
}

################################################################################

barplot_snp_counts_of_gene = function(gene, exp_group, df){
  
  this_gene      = gene
  this_exp_group = exp_group
  
  max_paternal = max(
                 ( df %>% 
                     filter(transcript == this_gene & group == this_exp_group)  %>%
                     group_by(experiment) %>%
                     mutate(paternal_sum = sum(paternal))
                    )$paternal_sum 
                 )
  
  max_maternal = max(
    ( df %>% 
        filter(transcript == this_gene & group == this_exp_group)  %>%
        group_by(experiment) %>%
        mutate(maternal_sum = sum(maternal))
    )$maternal_sum 
  )
  
  this_ymax = max( max_maternal, max_paternal )
  
  df = df %>% filter(transcript == this_gene & group == this_exp_group) %>% arrange( position )
  
  
  p1 = barplot_snp_counts_of_gene_helper(gene        = this_gene, 
                                         df          = df, 
                                         exp_group   = this_exp_group, 
                                         allele_type = "paternal",
                                         ymax        = this_ymax) 
  
  p2 = barplot_snp_counts_of_gene_helper(gene        = this_gene, 
                                         df          = df, 
                                         exp_group   = this_exp_group, 
                                         allele_type = "maternal",
                                         ymax        = this_ymax)
  p2 = p2 + scale_y_reverse(limits = (c(this_ymax, 0)))
  
  result_p = plot_grid( p1, p2, ncol = 1   )
  
  return(result_p)
}
 

barplot_snp_counts_of_gene(gene = "Wtap", exp_group = "8cell", df = detailed_riboseq_table_raw)

##############################################################

# 
# wtap_ribo_data = detailed_riboseq_table_raw %>% filter(transcript == "Wtap")
# 
# wtap_rna_data = detailed_rnaseq_table_raw %>% filter(transcript == "Wtap")
# 
# find_paternal_ratio_in_expgroup = function( df, exp_group ){
#   aggregated_counts = 
#     df %>% filter(group == exp_group) %>% 
#     summarise( paternal_total = sum(paternal), maternal_total = sum(maternal) ) 
#   
#   paternal_ratio = aggregated_counts$paternal_total / (aggregated_counts$paternal_total + aggregated_counts$maternal_total)
#   
#   return(paternal_ratio)
# }
# 
# 
# find_paternal_ratio_in_expgroup = function( df ){
#   paternal_ratios = 
#     df %>% group_by(group) %>% 
#     summarise( paternal_total = sum(paternal), maternal_total = sum(maternal) ) %>%
#     mutate(paternal_ratio = paternal_total / (paternal_total + maternal_total) ) %>%
#     filter( group == "2cell" | group == "4cell" | group == "8cell") %>% arrange( group )
#   
#   return(paternal_ratios)
# }
colorRampPalette(c("navy", "white", "red"))(3)
