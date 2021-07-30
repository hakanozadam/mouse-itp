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
library(RColorBrewer)

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

FOURCELL_BACKGROUND_COLOR = "#d9d9d9"

COUNT_NORMALIZATION_FACTOR = 10000

RATIO_RIBO_COLOUR = "orange"
RATIO_RNA_COLOUR  = "blue"


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
       x     = "Stage", 
       y     = "Percentage")


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
  
## Add counts per thousand to the dataframes
detailed_riboseq_table = 
  detailed_riboseq_table_raw %>% 
  group_by(experiment) %>%
  mutate( experiment_total = sum(paternal + maternal) ) %>%
  mutate( paternal_per_k = (paternal / experiment_total)* COUNT_NORMALIZATION_FACTOR,
          maternal_per_k = (maternal / experiment_total)*COUNT_NORMALIZATION_FACTOR ) %>%
  mutate(type = "ribo")


detailed_rnaseq_table = 
  detailed_rnaseq_table_raw %>% 
  group_by(experiment) %>%
  mutate( experiment_total = sum(paternal + maternal) ) %>%
  mutate( paternal_per_k = (paternal / experiment_total)* COUNT_NORMALIZATION_FACTOR,
          maternal_per_k = (maternal / experiment_total)*COUNT_NORMALIZATION_FACTOR ) %>%
  mutate(type = "rna")


detailed_table = rbind(detailed_riboseq_table, detailed_rnaseq_table)

###
# To handle graphs easier, we add dummy experiments to 2,4,8 cell groups
# so that they have equal (4) number of experiments.
# The new experiments have 0 counts.
detailed_riboseq_table_with_dummies =
  detailed_riboseq_table %>%
  filter(group == "2cell" | group == "4cell" | group == "8cell")

tmp_dummy_entry = detailed_riboseq_table %>% 
                    filter( experiment == "2cell-1" ) %>%
                    mutate( experiment = "2cell-4" ) %>%
                    mutate( paternal = 0, maternal = 0 , paternal_per_k = 0, maternal_per_k = 0 )

detailed_riboseq_table_with_dummies = rbind(detailed_riboseq_table_with_dummies, tmp_dummy_entry)


tmp_dummy_entry = detailed_riboseq_table %>% 
  filter( experiment == "4cell-1" ) %>%
  mutate( experiment = "4cell-4" ) %>%
  mutate( paternal = 0, maternal = 0 , paternal_per_k = 0, maternal_per_k = 0 )

detailed_riboseq_table_with_dummies = rbind(detailed_riboseq_table_with_dummies, tmp_dummy_entry)



detailed_rnaseq_table_with_dummies =
  detailed_rnaseq_table %>%
  filter(group == "2cell" | group == "4cell" | group == "8cell")

tmp_dummy_entry = detailed_rnaseq_table %>%
  filter(experiment == "4cell-1") %>%
  mutate(experiment = "4cell-3") %>%
  mutate( paternal = 0, maternal = 0 , paternal_per_k = 0, maternal_per_k = 0 )

detailed_rnaseq_table_with_dummies = rbind( detailed_rnaseq_table_with_dummies, tmp_dummy_entry )

tmp_dummy_entry = detailed_rnaseq_table %>%
  filter(experiment == "4cell-1") %>%
  mutate(experiment = "4cell-4") %>%
  mutate( paternal = 0, maternal = 0 , paternal_per_k = 0, maternal_per_k = 0 )

detailed_rnaseq_table_with_dummies = rbind( detailed_rnaseq_table_with_dummies, tmp_dummy_entry )

detailed_table_with_dummies = rbind(  detailed_riboseq_table_with_dummies, detailed_rnaseq_table_with_dummies )
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

### To be deleted
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
### TO be deleted 
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
 

barplot_snp_counts_of_gene(gene = "Nin", exp_group = "8cell", df = detailed_riboseq_table_raw)

##############################################################
##############################################################
##############################################################


experiment_snp_counts_of_gene_helper = function( gene, df, exp_group, allele_type, ymax  ){
  gene_data = df %>% filter(transcript == gene & group == exp_group ) %>% arrange( position  )
  
  barplot_colors = colorRampPalette(c("navy", "yellow", "red"))( length( unique(gene_data$position  ) ) )
  
  if(allele_type == "paternal") { 
    p =  ggplot(data=gene_data %>% filter(group==exp_group), 
                aes(x    = experiment, 
                    y    = paternal_per_k, 
                    fill = factor(position, levels = sort(unique( gene_data$position  ) )  )  )  )
  }
  else{
    p =  ggplot(data=gene_data %>% filter(group==exp_group), 
                aes(x    = experiment, 
                    y    = maternal_per_k, 
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
################################################################################
################################################################################

experiment_snp_counts_of_gene = function(gene, exp_group, df){
  
  this_gene      = gene
  this_exp_group = exp_group
  
  max_paternal = max(
    ( df %>% 
        filter(transcript == this_gene & group == this_exp_group)  %>%
        group_by(experiment) %>%
        mutate(paternal_sum = sum(paternal_per_k))
    )$paternal_sum 
  )
  
  max_maternal = max(
    ( df %>% 
        filter(transcript == this_gene & group == this_exp_group)  %>%
        group_by(experiment) %>%
        mutate(maternal_sum = sum(maternal_per_k))
    )$maternal_sum 
  )
  
  this_ymax = max( max_maternal, max_paternal )
  
  df = df %>% filter(transcript == this_gene & group == this_exp_group) %>% arrange( position )
  
  
  p1 = experiment_snp_counts_of_gene_helper(gene        = this_gene, 
                                         df          = df, 
                                         exp_group   = this_exp_group, 
                                         allele_type = "paternal",
                                         ymax        = this_ymax) 
  
  p2 = experiment_snp_counts_of_gene_helper(gene        = this_gene, 
                                         df          = df, 
                                         exp_group   = this_exp_group, 
                                         allele_type = "maternal",
                                         ymax        = this_ymax)
  p2 = p2 + scale_y_reverse(limits = (c(this_ymax, 0)))
  
  result_p = plot_grid( p1, p2, ncol = 1   )
  
  return(result_p)
}

################################################################################

#combined_snp_counts_of_gene = function(gene, ribo_df, rnaseq_df ){
#  
#}



p2_ribo = experiment_snp_counts_of_gene(gene = "Mysm1", exp_group = "2cell", df = detailed_riboseq_table)
p2_rna  = experiment_snp_counts_of_gene(gene = "Mysm1", exp_group = "2cell", df = detailed_rnaseq_table)

p4_ribo = experiment_snp_counts_of_gene(gene = "Mysm1", exp_group = "4cell", df = detailed_riboseq_table)
p4_rna  = experiment_snp_counts_of_gene(gene = "Mysm1", exp_group = "4cell", df = detailed_rnaseq_table)

p8_ribo = experiment_snp_counts_of_gene(gene = "Mysm1", exp_group = "8cell", df = detailed_riboseq_table)
p8_rna  = experiment_snp_counts_of_gene(gene = "Mysm1", exp_group = "8cell", df = detailed_rnaseq_table)


p4_ribo

pdf("scale_problem.pdf")
plot_grid( p2_ribo, p2_rna,
           p4_ribo, p4_rna,
           p8_ribo, p8_rna,
           nrow = 1   )
dev.off()

View(
  detailed_rnaseq_table %>% 
  filter(group == "4cell") %>%
  filter(transcript == "Wtap")
)


raw_p2_ribo = barplot_snp_counts_of_gene(gene = "Mysm1", exp_group = "2cell", df = detailed_riboseq_table)
raw_p2_rna = barplot_snp_counts_of_gene(gene = "Mysm1", exp_group = "2cell", df = detailed_rnaseq_table)

raw_p4_ribo = barplot_snp_counts_of_gene(gene = "Mysm1", exp_group = "4cell", df = detailed_riboseq_table)
raw_p4_rna = barplot_snp_counts_of_gene(gene = "Mysm1", exp_group = "4cell", df = detailed_rnaseq_table)

raw_p8_ribo = barplot_snp_counts_of_gene(gene = "Mysm1", exp_group = "8cell", df = detailed_riboseq_table)
raw_p8_rna = barplot_snp_counts_of_gene(gene = "Mysm1", exp_group = "8cell", df = detailed_rnaseq_table)

plot_grid( raw_p2_ribo, raw_p2_rna,
           raw_p4_ribo, raw_p4_rna,
           raw_p8_ribo,raw_p8_rna,
           nrow = 1   )


##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################

generate_blank_plot = function(bg = "white"){
  
  this_background = element_blank()
    
  if(bg != "white"){
    this_background = element_rect(fill = bg)
  }  
  
  df <- data.frame()
  p = ggplot(df) + geom_point() + 
    xlim(0, 4) + ylim(0, 100) + 
    theme(
      panel.border     = element_blank(),
      panel.grid       = element_blank(),
      panel.background = this_background,
      axis.title.x     = element_blank(),
      axis.title.y     = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x=element_blank(),
      plot.background = element_blank(),
      axis.ticks.y=element_blank()
    ) 
  return(p)
}

format_y_axis = function(x){
  return( sprintf("%.0f", x) )
}

make_y_axis = function(ymax){
  this_multiplier    = floor(ymax / 5)
  
  if( (ymax %% 5) >= 1 ) {
    this_multiplier = this_multiplier + 1
  } 
  
  yticks = seq(0, 5*this_multiplier, 5)
  
  return(yticks)
  
}

adjust_ymax = function( y ){
  this_multiplier    = floor(y / 5)
  
  if(y == 0){
    return(5)
  }
  
  if( (y %% 5) > 0 ) {
    this_multiplier = this_multiplier + 1
  } 
  
  return( 5*this_multiplier)
}

################################################################################


generate_legend_unit = function(gene, experiment_type){
  
  gene_data = 
    detailed_table_with_dummies %>%
    filter( transcript == gene & type == experiment_type ) %>%
    filter( group == "2cell" | group == "4cell" | group == "8cell")
  
  number_of_snps = length( unique(gene_data$position  ) )
  
  bar_palettes = list( ribo = "Oranges", rna = "Blues"  )
  title_colors = list( ribo = "orange", rna = "blue"  )
  title_texts  = list( ribo = "ribo", rna = "rna"  )
  
  barplot_colors = colorRampPalette(brewer.pal(9,bar_palettes[[experiment_type]]))( number_of_snps + 5 )[-seq(1,5)] 
  
  df = data.frame( x = rep(1,number_of_snps), y = rep(1,number_of_snps), position = seq(1,number_of_snps))
  
  p = ggplot(data=df, 
             aes(x    = x, 
                 y    = y, 
                 fill = factor(position, levels = sort(unique( df$position  ) )  )  )  )  + 
    geom_bar(stat="identity", width= 0.8, color = "black") + 
    theme(
      panel.border     = element_blank(),
      panel.grid       = element_blank(),
      panel.background = element_blank(),
      axis.title.x     = element_blank(),
      #axis.title.x     = element_text(family = "helvetica", face = "plain", size = 10),
      legend.position  = "none",
      axis.title.y     = element_blank(),
      axis.text.x = element_blank(),
      #axis.text.x = element_text(family = "helvetica", face = "plain", size = 10),
      axis.ticks.x=element_blank(),
      axis.ticks.y=element_blank(),
      plot.background = element_blank(),
      axis.text.y = element_blank(),
      plot.title = element_text(color = title_colors[[experiment_type]], size=10, face="bold", hjust = 0.5),
    ) + 
    scale_fill_manual(values = barplot_colors)  + 
    ggtitle(title_texts[[experiment_type]])
  
  #if(experiment_type == "rna"){
  #  p = p + scale_y_reverse()
  #}
  
  return(p)
}

generate_legend = function(gene){
  p_ribo = generate_legend_unit(gene, experiment_type = "ribo")
  p_rna  = generate_legend_unit(gene, experiment_type = "rna")
  
  this_blank_pad  = generate_blank_plot()
  this_legend_pre = plot_grid( this_blank_pad, p_ribo, this_blank_pad, p_rna, 
                           nrow = 1,
                           rel_widths = c(0.2, 1 , 0.2, 1 ) )
  this_legend     = plot_grid(this_blank_pad, this_legend_pre, this_blank_pad,
                              ncol = 1, rel_heights = c(1, 0.6, 1)) 
  
  return(this_legend)
}

snp_detailed_legend = generate_legend("Nin")

snp_detailed_legend
generate_legend("Nin")

################################################################################

################################################################################
### U N I T     P L O T 

unit_plot = function(gene, experiment_type, experiment_group, allele_type, ymax){
  
  bar_palettes = list( ribo = "Oranges", rna = "Blues"  )
  
  y_scale_expand = c(0,0)
  
  ymax = adjust_ymax(ymax)
  
  gene_data = 
    detailed_table_with_dummies %>%
    filter( transcript == gene & type == experiment_type ) %>%
    filter( group == experiment_group)
  
  number_of_snps = length( unique(gene_data$position  ) )
  
  paternal_color_palette = colorRampPalette(brewer.pal(9,bar_palettes[[experiment_type]]))( number_of_snps + 5 )[-seq( 1,  5) ]
  #paternal_color_palette = colorRampPalette(brewer.pal(9,bar_palettes[[experiment_type]]))( number_of_snps  ) 
  barplot_colors         = paternal_color_palette
  
  yticks = make_y_axis(ymax)
  
  if(experiment_group == "4cell"){
    background_canvas = element_rect(fill = FOURCELL_BACKGROUND_COLOR)
  }
  else{
    background_canvas = element_rect(fill = "white")
  }
  
  
  if(allele_type == "paternal"){
    plot_base = ggplot(data=gene_data, 
                       aes(x    = experiment, 
                           y    = paternal_per_k, 
                           fill = factor(position, levels = sort(unique( gene_data$position  ) )  )  )  ) 
  }else{
    plot_base = ggplot(data=gene_data, 
                       aes(x    = experiment, 
                           y    = maternal_per_k, 
                           fill = factor(position, levels = sort(unique( gene_data$position  ) )  )  )  )
  }
  
  this_plot = 
    plot_base + 
    geom_bar(stat="identity", width= 1) + 
    theme(
      panel.border     = element_blank(),
      panel.grid       = element_blank(),
      panel.background = element_blank(),
      axis.title.x     = element_blank(),
      legend.position  = "none",
      # plot.margin      = margin(5.5, 5.5, 0, 5.5),
      axis.title.y     = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x=element_blank(),
      plot.background = background_canvas,
      axis.text.y = element_text(family = "helvetica", face = "plain", size = 10),
      # axis.ticks.y=element_blank()
    ) +
    scale_fill_manual(values = barplot_colors) 
  
  if(experiment_type == "ribo"){
    if(allele_type == "maternal"){
      this_plot = this_plot + scale_y_reverse(position = "right", 
                                              labels   = format_y_axis, 
                                              breaks   = yticks,
                                              limits   = c(ymax, 0) , 
                                              expand   = y_scale_expand ) 
    }
    else{
      this_plot = this_plot + scale_y_continuous(position = "right", 
                                                 labels   = format_y_axis, 
                                                 breaks   = yticks,
                                                 limits   = c(0, ymax) , 
                                                 expand   = y_scale_expand )
    }
  } 
  else{
    if(allele_type == "maternal"){
      this_plot = this_plot + scale_y_reverse(labels = format_y_axis, 
                                              limits = c(ymax, 0), 
                                              breaks   = yticks,
                                              expand = y_scale_expand )
    }
    else{
      this_plot = this_plot + scale_y_continuous( labels = format_y_axis, 
                                                  limits = c(0, ymax),
                                                  breaks   = yticks,
                                                  expand = y_scale_expand  )
    }
  }
  
  return(this_plot)
}


paternal_maternal_unit_plot = function( gene, experiment_group, experiment_type ){
  
  gene_sums = 
    detailed_table_with_dummies %>%
    filter( transcript == gene & type == experiment_type ) %>%
    filter( group == experiment_group) %>%
    group_by(experiment) %>%
    mutate(paternal_sum = sum(paternal_per_k), maternal_sum = sum(maternal_per_k))

  ymax = max( max(gene_sums$paternal_sum), max( gene_sums$maternal_sum ) )
  
  paternal_plot= unit_plot(gene = gene, experiment_type = experiment_type, experiment_group = experiment_group, allele_type = "paternal", ymax = ymax)
  maternal_plot= unit_plot(gene = gene, experiment_type = experiment_type, experiment_group = experiment_group, allele_type = "maternal", ymax = ymax)
  
  this_plot = plot_grid(paternal_plot, maternal_plot, ncol = 1)
  
  return(this_plot)
}

paternal_maternal_unit_plot(gene = "Nin", experiment_group = "8cell", experiment_type = "rna")


plot_experiment_group = function(gene, experiment_group){
  ribo_plot = paternal_maternal_unit_plot(gene = gene, experiment_group = experiment_group, experiment_type = "ribo")
  rna_plot  = paternal_maternal_unit_plot(gene = gene, experiment_group = experiment_group, experiment_type = "rna")
  
  this_plot = plot_grid(  rna_plot, ribo_plot, ncol = 2 , align = "h", axis = "b" )
  
  return(this_plot)
}

plot_experiment_group(gene = "Nin", experiment_group = "8cell")


################################################################################

plot_snp_gene_detailed = function(gene){
# This is the main function that plots the snps in detail  

  
  plot_2cell_raw = plot_experiment_group(gene = gene, experiment_group = "2cell")
  
  title_2cell = ggdraw() + 
    draw_label(
      "2cell",
      fontfamily = 'helvetica',
      size = 10
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 0)
    )
  
  plot_2cell = plot_grid(
    title_2cell, plot_2cell_raw,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.05, 2)
  )
  
  
  plot_4cell_raw =  plot_experiment_group(gene = gene, experiment_group = "4cell")  

  title_4cell = ggdraw() + 
    draw_label(
      "4cell",
      fontfamily = 'helvetica',
      size = 10
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 0)
    )
  
  plot_4cell = plot_grid(
    title_4cell, plot_4cell_raw,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.05, 2)
  )  
  
  plot_8cell_raw = plot_experiment_group(gene = gene, experiment_group = "8cell")

  title_8cell = ggdraw() + 
    draw_label(
      "8cell",
      fontfamily = 'helvetica',
      size = 10
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 0)
    )
  
  plot_8cell = plot_grid(
    title_8cell, plot_8cell_raw,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.05, 2)
  )  
    
  separator_1 = generate_blank_plot()
  
  snp_detailed_legend = generate_legend(gene)
  
  raw_plot = plot_grid( plot_2cell, separator_1, plot_4cell, separator_1, plot_8cell, separator_1, snp_detailed_legend,
                        rel_widths = c(1, 0.2, 1, 0.2, 1, 0.1, 0.4)  , 
                        rel_heights =  c(1, 1, 1, 1, 1, 1, 0.01),
                        ncol       = 7  )
  
  title_main = ggdraw() + 
    draw_label(
      gene,
      fontface = 'bold',
      fontfamily = 'helvetica',
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 0)
    )
  
  this_plot = plot_grid(
    title_main, raw_plot,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.05, 1)
  )
  
  return(this_plot)
}

plot_snp_gene_detailed(gene = "Rpl38")  


################################################################################
################################################################################
################################################################################
################################################################################

compute_sd_of_ratios = function(paternal_counts, total_counts){
  paternal_mean = mean(paternal_counts)
  total_mean    = mean(total_counts)
  
  f = (paternal_mean / total_mean)
  
  this_cov = cov( paternal_counts, total_counts )
  
  paternal_sd = sd(paternal_counts)
  total_sd    = sd(total_counts)
  
  this_sd =  f * sqrt(  (paternal_sd / paternal_mean)**2 + 
                         (total_sd/ total_mean)**2 - 
                         2*(this_cov/ (paternal_mean* total_mean) )  )
  
  return(this_sd)
}


determine_propogated_sd = function(df){
  pre_data = df %>%  group_by(group, experiment) %>%
    summarise(experiment_total = sum(paternal_per_k + maternal_per_k) , experiment_paternal_total = sum(paternal_per_k))
  
  prop_sd = pre_data %>% group_by(group) %>%
    summarise( group_sd =  compute_sd_of_ratios( experiment_paternal_total, experiment_total )  ) %>%
    arrange(group)
  
  return(prop_sd$group_sd)
}



################################################################################

find_paternal_ratios = function(df){
 
  ratios = 
    df %>%
    group_by( group ) %>%
    summarise( paternal_total = sum(paternal), 
               maternal_total = sum(maternal), 
               overall_total  = paternal_total + maternal_total, 
               paternal_ratio = paternal_total / overall_total ) 
  
  return( ratios$paternal_ratio )
}




find_paternal_per_k_ratios = function(df){
  
  ratios = 
    df %>%
    group_by( group ) %>%
    summarise( paternal_total = sum(paternal_per_k), 
               maternal_total = sum(maternal_per_k), 
               overall_total  = paternal_total + maternal_total, 
               paternal_ratio = paternal_total / overall_total ) %>%
    arrange(group)
  
  return( ratios$paternal_ratio )
}

################################################################################

plot_snp_ratios = function(gene, ymax = 0){
  gene_data_main = 
    detailed_table %>% 
    filter( transcript == gene ) %>%
    filter ( group == "2cell" | group == "4cell" | group == "8cell"  )
  
  gene_data_ribo = 
    gene_data_main %>%
    filter( type == "ribo")
  
  gene_data_rna = 
    gene_data_main %>%
    filter( type == "rna")
  
  ribo_ratios = find_paternal_per_k_ratios( gene_data_ribo )
  rna_ratios  = find_paternal_per_k_ratios( gene_data_rna )
  
  ribo_sd = determine_propogated_sd( gene_data_ribo  )
  rna_sd  = determine_propogated_sd(  gene_data_rna )
  
  plot_df = data.frame(  paternal_ratios = c(ribo_ratios, rna_ratios), 
                         stage           = rep( c("2cell", "4cell", "8cell") , 2), 
                         sd              = c(ribo_sd, rna_sd), 
                         type =  c( rep("ribo", 3), rep("rna", 3)  )  )
  
  p = ggplot(data=plot_df, 
             aes(x    = stage, y = paternal_ratios, group = type )  )  +
    geom_point( aes(colour = type, shape = type), position = position_dodge(width = 0.1), size = 4 ) +
    geom_line(aes(linetype=type, colour = type), position = position_dodge(width = 0.1) ) + 
    scale_linetype_manual(values=c("solid", "dashed")) + 
    geom_errorbar(aes(ymin=paternal_ratios-sd, ymax=paternal_ratios+sd, color = type), position = position_dodge(width = 0.1), width = 0.1 ) + 
    scale_colour_manual(values = c(RATIO_RIBO_COLOUR, RATIO_RNA_COLOUR)   ) + 
    theme(
      panel.border      = element_blank(),
      panel.grid        = element_blank(),
      plot.title        = element_text(hjust = 0.5, family = "helvetica", face = "bold", size = 12),
      panel.background  = element_blank(),
      axis.text.y       = element_text(family = "helvetica", face = "plain", size = 8),
      axis.text.x       = element_text(family = "helvetica", face = "plain", size = 8),
      axis.title.y      = element_text(family = "helvetica", face = "plain", size = 8),
      axis.title.x      = element_text(family = "helvetica", face = "plain", size = 8),
      legend.title      = element_blank(),
      legend.text        = element_text(family = "helvetica", face = "plain", size = 8)
    ) +
    labs(title = gene, y = "paternal ratio") 
  
  if(ymax > 0){
    p = p + ylim(c(0, ymax))
  }
  
  return(p)
}

plot_snp_ratios("Rpl38")


