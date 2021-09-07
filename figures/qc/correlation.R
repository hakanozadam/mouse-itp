library(Seurat)
library(ribor)
library(reshape2)
library(edgeR)
library(smatr)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(cowplot)

source('./Ribo_Summary_Function.R')
source('./rename_experiments.R')

human_ribo_file = '../../../../itp/human-itp_v4.ribo'
mouse_ribo_file = '../../../mouse-itp_v5.ribo'
  
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
FIGURE_FONT    = "helvetica"


################################################################################

plot_pairwise_relationships = function (counts_w, 
                                        id1, id2, 
                                        xlab    = "", 
                                        ylab    = "", 
                                        num_bin = 52, 
                                        xrange  = 100000, 
                                        yrange  = 100000  ) { 
  
  sp = ggscatter(counts_w, x = id1, y = id2,
                  #                add = "reg.line", conf.int = FALSE,     
                  #                add.params = list(color = "blue", size = 0.5),
                  font.family = "Helvetica", 
                  size        = 0.2,
                  color       = "gray", 
                  alpha       = 0.3, 
                  ggtheme     = theme_bw()) 
  
  formatted =   sp +   
    scale_x_log10(labels = scales::label_number_si(), limits = c(0.3, xrange)) +   
    scale_y_log10(labels = scales::label_number_si(), limits = c(0.3, yrange)) + 
    labs (x=xlab, y = ylab) +
    stat_cor(method        = "spearman", 
             aes(label     = ..r.label..), 
             cor.coef.name = "rho", 
             digits        = 2)  + 
    geom_hex(bins= num_bin, aes(alpha=log10(..count..) ), fill="#bf5700" ) +
    theme( axis.text.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
           axis.title.x     = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
           axis.text.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
           axis.title.y     = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
           plot.title       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE)
           )
  return (formatted)  
}


################################################################################

human_ribo_rc <- get_region_counts(human,
                             range.lower = MOUSE_MIN_LENGTH,
                             range.upper = MOUSE_MAX_LENGTH,
                             length      = TRUE,
                             transcript  = FALSE,
                             tidy        = FALSE,
                             alias       = TRUE,
                             region      = c("CDS"), 
                             compact     = FALSE)

human_rcw = dcast(human_ribo_rc, transcript ~ experiment)  

################################################################################


# mouse_exp_name_mapper = make_name_mapper( mouse_ribo@experiments  )

#experiment_name_mapper = mouse_exp_name_mapper

human_exp_name_mapper = make_name_mapper( human@experiments  )

human_rename_experiments = function(experiment_name){
  return(human_exp_name_mapper[[experiment_name]])
}

human_rename_experiments = Vectorize(human_rename_experiments, USE.NAMES = FALSE)

colnames(human_rcw) = c( colnames(human_rcw)[1], human_rename_experiments( colnames(human_rcw)[-1]  )  )

human_experiment_names = colnames(human_rcw)[-1]

################################################################################
#### Replicate Clustering By Spearman
####  H E A T M A P

## Replicate clustering by Spearman
## Input is read counts; first column is gene names
replicate_clustering_spearman = function (rcw, 
                                          cpm_threshold = 1, 
                                          breaks_manual = seq(0,1,.01),  
                                          filter        = FALSE, 
                                          clustering    = T){ 
  if (filter) { 
    # Define set of expressed genes before calculating spearman
    expressed = rowSums( cpm (rcw[,-1]) > cpm_threshold) > (dim(rcw)[2] / 3)
    cor_matrix = cor(rcw[expressed,-1], method = "spearman")
  }
  else { 
    cor_matrix = cor(rcw[,-1], method = "spearman")
  }
  
  heatmap_correlation = pheatmap (
                          cor_matrix, 
                          color        =  colorRampPalette(brewer.pal( 9 ,"Greens"))(length(breaks_manual)),
                          #color       =  colorRampPalette( rev(grDevices::heat.colors(100) )  )(length(breaks_manual)),
                          cluster_rows = clustering, 
                          cluster_cols = F,
                          breaks       = breaks_manual,
                          fontsize          = 7,
                          fontsize_col      = 7,
                          fontsize_row      = 7)
  print(heatmap_correlation)
  return(list(cor_matrix = cor_matrix, heatmap_correlation = heatmap_correlation) )
}

#human_cors = replicate_clustering_spearman( human_rcw, breaks_manual = seq(0.6,1,.01), clustering = F )

################################################################################
##### Scatter Plots

## Example Pairwise comparisons. These are not the most beautiful but I don't we should try to make them better
sp_10M_1_vs_10M_2 = plot_pairwise_relationships(
                                 human_rcw, 
                                 human_experiment_names[1], 
                                 human_experiment_names[2], 
                                 xrange  = 2000, 
                                 yrange  = 2000, 
                                 num_bin = 50,
                                 xlab    = human_experiment_names[1], 
                                 ylab    = human_experiment_names[2])
sp_10M_1_vs_10M_2

sp_10M_2_vs_10M_3 = plot_pairwise_relationships(
  human_rcw, 
  human_experiment_names[2], 
  human_experiment_names[3], 
  xrange  = 2000, 
  yrange  = 2000, 
  num_bin = 50,
  xlab    = human_experiment_names[2], 
  ylab    = human_experiment_names[3])
sp_10M_2_vs_10M_3


sp_10M_1_vs_10M_3 = plot_pairwise_relationships(
  human_rcw, 
  human_experiment_names[1], 
  human_experiment_names[3], 
  xrange  = 2000, 
  yrange  = 2000, 
  num_bin = 50,
  xlab    = human_experiment_names[1], 
  ylab    = human_experiment_names[3])
sp_10M_1_vs_10M_3


sp_100_1_vs_100_2 = plot_pairwise_relationships(
  human_rcw, 
  human_experiment_names[4], 
  human_experiment_names[5], 
  xrange  = 2000, 
  yrange  = 2000, 
  num_bin = 50,
  xlab    = human_experiment_names[4], 
  ylab    = human_experiment_names[5])

sp_100_1_vs_100_2


sp_100_2_vs_100_3 = plot_pairwise_relationships(
  human_rcw, 
  human_experiment_names[5], 
  human_experiment_names[6], 
  xrange  = 2000, 
  yrange  = 2000, 
  num_bin = 50,
  xlab    = human_experiment_names[5], 
  ylab    = human_experiment_names[6])

sp_100_2_vs_100_3

sp_100_1_vs_100_3 = plot_pairwise_relationships(
  human_rcw, 
  human_experiment_names[4], 
  human_experiment_names[6], 
  xrange  = 2000, 
  yrange  = 2000, 
  num_bin = 50,
  xlab    = human_experiment_names[4], 
  ylab    = human_experiment_names[6])

sp_100_1_vs_100_3

sp_100_1_vs_100_3 + theme( legend.position = "none"  )


sp_10M_1_vs_100_1 = plot_pairwise_relationships(
  human_rcw, 
  human_experiment_names[1], 
  human_experiment_names[4], 
  xrange  = 2000, 
  yrange  = 2000, 
  num_bin = 50,
  xlab    = human_experiment_names[1], 
  ylab    = human_experiment_names[4])

sp_10M_1_vs_100_1

sp_10M_2_vs_100_1 = plot_pairwise_relationships(
  human_rcw, 
  human_experiment_names[2], 
  human_experiment_names[4], 
  xrange  = 2000, 
  yrange  = 2000, 
  num_bin = 50,
  xlab    = human_experiment_names[2], 
  ylab    = human_experiment_names[4])

sp_10M_2_vs_100_1

sp_10M_3_vs_100_1 = plot_pairwise_relationships(
  human_rcw, 
  human_experiment_names[3], 
  human_experiment_names[4], 
  xrange  = 2000, 
  yrange  = 2000, 
  num_bin = 50,
  xlab    = human_experiment_names[3], 
  ylab    = human_experiment_names[4])

sp_10M_3_vs_100_1

sp_grid_alt = 
plot_grid( sp_10M_1_vs_10M_2 + theme( legend.position = "none"  ),
           sp_10M_2_vs_10M_3 + theme( legend.position = "none"  ),
           sp_10M_1_vs_10M_3 + theme( legend.position = "none"  ),
           sp_100_1_vs_100_2 + theme( legend.position = "none"  ),
           sp_100_2_vs_100_3 + theme( legend.position = "none"  ),
           sp_100_1_vs_100_3 + theme( legend.position = "none"  ),
           sp_10M_1_vs_100_1 + theme( legend.position = "none"  ),
           sp_10M_2_vs_100_1 + theme( legend.position = "none"  ),
           sp_10M_3_vs_100_1 + theme( legend.position = "none"  ),
           align = "hv",
           ncol  = 3)

sp_grid_alt


cpm_threshold = 1
detected      = rowSums( cpm (human_rcw[,-1]) > cpm_threshold) > (dim(human_rcw)[2] / 3)

## Comparing the variance/mean relationship and ITP ~ Monosome similiarity
human_means_cpm = data.frame("Average_100" = apply( cpm(human_rcw[detected,-1])[,4:6]  , 1, mean),
                             "Average_10M" = apply( cpm(human_rcw[detected,-1])[,1:3] , 1, mean)  )

scatter_axis_range = 10000

human_100_vs_10M_means = plot_pairwise_relationships(
  human_means_cpm, 
  "Average_100", 
  "Average_10M", 
  xrange  = scatter_axis_range, 
  yrange  = scatter_axis_range, 
  num_bin = 50,
  xlab    = "100 Average", 
  ylab    = "10M Average")

human_100_vs_10M_means

sp_10M_2_vs_10M_3_new = plot_pairwise_relationships(
  human_rcw, 
  human_experiment_names[2], 
  human_experiment_names[3], 
  xrange  = scatter_axis_range, 
  yrange  = scatter_axis_range, 
  num_bin = 50,
  xlab    = human_experiment_names[2], 
  ylab    = human_experiment_names[3])

sp_100_1_vs_100_3_new = plot_pairwise_relationships(
  human_rcw, 
  human_experiment_names[4], 
  human_experiment_names[6], 
  xrange  = scatter_axis_range, 
  yrange  = scatter_axis_range, 
  num_bin = 50,
  xlab    = human_experiment_names[4], 
  ylab    = human_experiment_names[6])

sp_grid = 
  plot_grid( 
             sp_10M_2_vs_10M_3_new      + theme( legend.position = "none"  ),
             sp_100_1_vs_100_3_new      + theme( legend.position = "none"  ),
             human_100_vs_10M_means     + theme( legend.position = "none"  ),
             align = "hv",
             ncol  = 3)

sp_grid

################################################################################

cpm_threshold = 1
detected      = rowSums( cpm (human_rcw[,-1]) > cpm_threshold) > (dim(human_rcw)[2] / 3)

## Comparing the variance/mean relationship and ITP ~ Monosome similiarity
means_sd_cpm = data.frame(Mean_log2CPM_100 = apply( log2(cpm(human_rcw[detected,-1])[,4:6] + 0.5) , 1, mean),
                          Mean_log2CPM_Monosome = apply( log2(cpm(human_rcw[detected,-1])[,1:3] + 0.5 ), 1, mean),
                          Sd_log2CPM_100 = apply( log2(cpm(human_rcw[detected,-1])[,4:6] + 0.5) , 1, sd),
                          Sd_log2CPM_Monosome = apply( log2(cpm(human_rcw[detected,-1])[,1:3] + 0.5 ), 1, sd)
)
plot(means_sd_cpm$Mean_log2CPM_100, sqrt(means_sd_cpm$Sd_log2CPM_100), pch=19, cex = 0.4)
plot(means_sd_cpm$Mean_log2CPM_Monosome, sqrt(means_sd_cpm$Sd_log2CPM_Monosome), pch=19, cex = 0.4)

average_correlation = cor.test(means_sd_cpm$Mean_log2CPM_100, means_sd_cpm$Mean_log2CPM_Monosome, 
         method = "spearman")
# plot(means_sd_cpm$Mean_log2CPM_100, means_sd_cpm$Mean_log2CPM_Monosome , pch = 19, cex = .5) 

## Spearman Correction -> True correlation Coefficient
# https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005206

# We might want to think of how to estimate reliability better. 
# I am going to use the highest observed correlation coefficient (Change: I pciked 0.94 from mnosome data)

r_true = average_correlation$estimate / 0.94

# We can use the SMA residuals to plot against txn length, GC-content, etc
s1 = sma (Mean_log2CPM_100  ~ Mean_log2CPM_Monosome, 
          data = means_sd_cpm) 
plot(s1, pch = 19, cex = 0.3, las = 1, xlab = "Conventional\n Mean log2 CPM", ylab = "PAC-ITPl\n Mean log2 CPM")
sma_res = residuals(s1)


################################################################################
################################################################################
################################################################################
################################################################################
###                       M O U S E   F I G U R E S                          ###

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

################################################################################
#########         M o u s e     S c a t t e r     P l o t s           ##########

mouse_scatter_axis_range = 1000

sp_1cell_2_5 = plot_pairwise_relationships(
  rcw, 
  "1cell-2", 
  "1cell-5", 
  xrange  = mouse_scatter_axis_range, 
  yrange  = mouse_scatter_axis_range, 
  num_bin = 50,
  xlab    = "1cell-2", 
  ylab    = "1cell-5")

sp_1cell_2_5 + theme( legend.position = "none"  )

sp_mii_3_4 = plot_pairwise_relationships(
  rcw, 
  "MII-3", 
  "MII-4", 
  xrange  = mouse_scatter_axis_range, 
  yrange  = mouse_scatter_axis_range, 
  num_bin = 50,
  xlab    = "MII-3", 
  ylab    = "MII-4")

sp_mii_3_4 + theme( legend.position = "none"  )


sp_gv_3_4 = plot_pairwise_relationships(
  rcw, 
  "GV-3", 
  "GV-4", 
  xrange  = mouse_scatter_axis_range, 
  yrange  = mouse_scatter_axis_range, 
  num_bin = 50,
  xlab    = "GV-3", 
  ylab    = "GV-4")

sp_gv_3_4 + theme( legend.position = "none"  )

mouse_sp_grid = 
  plot_grid( 
    sp_gv_3_4        + theme( legend.position = "none"  ),
    sp_mii_3_4       + theme( legend.position = "none"  ),
    sp_1cell_2_5     + theme( legend.position = "none"  ),
    align = "hv",
    ncol  = 3)

mouse_sp_grid

################################################################################

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
sub_0   = subsampleMatrix( selected_rcw, 0)

# Plot number of genes detected as a function of UMI count
# We should plot a line graph that shows cell to cell relationship as well the sum across cells. 
genes_detected_df = data.frame( 
  sample_depth = c( colSums(selected_rcw ), colSums(sub_40k), colSums(sub_30k), 
                    colSums(sub_20k), colSums(sub_10k), colSums(sub_5k), colSums(sub_0)), 
  genes_detected = c(   apply (selected_rcw, 2, function(x){sum(x > 0)}), 
                        apply (sub_40k, 2, function(x){sum(x > 0)}), 
                        apply (sub_30k, 2, function(x){sum(x > 0)}), 
                        apply (sub_20k, 2, function(x){sum(x > 0)}), 
                        apply (sub_10k, 2, function(x){sum(x > 0)}),
                        apply (sub_5k, 2, function(x){sum(x > 0)}), 
                        apply (sub_0, 2, function(x){sum(x > 0)})
                    ), 
  sample_id =  rep(colnames(selected_rcw), 7)
                  
  )

## Ggplot will be finalized
ggplot(data=genes_detected_df, aes(x=sample_depth, y=genes_detected, group=sample_id) ) +
  geom_line(linetype = "dashed", aes(col = sample_id) )  +  ylim(c( 0, 7000) ) +
  geom_point() + xlab ("Number of CDS Mapping UMIs") + ylab("Genes Detected") + theme_bw()



## Mouse correlations
rcw_columns_reordered = rcw[, c(1, 17:26, 2:16)]
mouse_cors = replicate_clustering_spearman( rcw_columns_reordered, breaks_manual = seq(0,1,.01), filter = T, clustering = F)

## We can add a few example pairwise comparisons: 
plot_pairwise_relationships(rcw, 
                            "GV-3", 
                            "GV-5", 
                            xrange = 3000, yrange = 3000, num_bin = 50,
                            xlab = "GV-3", ylab = "GV-5")


plot_pairwise_relationships(rcw, 
                            "1cell-3", 
                            "1cell-5", 
                            xrange = 1000, yrange = 1000, num_bin = 50,
                            xlab = "1cell-3", ylab = "1cell-5")

# my matrix subsampling gives very similar results
set.seed(3)
test_normalize = LogNormalize(as.matrix(rcw_columns_reordered[,-1]), scale.factor = 30000, verbose = TRUE)
variables = FindVariableFeatures(test_normalize, selection.method = "vst")
rcw_columns_reordered[variables$vst.variance.standardized > 5,] 
## Spin1 is a great one.  -> https://journals.biologists.com/dev/article/124/2/493/39566/Spindlin-a-major-maternal-transcript-expressed-in
#Rfpl4 seems very interesting as well
## AU022751 -> Ovary specific

## Zbed3 protein accumulates by 8-cell but translation is over by 4-cell stage
# https://www.cell.com/cell-reports/fulltext/S2211-1247(17)31795-3
variance_threshold = 4.5
pheatmap(test_normalize[variables$vst.variance.standardized> variance_threshold, ], 
         labels_row = strip_extension(rcw_columns_reordered[variables$vst.variance.standardized > variance_threshold,1]), 
         cluster_cols = F, 
         cutree_rows  = 6 )


################################################################################
## Mouse RNA-Seq
mouse_rna_cds = read.csv('./rnaseq/cds_counts.csv')
colnames(mouse_rna_cds)  = gsub(colnames(mouse_rna_cds), pattern = ".", fixed = T, replacement = "-")

for (id in 1:length(mouse_rna_cds$transcript)) { 
  mouse_rna_cds$transcript[id] = rename_default(mouse_rna_cds$transcript[id])
}
## Remove 20210607.RNAseq.4cell.cross.B
mouse_rna_cds = mouse_rna_cds[,-11]

## RNA-Seq correlation
mouse_rna_cors = replicate_clustering_spearman( mouse_rna_cds, breaks_manual = seq(0,1,.01), filter = T)

## RNA-Seq clustering by variable genes
## Based on the above we might want to change this as well
test_normalize_rna = LogNormalize(mouse_rna_cds[,-1], scale.factor = 250000, verbose = TRUE)
variables = FindVariableFeatures(test_normalize_rna, selection.method = "vst")
mouse_rna_cds[variables$vst.variance.standardized > 7, 1] 

pheatmap(log10 (mouse_rna_cds[variables$vst.variance.standardized > 7,-1] + 1 ) , 
         labels_row = mouse_rna_cds[variables$vst.variance.standardized > 7,1],
)

## Combined data correlation
all_counts = merge(mouse_rna_cds, rcw, by = "transcript")
## Might want to think about which genes are included in the subsampling. 
## It might be good to filter based on expression before subsampling. 
allcounts_cds_subsample = subsampleMatrix (all_counts[,-1], desiredDepth = 28000)

combined_cors = replicate_clustering_spearman( allcounts_cds_subsample, breaks_manual = seq(0,1,.01), filter = T)

plot_pairwise_relationships(allcounts_cds_subsample, "20210614-ITP-GV-50-C", 
                            "20210614-ITP-GV-50-E", 
                            xrange = 2000, yrange = 2000, num_bin = 50,
                            xlab = "GV_RepC", ylab = "GV_RepE")

plot_pairwise_relationships(allcounts_cds_subsample, "20210614-ITP-GV-50-C", 
                            "X20210607-RNAseq-GV-B", 
                            xrange = 2000, yrange = 2000, num_bin = 50,
                            xlab = "MII_Ribo", ylab = "MII_RNASeq")

## A primitive plotting function for single genes. 
## HO: I think this does warrant some effort to make it prettier. 

plot_cpm_across_conditions = function(counts, gene, stages = c("all")) { 
  # counts is a numeric matrix of raw read counts + gene_ids
  # gene to be plotted. Formatted as "Obox2"
  # Stages define which samples to group
  if (stages[1] == "all" )  {  
    cpms = cpm(counts[,-1])
  } else { 
    selected_samples = grep(paste(stages, collapse  = "|"), colnames(counts))
    cpms = cpm(counts[,selected_samples])
  } 
  tidy_cpm = melt(cpms[ strip_extension(as.character(counts$transcript)) %in% gene,  ] )
  tidy_cpm$stage = sapply(strsplit(as.character(colnames(cpms)), split = "-"), "[[", 3) 
  tidy_cpm$method = sapply(strsplit(as.character(colnames(cpms)), split = "-"), "[[", 2) 
  
  ggplot(tidy_cpm, 
         aes( x =  strip_extension(as.character(gene)), y = value) ) + 
    geom_jitter(aes(color = stage, shape = method),
                position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.2),
                size = 1 ) +
    # stat_summary( aes(color = variable),
    #  geom="line", lwd=2, fun = mean, position = position_dodge(0.4) ) + 
    # scale_color_manual(values =  c('#089099', '#d12959')) + 
    theme_bw()+
    ylab("CPM") + xlab("") 
}

## I have a list with a large number of genes that might be useful to include. 
## Here is one example while we work on the aesthetics. 

plot_cpm_across_conditions(all_counts, "Nfrkb", stages = c("RNAseq-1cell", "RNAseq-MII", "ITP-1cell", "ITP-MII", "RNAseq-2cell", "ITP-2cell") ) 
plot_cpm_across_conditions(all_counts, "Ppig", stages = c("RNAseq-MII", "RNAseq-1cell", "RNAseq-2cell", "RNAseq-4cell", "ITP-MII", "ITP-1cell","ITP-2cell" , "ITP-4cell") ) 
plot_cpm_across_conditions(all_counts, "Nop58", stages = c("RNAseq-1cell", "RNAseq-2cell", "RNAseq-4cell", "ITP-1cell","ITP-2cell" , "ITP-4cell") ) ## Proteomics


## Still Need to Proportionality; Differential TE analysis; Proteomics from Can's main script. 

################################################################################
#######                P D F    P R O D U C T I O N                      #######

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

save_plot_pdf("supp_pairwise_correlation.pdf", sp_grid, width = unit(7.05, "in"), height = unit(2.35, "in"))

save_plot_pdf("mouse_correlation_map.pdf", mouse_cors$heatmap_correlation[[4]], width = unit(4.7, "in"), height = unit(4.7, "in"))

save_plot_pdf("supp_mouse_pairwise_correlation.pdf", mouse_sp_grid, width = unit(7.05, "in"), height = unit(2.35, "in"))


