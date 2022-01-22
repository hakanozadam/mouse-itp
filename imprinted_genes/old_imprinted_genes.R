
library(cowplot)
library(ggplot2)
library(dplyr)
library(Cairo)

setwd("/data/projects/mice_cell_development/mouse-itp/imprinted_genes/")

## Supp. Data 1 from
# https://www.nature.com/articles/s41467-021-23510-4#MOESM4

CONSTANT_ADDITION = 1
# Note that in the data frames, x refers to ribosome profioling and y refers to RNA-Seq data

published_imprinted_gene_list_file = "./imprinted_gene_list_41467_2021_23510_MOESM4_ESM.csv"


# Another imprinted gene list from
# https://www.cell.com/cell/fulltext/S0092-8674(19)30106-0?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867419301060%3Fshowall%3Dtrue#%20
tucci_et_al_imprinted_gene_list_file = "./tucci_et_al_2019_mouse_imprinted_genes.csv"

  
four_cell_snp_counts_file  = "./fourcell_snp_counts.csv.gz"
eight_cell_snp_counts_file = "./eightcell_snp_counts.csv.gz" 

four_cell_all_snp_counts_file  = "./fourcell_all_snp_counts.csv.gz"
eight_cell_all_snp_counts_file = "./eightcell_all_snp_counts.csv.gz" 

deseq_2_to_4_differential_file = "../figures/proteomics/deseq2_results/cell2_4_deseq2.csv"
deseq_4_to_8_differential_file = "../figures/proteomics/deseq2_results/cell4_8_deseq2.csv"



################################################################################

differential_table_2_to_4 = read.csv(deseq_2_to_4_differential_file)
differential_table_4_to_8 = read.csv(deseq_4_to_8_differential_file)

differential_table_2_to_4$transcript =  unlist(lapply( strsplit(differential_table_2_to_4$transcript, split = "-"), "[[", 1) )
differential_table_4_to_8$transcript =  unlist(lapply( strsplit(differential_table_4_to_8$transcript, split = "-"), "[[", 1) )



################################################################################


published_imprinted_gene_table = read.csv(published_imprinted_gene_list_file, header = 1)

published_imprinted_genes = published_imprinted_gene_table$Gene.symbol

tucci_et_al_imprinted_gene_list = read.csv(tucci_et_al_imprinted_gene_list_file, header = 1)

four_cell_counts  = read.csv(four_cell_snp_counts_file)
eight_cell_counts = read.csv(eight_cell_snp_counts_file)


four_cell_all_snp_counts  = read.csv(four_cell_all_snp_counts_file)
eight_cell_all_snp_counts = read.csv(eight_cell_all_snp_counts_file)


compute_ratios = function(count_df){
  result_df = 
    count_df %>% mutate(
      ribo_ratio = (paternal_x + CONSTANT_ADDITION) / (maternal_x + CONSTANT_ADDITION),
      rna_ratio  = (paternal_y + CONSTANT_ADDITION) / (maternal_y + CONSTANT_ADDITION),
      ratio_diff = ribo_ratio - rna_ratio
    )
    
  return(result_df)
}

four_cell_ratios_df  = compute_ratios(four_cell_counts)
eight_cell_ratios_df = compute_ratios(eight_cell_counts) 

four_cell_all_snp_ratios_df = compute_ratios( four_cell_all_snp_counts )
eight_cell_all_snp_ratios_df = compute_ratios( eight_cell_all_snp_counts )
  
hist( log2(four_cell_ratios_df$ribo_ratio), breaks = 30, 
      main = "Four Cell Ribo log2(paternal / Maternal)") 

hist( log2(four_cell_ratios_df$rna_ratio), breaks = 30, 
      main = "Four Cell RNA log2(paternal / Maternal)")

hist( log2(four_cell_all_snp_ratios_df$rna_ratio), breaks = 30,
      main = "ALL Four Cell RNA log2(paternal / Maternal)")

hist( log2(eight_cell_all_snp_ratios_df$rna_ratio), breaks = 30,
      main = "ALL Eight Cell RNA log2(paternal / Maternal)")

hist( log2(eight_cell_ratios_df$ribo_ratio), breaks = 30, main = "Eight Cell Ribo log2(paternal / Maternal)") 
hist( log2(eight_cell_ratios_df$rna_ratio), breaks = 30, main = "Eight Cell RNA log2(paternal / Maternal)")


fourcell_intersect_published  = four_cell_ratios_df[ four_cell_ratios_df$gene %in% published_imprinted_genes, ]
eightcell_intersect_published = eight_cell_ratios_df[ eight_cell_ratios_df$gene %in% published_imprinted_genes, ]

hist( log2(fourcell_intersect_published$ribo_ratio),  breaks = 30, main = "Four Cell Ribo log2(paternal / Maternal)") 

hist( log2(eightcell_intersect_published$ribo_ratio), breaks = 30, main = "Eight Cell Ribo log2(paternal / Maternal)") 


################################################################################

################################################################################
### Comparison of imprinted genes with published data

View(tucci_et_al_imprinted_gene_list)

tucci_et_al_imprinted_gene_list_mp = 
tucci_et_al_imprinted_gene_list[
  tucci_et_al_imprinted_gene_list$ExpressedAllele %in% c("M", "P"), ]



fourcells_restricted_to_tucci = 
  merge(four_cell_all_snp_ratios_df, tucci_et_al_imprinted_gene_list_mp, 
        by.x = "gene",
        by.y = "Gene")

eightcells_restricted_to_tucci = 
  merge(eight_cell_all_snp_ratios_df, tucci_et_al_imprinted_gene_list_mp, 
        by.x = "gene",
        by.y = "Gene")


setdiff( tucci_et_al_imprinted_gene_list_mp$Gene, 
         fourcells_restricted_to_tucci$gene)


eight_cell_agree = 
(eightcells_restricted_to_tucci$maternal_y + eightcells_restricted_to_tucci$paternal_y ) > 4 &
  (
    (eightcells_restricted_to_tucci$rna_ratio > 1 & eightcells_restricted_to_tucci$ExpressedAllele =="P") |
      (eightcells_restricted_to_tucci$rna_ratio < 1 & eightcells_restricted_to_tucci$ExpressedAllele =="M")
    
  )

sum(eight_cell_agree)

sum((eightcells_restricted_to_tucci$maternal_y + eightcells_restricted_to_tucci$paternal_y ) > 4 )

### / Comparison of imprinted genes with published data
################################################################################


quants_eightcell_rna = quantile(eight_cell_ratios_df$rna_ratio, probs = c(0.05, 0.95))
quants_fourcell_rna = quantile(four_cell_ratios_df$rna_ratio, probs = c(0.05, 0.95))

eight_cell_ribo_maternal_threshold = quants_eightcell_rna[1]
eight_cell_ribo_paternal_threshold = quants_eightcell_rna[2]

eight_cell_rna_maternal_index = eight_cell_ratios_df$rna_ratio <= eight_cell_ribo_maternal_threshold 
eight_cell_rna_paternal_index = eight_cell_ratios_df$rna_ratio >= eight_cell_ribo_paternal_threshold 
eight_cell_rna_background_index = !(eight_cell_rna_maternal_index | eight_cell_rna_paternal_index)


four_cell_ribo_maternal_threshold = quants_fourcell_rna[1]
four_cell_ribo_paternal_threshold = quants_fourcell_rna[2]

four_cell_rna_maternal_index = four_cell_ratios_df$rna_ratio <= four_cell_ribo_maternal_threshold
four_cell_rna_paternal_index = four_cell_ratios_df$rna_ratio >= four_cell_ribo_paternal_threshold
four_cell_rna_background_index = !( four_cell_rna_maternal_index | four_cell_rna_paternal_index )


get_fold_changes = function(diff_df, ratio_df, idx){
  this_index = diff_df$transcript %in% (ratio_df[idx, ])$gene
  return( (diff_df[this_index,])$TE_log2FoldChange )
}


### 4cell RNA Imprinted Genes

# Gene Labels
four_cell_labels_rna = c(
  rep("BG", sum(four_cell_rna_background_index)),
  rep("paternal", sum(four_cell_rna_paternal_index)),
  rep("maternal", sum(four_cell_rna_maternal_index))
)

# Fold Changes from 4 to 2 
TE_foldchanges_2_4 = c(
  get_fold_changes(differential_table_2_to_4, four_cell_ratios_df ,four_cell_rna_background_index),
  get_fold_changes(differential_table_2_to_4, four_cell_ratios_df ,four_cell_rna_paternal_index),
  get_fold_changes(differential_table_2_to_4, four_cell_ratios_df ,four_cell_rna_maternal_index)
)

# Fold Changes from 8 to 4
TE_foldchanges_4_8 = c(
  get_fold_changes(differential_table_4_to_8, four_cell_ratios_df ,four_cell_rna_background_index),
  get_fold_changes(differential_table_4_to_8, four_cell_ratios_df ,four_cell_rna_paternal_index),
  get_fold_changes(differential_table_4_to_8, four_cell_ratios_df ,four_cell_rna_maternal_index)
)




baxplot_df_2_to_4_rna_TE = data.frame( type = four_cell_labels_rna, fold_change =  TE_foldchanges_2_4   )

boxplot_2_4_from_4_imprinted_rna = 
ggplot(baxplot_df_2_to_4_rna_TE, aes(x=type, y=fold_change)) + 
  geom_boxplot() + 
  ggtitle("TE From 4 to 2, RNA imprinted genes in 4cell")

ggsave("boxplot_2_4_from_4_imprinted_rna.pdf", 
       plot   = boxplot_2_4_from_4_imprinted_rna, 
       device = cairo_pdf, 
       width  = 5,
       height = 6,
       dpi    = 600 )

## We test if paternal and maternal TEs are different
t.test( get_fold_changes(differential_table_2_to_4, four_cell_ratios_df ,four_cell_rna_paternal_index),
        get_fold_changes(differential_table_2_to_4, four_cell_ratios_df ,four_cell_rna_maternal_index),
        alternative = "g")

baxplot_df_4_to_8_rna_TE = data.frame( type = four_cell_labels, fold_change =  TE_foldchanges_4_8   )

ggplot(baxplot_df_4_to_8_rna_TE, aes(x=type, y=fold_change)) + 
  geom_boxplot() + 
  ggtitle("TE From 8 to 4, RNA imprinted genes in 4cell")

t.test( get_fold_changes(differential_table_4_to_8, four_cell_ratios_df ,four_cell_rna_paternal_index),
        get_fold_changes(differential_table_4_to_8, four_cell_ratios_df ,four_cell_rna_maternal_index),
        alternative = "l")

# End of 4cell imprinted genes
###############################

### 8cell RNA Imprinted Genes

# Gene Labels
eight_cell_labels_rna = c(
  rep("BG", sum(eight_cell_rna_background_index)),
  rep("paternal", sum(eight_cell_rna_paternal_index)),
  rep("maternal", sum(eight_cell_rna_maternal_index))
)

# Fold Changes from 8 to 4 
TE_foldchanges_4_8_8cell_imprinted_genes = c(
  get_fold_changes(differential_table_4_to_8, eight_cell_ratios_df ,eight_cell_rna_background_index),
  get_fold_changes(differential_table_4_to_8, eight_cell_ratios_df ,eight_cell_rna_paternal_index),
  get_fold_changes(differential_table_4_to_8, eight_cell_ratios_df ,eight_cell_rna_maternal_index)
)


baxplot_df_4_to_8_rna_TE_8cell_imprinted_genes = 
  data.frame( type = eight_cell_labels_rna, fold_change =  TE_foldchanges_4_8_8cell_imprinted_genes   )


ggplot(baxplot_df_4_to_8_rna_TE_8cell_imprinted_genes, aes(x=type, y=fold_change)) + 
  geom_boxplot() + 
  ggtitle("TE From 8 to 4, RNA imprinted genes in 8cell")

t.test( get_fold_changes(differential_table_4_to_8, eight_cell_ratios_df ,eight_cell_rna_paternal_index),
        get_fold_changes(differential_table_4_to_8, eight_cell_ratios_df ,eight_cell_rna_maternal_index),
        alternative = "g")

