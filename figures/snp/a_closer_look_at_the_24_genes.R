
library(dplyr)
library(reshape2)
library(ggplot2)

#########################################################################

snp_file = "./selected_snps.csv.gz"

#########################################################################

cluster_1_genes = c(  'Nop14',
                     'Slc13a2'
)

cluster_2_genes = c('Cbx3',
                   'Srpk1',
                   'Umps',
                   'Mysm1',
                   'Ppp2ca',
                   'Bcat1'
)

cluster_3_genes = c('Folr1',
                   'Zfp296',
                   'Nin',
                   'Ddx21',
                   'Eif3d',
                   'Pa2g4',
                   'Mrps9',
                   'Tsen2'
)

cluster_4_genes = c('Cdk1',
                   'Baz1a',
                   'Lclat1',
                   'Ncoa3',
                   'Lyar',
                   'Dyrk3',
                   'Ccnh',
                   'Tmppe' 
)

##############################################################################

snp_table_raw = read.csv(snp_file)

snp_table_summary_2 = 
  snp_table_raw %>% group_by(gene) %>% 
  mutate(cds_total = sum(region == "CDS"),
         utr3_total = sum(region == "UTR3"),
         utr5_total = sum(region == "UTR5")) %>%
  mutate( gene_total = cds_total + utr3_total + utr5_total ) %>%
  distinct(gene_total)

hist(snp_table_summary_2$gene_total)

snp_table_summary = 
snp_table_raw %>% group_by(gene) %>% 
  mutate(cds_total = sum(region == "CDS"),
         utr3_total = sum(region == "UTR3"),
         utr5_total = sum(region == "UTR5")) %>%
  mutate( gene_total = cds_total + utr3_total + utr5_total ) %>%
  summarise( cds_percent  = (cds_total / gene_total) * 100,
          utr5_percent = (utr5_total / gene_total) * 100,
          utr3_percent = (utr3_total / gene_total) * 100) %>%
  distinct()

snp_table_long = melt(snp_table_summary, id.vars = "gene", 
                      measure.vars = c("cds_percent", "utr5_percent", "utr3_percent"),
                      value.name = "percentage")

snp_table_long$variable = factor(snp_table_long$variable,
                                 levels = c("utr5_percent",
                                            "cds_percent",
                                            "utr3_percent"),
                                 labels = c("5' utr",
                                            "cds",
                                            "3' utr"))

percent_barplot = function(df, this_title = ""){
  p<-ggplot(df, aes(x=variable, y=percentage, color=variable)) +
    geom_boxplot() + 
    labs(title = this_title, x= "region", y = "percentage") + 
    theme(legend.position = "none")
  
  return(p)
}

percent_barplot(snp_table_long, this_title = "All Genes") 

cluster_4_percentage = 
snp_table_long %>% filter(gene %in% cluster_4_genes)

percent_barplot(cluster_4_percentage, this_title = "Cluster 4")

cluster_3_percentage = 
  snp_table_long %>% filter(gene %in% cluster_3_genes)

percent_barplot(cluster_3_percentage, this_title = "Cluster 3")

cluster_2_percentage = 
  snp_table_long %>% filter(gene %in% cluster_2_genes)

percent_barplot(cluster_2_percentage, this_title = "Cluster 2")

cluster_1_percentage = 
  snp_table_long %>% filter(gene %in% cluster_1_genes)

percent_barplot(cluster_1_percentage , this_title = "Cluster 1")


#######################################################################

### Stop codon analysis

stop_codons = c("TAG", "TAA", "TGA")

snp_table_raw %>% filter(paternal_codon %in% stop_codons ) 

snp_table_raw %>% filter(maternal_codon %in% stop_codons ) 

# We couldn't find stop codons at the SNP sites.

