
library(plotly)
library(ribor)
library(Seurat)
library(reshape2)
library(edgeR)
library(EnhancedVolcano)
library(M3C, include.only = "umap")
## We can use DESeq2 to take advantage of their effect size shrinkage. 
library(DESeq2)

source('../qc/Ribo_Summary_Function.R')
source("../qc/rename_experiments.R")

mouse_ribo_file = '../../../mouse-itp_v5.ribo'
mouse_rna_file  = '../../../cds_counts.csv.gz'
proteomics_file = '../../../Protein_MouseEmbryo_PMID_29281840.csv'

mm            = Ribo(mouse_ribo_file, rename = rename_default)
mouse_rna_cds = read.csv(mouse_rna_file)
prots         = read.csv (proteomics_file)

## We'll manually fix the col names
colnames(prots)
raw_prots_names = colnames(prots)
colnames(prots) = c(raw_prots_names[1], "1cell", "2cell", "4cell", "8cell", raw_prots_names[-c(1,2,3,4,5)])
colnames(prots)

min_len = 29
max_len = 35

ribo_rc <- get_region_counts(mm,
                             range.lower = min_len,
                             range.upper = max_len,
                             length      = TRUE,
                             transcript  = FALSE,
                             tidy        = F,
                             alias       = TRUE,
                             region      = c("CDS"), 
                             compact     = F)



mouse_exp_ribo_name_mapper = make_name_mapper( mm@experiments  )

rename_experiments_ribo = function(experiment_name){
  conversion = mouse_exp_ribo_name_mapper[[experiment_name]]
  result = paste( conversion, "ribo", sep = "_" )
  return(result)
}

rename_experiments_ribo = Vectorize(rename_experiments_ribo, USE.NAMES = FALSE)

ribo_rc$experiment = rename_experiments_ribo(ribo_rc$experiment)   

rcw = dcast(ribo_rc, transcript ~ experiment)


colnames(mouse_rna_cds)   = gsub(colnames(mouse_rna_cds), pattern = ".", fixed = T, replacement = "-")
mouse_exp_rna_name_mapper = make_name_mapper(colnames(mouse_rna_cds)[-1] )
rename_experiments_rna    = function(experiment_name){
  conversion = mouse_exp_rna_name_mapper[[experiment_name]]
  result     = result = paste( conversion, "rna", sep = "_" ) 
  return(result)
}
rename_experiments_rna    = Vectorize(rename_experiments_rna)
colnames(mouse_rna_cds)   = c( colnames(mouse_rna_cds)[1], rename_experiments_rna(colnames(mouse_rna_cds)[-1])   )

for (id in 1:length(mouse_rna_cds$transcript)) { 
  mouse_rna_cds$transcript[id] = rename_default(mouse_rna_cds$transcript[id])
}
## Remove 20210607.RNAseq.4cell.cross.B
mouse_rna_cds = mouse_rna_cds[,-11]

# Mouse_RNA_CDS as before and doesn't include "X20210607-RNAseq-4cell-cross-B"
all_counts = merge(mouse_rna_cds, rcw, by = "transcript")

# Convert ribo and rna-seq into density 
cds_len = get_region_lengths(mm, alias = T)[,c(1,4)]
# Density corresponds to Reads per 1kb 

allcounts_density = all_counts
for (col_id in 2:dim(all_counts)[2] ) { 
  allcounts_density[,col_id] = 1000* allcounts_density[,col_id] / 
    cds_len$CDS[match( allcounts_density$transcript, cds_len$transcript )]  
}

# CLR normalize density and average across replicates
set.seed(3)
allcounts_cds_clr_density = NormalizeData(allcounts_density[,-1], scale.factor = 10000, normalization.method= "CLR",verbose = TRUE)
rownames(allcounts_cds_clr_density) = all_counts[,1]

stage = unlist(lapply ( strsplit(colnames(all_counts)[-1], split = "-" ) , "[[", 1 )  ) 
method = unlist(lapply ( strsplit(colnames(all_counts)[-1], split = "_" ) , "[[", 2 )  ) 
groups = as.factor(paste(method, stage,  sep= "_") ) 

mean_clr_density = t(apply(allcounts_cds_clr_density, 1, FUN = function (x){ tapply (as.numeric(x), groups, mean)}) ) 


## Add protein expression
prots_rna_ribo = mean_clr_density[strip_extension(rownames(mean_clr_density) ) %in% prots$Gene.Name, ] 
prots_rna_ribo = cbind( prots_rna_ribo, 
                        prots[match(strip_extension(row.names(prots_rna_ribo)) , prots$Gene.Name), ]  ) 

## Remove annotations/blastocyst from the proteomics data
prots_rna_ribo = prots_rna_ribo[, -c(19:22)]

# Note that 451 ids in proteomics data did not match ids in ribo/rna
dim(prots_rna_ribo)

# Subset to non-zero expression; 29 transcripts with no detectable rna or ribo 
prots_rna_ribo = prots_rna_ribo[apply(prots_rna_ribo[,1:12], 1, sum) > 0, ]
dim(prots_rna_ribo)


## We should be using Spearman's correction before reporting these. 
## Estimation of the reliability of the assays based on replicate correlation: 
rna_reliability  = 0.79
ribo_reliability = 0.71
prot_reliability = 0.8


## Note that there are some negative correlation that will get replaced with zero in the following figure
## If using this version might either change color scheme or label as less than zero
## Will split this into RNA - Prot and Ribo-Prot 
# Note that col 13 is genhe name!
rna_prot = prots_rna_ribo[ , c(13, 11, 12, 7:10, 14:18)]
ribo_prot = prots_rna_ribo[, c(13, 5, 6, 1:4,14:18 )]

rna_prot_cors = cor(rna_prot[,-1], method = "spearman") / sqrt(prot_reliability * rna_reliability)
ribo_prot_cors = cor(ribo_prot[,-1], method = "spearman") / sqrt(prot_reliability * ribo_reliability)

rna_prot_cors = rna_prot_cors[1:6, 7:11]
ribo_prot_cors = ribo_prot_cors[1:6, 7:11]

rna_prot_cors
ribo_prot_cors

prots_rna_ribo_density_ordered = prots_rna_ribo[,c(5,6, 1:4, 14:18, 11, 12, 7:10)]

## UMAP clustering of the rank based data

ranking = apply ( prots_rna_ribo[, -c(13,19)], 2, rank)
umap(ranking, dotsize = 2, axistextsize = 12, legendtextsize = 12, textlabelsize = 4, 
     legendtitle = "", text = c("", "", "4 cell", "8 cell", rep("",12), "Morula"),
     labels = as.factor(c(rep("Ribo", 6) , rep( "RNA", 6), rep ("Prot", 5) )), 
) 

plotMDS(ranking)
cluster::pam(t(ranking), k = 5)


## Sankey Diagram of the Correlations
exp_stages  = c("GV", "MII", "1cell", "2cell", "4cell", "8cell")
prot_stages = c("1cell", "2cell", "4cell", "8cell", "morula")



# Remove 8-cell protein add blastocyst
ribo_orange = rgb(228,88,10 , maxColorValue = 255)
rna_blue   = rgb(55,135,192, maxColorValue = 255)

# All columns
selected_ribo = 1:6
selected_rna = 13:17
selected_protein = c(7:9, 11)

# Split into two graphs. One for GV - 1cell; One drawback is the ratio between two is warped. 
selected_ribo = 1:3
selected_rna = 13:15
selected_protein = c(7:9)

## 
# The other for 4,8 cell to morula
selected_ribo = 5:6
selected_rna = 16:17
selected_protein = 11

selected_cols = c(selected_ribo, selected_protein, selected_rna)

## We can create two link_lists and make two diagram to get better spacing
link_list = list( 
      source = c(), 
      target = c(), 
      value = c()
)
  # We will only calculate ITP-Prot; RNA-Prot correlation with Spearman's corrrection
  # We will update link list with >0 correlations
for (column_id in c(selected_ribo, selected_rna)  ) { 
  for (prot_id in selected_protein) { 
    if (column_id < 7) { 
        spearman = cor(prots_rna_ribo_density_ordered[,column_id], prots_rna_ribo_density_ordered[,prot_id], method = "spearman") / sqrt(prot_reliability * ribo_reliability)
      } else { 
        spearman = cor(prots_rna_ribo_density_ordered[,column_id], prots_rna_ribo_density_ordered[,prot_id], method = "spearman") / sqrt(prot_reliability * rna_reliability)
      }
    if (spearman >  0 ) {
      if (column_id < 7) { 
        link_list[["source"]] = c(link_list[["source"]], column_id-1)
        link_list[["target"]] = c(link_list[["target"]], prot_id-1)
        link_list[["value"]] = c(link_list[["value"]], round(spearman, 2) ) 
      } else { 
        link_list[["source"]] = c(link_list[["source"]], prot_id-1)
        link_list[["target"]] = c(link_list[["target"]], column_id-1)
        link_list[["value"]] = c(link_list[["value"]], round(spearman, 2) ) 
      }
    }
      
  }
}

cor_breaks = seq(0,1, by=.1)
cor_breaks = c(0,0.2, 0.4, 0.5, 0.6, 0.7, 0.8)
# Based on quantile of link_list value
cor_breaks = c(0,0.18, 0.32, 0.46, 0.54, 0.731)

# Based on hand selected 
cor_breaks = c(0,0.35, 0.42, 0.47, 0.52, 0.58, 0.731)


link_list$color = cut(link_list$value, breaks = cor_breaks, 
                    labels = colorRampPalette(
                      colors = c("lightgray", "darkblue"))(length(cor_breaks)-1 ) )


node_list_all = list(
  label = colnames(prots_rna_ribo_density_ordered),
  color = c(rep(ribo_orange, 6), rep("green", 5), rep(rna_blue, 6) ),
  y = c(0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 
        0.3, 0.4, 0.6, 0.8,
        0.1, 0.2, 0.3, 0.4, 0.6, 0.7), 
  x = c(rep(0.1, 6), 
        rep(0.3, 4), 
        rep(0.6, 6)), 
  pad = 40,
  thickness = 20,
  line = list(
    color = "black",
    width = 0.5
  )
)


node_list_48 = list(
  label = colnames(prots_rna_ribo_density_ordered),
  color = c(rep(ribo_orange, 6), rep("green", 5), rep(rna_blue, 6) ),
  y = c(0.1, 0.5, 0.55, 0.1, 0.5), 
  x = c(0.1, 0.1, 
        0.3, 
        0.6, 0.6), 
  pad = 40,
  thickness = 16,
  line = list(
    color = "black",
    width = 0.5
  )
)

node_list_GV_2 = list(
  label = colnames(prots_rna_ribo_density_ordered),
  color = c(rep(ribo_orange, 6), rep("green", 5), rep(rna_blue, 6) ),
  y = c(0.1, 0.34, 0.55, 
        0.3, 0.56, 0.77, 
        0.05, 0.24, 0.45), 
  x = c(0.1, 0.1, 0.1, 
        0.3, 0.3, 0.3,
        0.6, 0.6, 0.6), 
  pad = 40,
  thickness = 16,
  line = list(
    color = "black",
    width = 0.5
  )
)


## Node color added manually
fig <- plot_ly(
  type = "sankey",
  orientation = "v",
  textfont = list( family = "Arial"),

#  node = node_list_all , 
  node = node_list_48 , 
#  node =node_list_GV_2,
  link = link_list
)

fig <- fig %>% layout(
  title = "Spearman Correlation",
  font = list(
    size = 10
  )
)

fig

## Orca requires the commandline-tools separately installed
orca(fig, "sankey.pdf")



strip_last_3 = function(genename) {
  return(substr( genename, 1, nchar(genename)-4) ) 
}

## We will add the pairwise differential expression using our typical method and compare to propr
# Let's split the count table to combinations of stages before going through the normalization process

## EDGER
# This is going to assume that we will give ITP and RNA-Seq from two conditions as input
# calculate_differential_RNA_TE = function (count_table, gene_list = all_counts[,1]) { 
#   
#   stage = unlist(lapply ( strsplit(colnames(count_table), split = "-" ) , "[[", 3 )  ) 
#   method = unlist(lapply ( strsplit(colnames(count_table), split = "-" ) , "[[", 2 )  ) 
#   groups = as.factor(paste(method, stage,  sep= "_") ) 
#   
#   y <- DGEList(counts=count_table, group=groups, genes = gene_list)
#   keep <- filterByExpr(y)
#   
#   y <- y[keep,,keep.lib.sizes=FALSE]
#   ## TMMwsp is a bit better for asymmetric feature sets
#   y <- calcNormFactors(y, method = "TMMwsp")
#   design <- model.matrix(~0+ groups)
#   colnames(design) <- levels(groups)
#   y <- estimateDisp(y,design)
#   
#   plotBCV(y)
#   plotMDS(y, labels = groups)
#   
#   fit <- glmQLFit(y,design)
#   
#   rna_samples = colnames(design)[grep("RNA", colnames(design))]
#   ribo_samples = colnames(design)[grep("ITP", colnames(design))]
#   rna_string <<- paste(rna_samples[1], rna_samples[2], sep = " - ")
# #  print(rna_string)
#   dif1 = paste ("(", paste (ribo_samples[1], rna_samples[1], sep = " - "), ")" , sep = "") 
#   dif2 = paste ("(", paste (ribo_samples[2], rna_samples[2], sep = " - "), ")" , sep = "") 
#   te_string <<- paste(dif1, dif2, sep = " - ")
# #  print(te_string)
# # We are going to assume that the rnas and itps are grouped correctly. 
#   my.contrasts = makeContrasts(
#     TE = te_string , 
#     RNA = rna_string, 
#     levels = design
#   )
#   
#   qlf_TE <- glmQLFTest(fit, contrast=my.contrasts[,1])
#   qlf_RNA <- glmQLFTest(fit, contrast=my.contrasts[,2])
#   
#   print(summary(decideTests(qlf_TE, p.value = 0.05, adjust.method = "fdr")))
#   print(summary(decideTests(qlf_RNA, p.value = 0.05, adjust.method = "fdr")))
#   
#   print (table (decideTests(qlf_TE, p.value = 0.05, adjust.method = "fdr"), decideTests(qlf_RNA, p.value = 0.05, adjust.method = "fdr")  ))
#   
#   return(list(qlf_TE, qlf_RNA))  
# }
# 
# all_counts_reordered = all_counts[, c(16:23 , 2:15, 43:47, 48, 24:42)]
# 
# GV_MII = calculate_differential_RNA_TE(all_counts_reordered[,c(1:8, 23:32 )])
# MII_1cell = calculate_differential_RNA_TE(all_counts_reordered[,c(5:12, 28:37 )])
# cell1_2 = calculate_differential_RNA_TE(all_counts_reordered[,c(9:16, 33:40 )])
# cell2_4 = calculate_differential_RNA_TE(all_counts_reordered[,c(13:18, 38:43 )])
# cell4_8 = calculate_differential_RNA_TE(all_counts_reordered[,c(17:22, 41:47 )])
# 


# plot_cpm_across_conditions(all_counts, "Bub1b", stages = c("RNAseq-2cell", "ITP-2cell", "RNAseq-1cell", "ITP-1cell") ) 
# 

################################################################################
## DIFF Analaysis with DESEQ2

#all_counts_reordered = all_counts[, c(16:23 , 2:15, 43:47, 48, 24:42)]

all_counts_reordered = all_counts[, c(16:23 , 2:15, 39:48, 24:38)]


calculate_interaction_based_TE = function (count_table) { 
  stage = unlist(lapply ( strsplit(colnames(count_table), split = "-" ) , "[[", 1 )  ) 
  method = unlist(lapply ( strsplit(colnames(count_table), split = "_" ) , "[[", 2 )  ) 
  
  dds <- DESeqDataSetFromMatrix(countData = count_table,
                              colData = data.frame(stage = as.factor(stage), method = as.factor(method)),
                              design =  ~ stage + method + stage:method)
  row.names(dds) = all_counts[,1]                       

  dds <- DESeq(dds)
# counts(dds)
  plotDispEsts(dds, CV = T)
  interaction_term <<- tail ( resultsNames(dds) , n =1 ) 
  res <- results(dds, name = interaction_term)
  resOrdered <- res[order(res$pvalue),]

  sum(res$padj < 0.01, na.rm=TRUE)
  summary(results(dds, alpha=0.01))

  resLFC <- lfcShrink(dds, coef= interaction_term)
  resLFC <- resLFC[order(resLFC$pvalue),]
  summary(resLFC, alpha = 0.01)

  return(resLFC)
}

GV_MII_deseq2    = calculate_interaction_based_TE(all_counts_reordered[,c(1:8, 23:32 )])
MII_1cell_deseq2 = calculate_interaction_based_TE(all_counts_reordered[,c(5:12, 28:37 )])
cell1_2_deseq2   = calculate_interaction_based_TE(all_counts_reordered[,c(9:16, 33:40 )])
cell2_4_deseq2   = calculate_interaction_based_TE(all_counts_reordered[,c(13:18, 38:43 )])
cell4_8_deseq2   = calculate_interaction_based_TE(all_counts_reordered[,c(17:22, 41:47 )])


calculate_RNADE = function (count_table) { 
  stage = unlist(lapply ( strsplit(colnames(count_table), split = "-" ) , "[[", 1 )  ) 
  
  dds <- DESeqDataSetFromMatrix(countData = count_table,
                                colData = data.frame(stage = as.factor(stage)),
                                design =  ~ stage)
  row.names(dds) = all_counts[,1]                       
  
  dds <- DESeq(dds)
  # counts(dds)
  plotDispEsts(dds, CV = T)
  
  res <- results(dds)
  resOrdered <- res[order(res$pvalue),]
  
  sum(res$padj < 0.01, na.rm=TRUE)
  summary(results(dds, alpha=0.01))
  
  resLFC <- lfcShrink(dds, res= res, coef = 2)
  resLFC <- resLFC[order(resLFC$pvalue),]
  summary(resLFC, alpha = 0.01)
  
  return(resLFC)
}

GV_MII_RNA = calculate_RNADE(all_counts_reordered[,c(1:8)])
MII_1cell_RNA = calculate_RNADE(all_counts_reordered[,c(5:12 )])
cell1_2_RNA = calculate_RNADE(all_counts_reordered[,c(9:16 )])
cell2_4_RNA = calculate_RNADE(all_counts_reordered[,c(13:18)])
cell4_8_RNA = calculate_RNADE(all_counts_reordered[,c(17:22)])


DESeq_results_RNA = as.data.frame(cell1_2_RNA[, c(2,5)] ) 
DESeq_results_TE = as.data.frame(cell1_2_deseq2[, c(2,5)] ) 

combined_table = merge(DESeq_results_TE, DESeq_results_RNA, by = "row.names", all.x = T)

combined_table$TE_significance = ifelse(combined_table$log2FoldChange.x < 0 & combined_table$padj.x < 0.01, -1,
                                        ifelse(combined_table$log2FoldChange.x > 0 & combined_table$padj.x < 0.01, 1, 0) ) 
combined_table$RNA_significance = ifelse(combined_table$log2FoldChange.y < 0 & combined_table$padj.y < 0.01, -1,
                                        ifelse(combined_table$log2FoldChange.y > 0 & combined_table$padj.y < 0.01, 1, 0) ) 
table(combined_table$TE_significance, combined_table$RNA_significance)

## MII - 1 cell is a very interesting stage as TE changes are not accompanied by RNA changes
## However, there is a distinct set of genes with increased translation efficiency

## There are dramatic changes in RNA expression from GV to MII and 1cell- 2cell - 4cell
## We detect almost no changes in TE between 2-4 cell stage; yet 1cell- 2stage there  a large number genes that are TE down despite no apparent RNA change

## Main focus should be on the MII-1cell + 1cell-2cell TE changes. 
## There is much more change at 1-2 cell than 1-2 TE yet there is a sizable unique fraction in TE
## 2-4 is completely driven by RNA changes
## GV has a lot of RNA change but a lot of TE changes as well

## We can probably change the label to highTE in XX; vs highTE in YY
keyvals = ifelse(DESeq_results_TE$log2FoldChange < -1 & DESeq_results_TE$padj < 0.01, '#6e005f',
                 ifelse(DESeq_results_TE$log2FoldChange > 1 & DESeq_results_TE$padj < 0.01, '#045275', 'grey60') )
names(keyvals)[keyvals == '#6e005f'] <- 'lowTE'
names(keyvals)[keyvals == '#045275'] <- 'highTE'
names(keyvals)[keyvals == 'grey60']  <- 'NS'


gv_mii_volcano_data             = as.data.frame(GV_MII_deseq2[, c(2,5)] )
mii_onecell_volcano_data        = as.data.frame(MII_1cell_deseq2[, c(2,5)] )
onecell_twocell_volcano_data    = as.data.frame(cell1_2_deseq2[, c(2,5)] )
twocell_fourcell_volcano_data   = as.data.frame(cell2_4_deseq2[, c(2,5)] )
fourcell_eightcell_volcano_data = as.data.frame(cell4_8_deseq2[, c(2,5)] )

### Volcano Plot Parameters

axisLabSize_main     = 7
titleLabSize_main    = 9
subtitleLabSize_main = 7
captionLabSize_main  = 7
legendIconSize_main  = 1.0
legendLabSize_main   = 8
encircleSize_main    = 1
cutoffLineWidth_main = 0.3
pointSize_main       = 0.2
vlineWidth_main      = 0.3
colAlpha_main        = 0.4
borderWidth_main     = 0.3
gridlines.minor_main = FALSE
gridlines.major_main = FALSE
legendPosition_main  = "none"


#####  GV to MII

keyvals = ifelse(gv_mii_volcano_data$log2FoldChange < -1 & gv_mii_volcano_data$padj < 0.01, '#6e005f',
                 ifelse(gv_mii_volcano_data$log2FoldChange > 1 & gv_mii_volcano_data$padj < 0.01, '#045275', 'grey60') )

names(keyvals)[keyvals == '#6e005f'] <- 'lowTE'
names(keyvals)[keyvals == '#045275'] <- 'highTE'
names(keyvals)[keyvals == 'grey60']  <- 'NS'

gv_mii_volcano_data_modified = 
  gv_mii_volcano_data %>%
  mutate( log2FoldChange = ifelse(log2FoldChange >= 5, 5, log2FoldChange) ) %>%
  mutate( log2FoldChange = ifelse(log2FoldChange <=  -5, -5, log2FoldChange) ) %>%
  mutate( padj = ifelse( -1 * log10(padj) >=  50, 10**(-25), padj) )

volcano_gv_vs_mii = 
EnhancedVolcano(gv_mii_volcano_data_modified[!is.na(keyvals),], 
                lab = rownames(gv_mii_volcano_data_modified)[!is.na(keyvals)], 
                x               = "log2FoldChange" , 
                y               = "padj", 
                title           = "", 
                subtitle        = "GV vs MII",
                selectLab       = c(""), 
                axisLabSize     = axisLabSize_main,
                titleLabSize    = titleLabSize_main,
                subtitleLabSize = subtitleLabSize_main,
                captionLabSize  = captionLabSize_main,
                legendIconSize  = legendIconSize_main,
                legendLabSize   = legendLabSize_main,
                encircleSize    = 1,
                cutoffLineType  = 'blank',
                vline           = c(0),
                #vlineType       = 'solid',
                cutoffLineWidth = cutoffLineWidth_main,
                pointSize       = pointSize_main,
                vlineWidth      = vlineWidth_main,
                caption         = "",
                legendLabels    = c("NS", "NS", "NS", "1% FDR"),
                colCustom       = keyvals[!is.na(keyvals)], 
                colAlpha        = colAlpha_main,
                gridlines.minor = gridlines.minor_main,
                gridlines.major = gridlines.major_main,
                borderWidth     = borderWidth_main
                #drawConnectors = TRUE,
                #pointSize = c(ifelse(qlf_rna$table$PValue< fdr5_nominalp_RNA, 1, 0.2)),
                #widthConnectors = 0.5,
                #colConnectors = 'black',
                )

volcano_gv_vs_mii_main = 
  volcano_gv_vs_mii + 
    theme(axis.ticks      = element_line(colour = "black", size = 0.3),
          legend.position = "none",
          plot.subtitle   = element_text(hjust = 0.5)) + 
    scale_x_continuous(limits = c(-5,5), 
                       breaks = c(-5,0,5),
                       labels = c("<= -5", "0", ">= 5")) + 
    scale_y_continuous(limits = c(0, 25) , 
                       breaks = c(0, 10, 25),
                       labels = c("0", "10", ">= 25"))

volcano_gv_vs_mii_main

#### MII to 1Cell

mii_onecell_volcano_data_modified = 
  mii_onecell_volcano_data %>%
  mutate( log2FoldChange = ifelse(log2FoldChange >= 5, 5, log2FoldChange) ) %>%
  mutate( log2FoldChange = ifelse(log2FoldChange <=  -5, -5, log2FoldChange) ) %>%
  mutate( padj = ifelse( -1 * log10(padj) >=  10, 10**(-10), padj) )

keyvals = ifelse(mii_onecell_volcano_data$log2FoldChange < -1 & mii_onecell_volcano_data$padj < 0.01, '#6e005f',
                 ifelse(mii_onecell_volcano_data$log2FoldChange > 1 & mii_onecell_volcano_data$padj < 0.01, '#045275', 'grey60') )

names(keyvals)[keyvals == '#6e005f'] <- 'lowTE'
names(keyvals)[keyvals == '#045275'] <- 'highTE'
names(keyvals)[keyvals == 'grey60']  <- 'NS'

volcano_mii_vs_1 = 
  EnhancedVolcano(mii_onecell_volcano_data_modified[!is.na(keyvals),], 
                  lab = rownames(mii_onecell_volcano_data_modified)[!is.na(keyvals)], 
                  x               = "log2FoldChange" , 
                  y               = "padj", 
                  title           = "", 
                  subtitle        = "MII vs 1Cell",
                  selectLab       = c(""), 
                  axisLabSize     = axisLabSize_main,
                  titleLabSize    = titleLabSize_main,
                  subtitleLabSize = subtitleLabSize_main,
                  captionLabSize  = captionLabSize_main,
                  legendIconSize  = legendIconSize_main,
                  legendLabSize   = legendLabSize_main,
                  encircleSize    = 1,
                  cutoffLineType  = 'blank',
                  vline           = c(0),
                  #vlineType       = 'solid',
                  cutoffLineWidth = cutoffLineWidth_main,
                  pointSize       = pointSize_main,
                  vlineWidth      = vlineWidth_main,
                  caption         = "",
                  legendLabels    = c("NS", "NS", "NS", "1% FDR"),
                  colCustom       = keyvals[!is.na(keyvals)], 
                  colAlpha        = colAlpha_main,
                  gridlines.minor = gridlines.minor_main,
                  gridlines.major = gridlines.major_main,
                  borderWidth     = borderWidth_main,
                  #drawConnectors = TRUE,
                  #pointSize = c(ifelse(qlf_rna$table$PValue< fdr5_nominalp_RNA, 1, 0.2)),
                  #widthConnectors = 0.5,
                  #colConnectors = 'black',
  )


volcano_mii_vs_1_main = 
  volcano_mii_vs_1 + 
  theme(axis.ticks      = element_line(colour = "black", size = 0.3),
        legend.position = "none",
        plot.subtitle   = element_text(hjust = 0.5)) + 
  scale_x_continuous(limits = c(-5, 5), 
                     breaks = c(-5, 0, 5),
                     labels = c("<= -5", "0", ">= 5")) + 
  scale_y_continuous(limits = c(0, 10), 
                     breaks = c(0, 5, 10),
                     labels = c("0", "5", ">= 10"))

volcano_mii_vs_1_main


#### 1Cell to 2Cell

onecell_twocell_volcano_data_modified = 
  onecell_twocell_volcano_data %>%
  mutate( log2FoldChange = ifelse(log2FoldChange >= 6, 6, log2FoldChange) ) %>%
  mutate( log2FoldChange = ifelse(log2FoldChange <=  -6, -6, log2FoldChange) ) %>%
  mutate( padj = ifelse( -1 * log10(padj) >=  10, 10**(-10), padj) )

keyvals = ifelse(onecell_twocell_volcano_data$log2FoldChange < -1 & onecell_twocell_volcano_data$padj < 0.01, '#6e005f',
                 ifelse(onecell_twocell_volcano_data$log2FoldChange > 1 & onecell_twocell_volcano_data$padj < 0.01, '#045275', 'grey60') )

names(keyvals)[keyvals == '#6e005f'] <- 'lowTE'
names(keyvals)[keyvals == '#045275'] <- 'highTE'
names(keyvals)[keyvals == 'grey60']  <- 'NS'


volcano_1_vs_2 = 
  EnhancedVolcano(onecell_twocell_volcano_data_modified[!is.na(keyvals),], 
                  lab = rownames(onecell_twocell_volcano_data_modified)[!is.na(keyvals)], 
                  x               = "log2FoldChange" , 
                  y               = "padj", 
                  title           = "", 
                  subtitle        = "1Cell vs 2Cell",
                  selectLab       = c(""), 
                  axisLabSize     = axisLabSize_main,
                  titleLabSize    = titleLabSize_main,
                  subtitleLabSize = subtitleLabSize_main,
                  captionLabSize  = captionLabSize_main,
                  legendIconSize  = legendIconSize_main,
                  legendLabSize   = legendLabSize_main,
                  encircleSize    = 1,
                  cutoffLineType  = 'blank',
                  vline           = c(0),
                  #vlineType       = 'solid',
                  cutoffLineWidth = cutoffLineWidth_main,
                  pointSize       = pointSize_main,
                  vlineWidth      = vlineWidth_main,
                  caption         = "",
                  legendLabels    = c("NS", "NS", "NS", "1% FDR"),
                  colCustom       = keyvals[!is.na(keyvals)], 
                  colAlpha        = colAlpha_main,
                  gridlines.minor = gridlines.minor_main,
                  gridlines.major = gridlines.major_main,
                  borderWidth     = borderWidth_main,
                  #drawConnectors = TRUE,
                  #pointSize = c(ifelse(qlf_rna$table$PValue< fdr5_nominalp_RNA, 1, 0.2)),
                  #widthConnectors = 0.5,
                  #colConnectors = 'black',
  )

volcano_1_vs_2

volcano_1_vs_2_main = 
  volcano_1_vs_2 + 
  theme(axis.ticks      = element_line(colour = "black", size = 0.3),
        legend.position = "none",
        plot.subtitle   = element_text(hjust = 0.5)) + 
  scale_x_continuous(limits = c(-6, 6), 
                     breaks = c(-6, 0, 6),
                     labels = c("<= -6", "0", ">= 6")) + 
  scale_y_continuous(limits = c(0, 10), 
                     breaks = c(0, 5, 10),
                     labels = c("0", "5", ">= 10"))

volcano_1_vs_2_main



## 2 to 4

twocell_fourcell_volcano_data_modified = 
  twocell_fourcell_volcano_data %>%
  mutate( log2FoldChange = ifelse(log2FoldChange >= 5, 5, log2FoldChange) ) %>%
  mutate( log2FoldChange = ifelse(log2FoldChange <=  -5, -5, log2FoldChange) ) %>%
  mutate( padj = ifelse( -1 * log10(padj) >=  10, 10**(-6), padj) )

keyvals = ifelse(twocell_fourcell_volcano_data$log2FoldChange < -1 & twocell_fourcell_volcano_data$padj < 0.01, '#6e005f',
                 ifelse(twocell_fourcell_volcano_data$log2FoldChange > 1 & twocell_fourcell_volcano_data$padj < 0.01, '#045275', 'grey60') )

names(keyvals)[keyvals == '#6e005f'] <- 'lowTE'
names(keyvals)[keyvals == '#045275'] <- 'highTE'
names(keyvals)[keyvals == 'grey60']  <- 'NS'

volcano_2_vs_4 = 
  EnhancedVolcano(twocell_fourcell_volcano_data_modified[!is.na(keyvals),], 
                  lab = rownames(twocell_fourcell_volcano_data_modified)[!is.na(keyvals)], 
                  x               = "log2FoldChange" , 
                  y               = "padj", 
                  title           = "", 
                  subtitle        = "2Cell vs 4Cell",
                  selectLab       = c(""), 
                  axisLabSize     = axisLabSize_main,
                  titleLabSize    = titleLabSize_main,
                  subtitleLabSize = subtitleLabSize_main,
                  captionLabSize  = captionLabSize_main,
                  legendIconSize  = legendIconSize_main,
                  legendLabSize   = legendLabSize_main,
                  encircleSize    = 1,
                  cutoffLineType  = 'blank',
                  vline           = c(0),
                  #vlineType       = 'solid',
                  cutoffLineWidth = cutoffLineWidth_main,
                  pointSize       = pointSize_main,
                  vlineWidth      = vlineWidth_main,
                  caption         = "",
                  legendLabels    = c("NS", "NS", "NS", "1% FDR"),
                  colCustom       = keyvals[!is.na(keyvals)], 
                  colAlpha        = colAlpha_main,
                  gridlines.minor = gridlines.minor_main,
                  gridlines.major = gridlines.major_main,
                  borderWidth     = borderWidth_main,
                  #drawConnectors = TRUE,
                  #pointSize = c(ifelse(qlf_rna$table$PValue< fdr5_nominalp_RNA, 1, 0.2)),
                  #widthConnectors = 0.5,
                  #colConnectors = 'black',
  )

volcano_2_vs_4_main = 
  volcano_2_vs_4 + 
  theme(axis.ticks      = element_line(colour = "black", size = 0.3),
        legend.position = "none",
        plot.subtitle   = element_text(hjust = 0.5)) + 
  scale_x_continuous(limits = c(-5, 5), 
                     breaks = c(-5, 0, 5), 
                     labels = c("<= -5", "0" , ">= 5") ) + 
  scale_y_continuous(limits = c(0, 6), 
                     breaks = c(0, 3, 6),
                     labels = c("0", "3", ">= 6"))

volcano_2_vs_4_main


## 4 to 8



fourcell_eightcell_volcano_data_modified =
  fourcell_eightcell_volcano_data %>%
  mutate( log2FoldChange = ifelse(log2FoldChange >= 2, 2, log2FoldChange) ) %>%
  mutate( log2FoldChange = ifelse(log2FoldChange <=  -2, -2, log2FoldChange) ) %>%
  mutate( padj = ifelse( -1 * log10(padj) >=  10, 10**(-1), padj) )

keyvals = ifelse(fourcell_eightcell_volcano_data_modified$log2FoldChange < -1 & fourcell_eightcell_volcano_data_modified$padj < 0.01, '#6e005f',
                 ifelse(fourcell_eightcell_volcano_data_modified$log2FoldChange > 1 & fourcell_eightcell_volcano_data_modified$padj < 0.01, '#045275', 'grey60') )

names(keyvals)[keyvals == '#6e005f'] <- 'lowTE'
names(keyvals)[keyvals == '#045275'] <- 'highTE'
names(keyvals)[keyvals == 'grey60']  <- 'NS'

volcano_4_vs_8 = 
  EnhancedVolcano(fourcell_eightcell_volcano_data_modified[!is.na(keyvals),], 
                  lab = rownames(fourcell_eightcell_volcano_data_modified)[!is.na(keyvals)], 
                  x               = "log2FoldChange" , 
                  y               = "padj", 
                  title           = "", 
                  subtitle        = "4Cell vs 8Cell",
                  selectLab       = c(""), 
                  axisLabSize     = axisLabSize_main,
                  titleLabSize    = titleLabSize_main,
                  subtitleLabSize = subtitleLabSize_main,
                  captionLabSize  = captionLabSize_main,
                  legendIconSize  = legendIconSize_main,
                  legendLabSize   = legendLabSize_main,
                  encircleSize    = 1,
                  cutoffLineType  = 'blank',
                  vline           = c(0),
                  #vlineType       = 'solid',
                  cutoffLineWidth = cutoffLineWidth_main,
                  pointSize       = pointSize_main,
                  vlineWidth      = vlineWidth_main,
                  caption         = "",
                  legendLabels    = c("NS", "NS", "NS", "1% FDR"),
                  colCustom       = keyvals[!is.na(keyvals)], 
                  colAlpha        = colAlpha_main,
                  gridlines.minor = gridlines.minor_main,
                  gridlines.major = gridlines.major_main,
                  borderWidth     = borderWidth_main,
                  #drawConnectors = TRUE,
                  #pointSize = c(ifelse(qlf_rna$table$PValue< fdr5_nominalp_RNA, 1, 0.2)),
                  #widthConnectors = 0.5,
                  #colConnectors = 'black',
  )

volcano_4_vs_8


volcano_4_vs_8_main = 
  volcano_4_vs_8 + 
  theme(axis.ticks      = element_line(colour = "black", size = 0.3),
        legend.position = "none",
        plot.subtitle   = element_text(hjust = 0.5)) + 
  scale_x_continuous(limits = c(-2, 2), 
                     breaks = c(-2, 0, 2), 
                     labels = c("<= -2", "0" , ">= 2") ) + 
  scale_y_continuous(limits = c(0, 1), 
                     breaks = c(0, 1),
                     labels = c("0", ">= 1"))

volcano_4_vs_8_main

################################################################################

keyvals = ifelse(mii_onecell_volcano_data$log2FoldChange < -1 & mii_onecell_volcano_data$padj < 0.01, '#6e005f',
                 ifelse(mii_onecell_volcano_data$log2FoldChange > 1 & mii_onecell_volcano_data$padj < 0.01, '#045275', 'grey60') )

names(keyvals)[keyvals == '#6e005f'] <- 'lowTE'
names(keyvals)[keyvals == '#045275'] <- 'highTE'
names(keyvals)[keyvals == 'grey60']  <- 'NS'

EnhancedVolcano(mii_onecell_volcano_data[!is.na(keyvals),], 
                lab = rownames(mii_onecell_volcano_data_modified)[!is.na(keyvals)], 
                x               = "log2FoldChange" , 
                y               = "padj", 
                title           = "", 
                subtitle        = "MII vs 1Cell",
                #selectLab       = c(""), 
                axisLabSize     = axisLabSize_main,
                titleLabSize    = titleLabSize_main,
                subtitleLabSize = subtitleLabSize_main,
                captionLabSize  = captionLabSize_main,
                legendIconSize  = legendIconSize_main,
                legendLabSize   = legendLabSize_main,
                encircleSize    = 1,
                cutoffLineType  = 'blank',
                vline           = c(0),
                #vlineType       = 'solid',
                cutoffLineWidth = cutoffLineWidth_main,
                pointSize       = pointSize_main,
                vlineWidth      = vlineWidth_main,
                caption         = "",
                legendLabels    = c("NS", "NS", "NS", "1% FDR"),
                colCustom       = keyvals[!is.na(keyvals)], 
                colAlpha        = colAlpha_main,
                gridlines.minor = gridlines.minor_main,
                gridlines.major = gridlines.major_main,
                borderWidth     = borderWidth_main,
                #drawConnectors = TRUE,
                #pointSize = c(ifelse(qlf_rna$table$PValue< fdr5_nominalp_RNA, 1, 0.2)),
                #widthConnectors = 0.5,
                #colConnectors = 'black',
)


## We should pick some example genes; Group them by RNA expression changes
## We should do GO enrichment of these classes. 

volcano_main_plot = plot_grid(volcano_mii_vs_1_main, 
                              volcano_1_vs_2_main,
                              ncol = 2)
