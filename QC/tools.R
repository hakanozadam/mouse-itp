## FROM Can 

summarize_ribo_results = function (ribo_file, min_len = 28, max_len = 35) { 
  library(ribor)
  library(edgeR)
  library(pheatmap)
  p2 = plot_length_distribution(x           = ribo_file,
                           region      = "CDS",
                           range.lower = 15,
                           range.upper = 40,
                           fraction    = TRUE)
  print(p2)
  
  for(seq_len in 15:40) {
    p1 = plot_region_counts(x           = ribo_file,
                            range.lower = seq_len,
                            range.upper = seq_len, 
                            title= seq_len)
    print (p1)
  }
  
  
  p3 = plot_region_counts(x           = ribo_file,
                     range.lower = min_len,
                     range.upper = max_len)
  print(p3)
  
  p4 = plot_metagene(ribo_file,
                site        = "stop",
                normalize   = TRUE,
                title       = "Stop Site Coverage",
                range.lower = min_len,
                range.upper = max_len)
  print(p4)
  
  p5 = plot_metagene(ribo_file,
                site        = "start",
                normalize   = TRUE,
                title       = "Start Site Coverage",
                range.lower = min_len,
                range.upper = max_len)
  print(p5)
  
  ribo_rc <- get_region_counts(ribo_file,
                               range.lower = min_len,
                               range.upper = max_len,
                               length      = TRUE,
                               transcript  = FALSE,
                               tidy = F,
                               alias       = TRUE,
                               region      = c("CDS"), 
                               compact = F)
  
  rcw = dcast(ribo_rc, transcript ~ experiment)  
  expressed_ribo = rowSums ( cpm(rcw[,-1]) > 1) > 2
  c3 = cor(rcw[expressed_ribo,-1], method = "spearman")
  pheatmap(c3, main  = "Ribosome Profiling")
  
  rnaseq <- get_rnaseq(ribo.object = ribo_file,
                       tidy        = F,
                       compact = F, 
                       region = "CDS",
                       alias       = TRUE)
  rnaseq_w = dcast(rnaseq, transcript ~ experiment)
  
}
