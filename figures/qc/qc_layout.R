
#### QC FIGURE LAYOUT

source("./correlation.R")

source('./metagene.R')

source('./region_counts.R')

## Scatter plots from human data
page_top      = sp_grid

bottom_legend = get_legend(ppig_stages_plot)



schematic_label = ggdraw() + 
  draw_label(
    "Schematic",
    fontface   = 'plain',
    fontfamily = 'helvetica',
    size       = FONT_LABEL_SIZE,
    angle      = 45
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )


metagene_and_region_counts = plot_grid(mouse_metagene_plot ,  mouse_region_counts_comperative_plot, nrow = 1 )

figure_layout = plot_grid(page_top, schematic_label, 
                          metagene_and_region_counts, 
                          mouse_sp_grid ,page_bottom, 
                          ncol = 1,
                          rel_heights = c(1,1,2,1,1))

#figure_layout



save_plot_pdf("qc_figure_layout.pdf", figure_layout, width = unit(7.20, "in"), height = unit(9, "in"))
