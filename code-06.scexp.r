rm(list = ls())
library(Seurat)
library(IOBR)
library(tidydr)

sc.obj.merge = readRDS("interm/annotated.rds")
gs = read.delim2("03.Model/08.Selected.xls")
gs = gs$gene

sc.sub.uc = subset(sc.obj.merge, Group == "UC")
sc.sub.ct = subset(sc.obj.merge, Group == "NC")

cell.cols = rev(pal_d3("category20c")(16))
exp.cols = colorRampPalette(c('#253494','#41b6c4','#f7fcf0','#ffffb2','#fecc5c','#e31a1c'))

for(g in gs){
  p1 = VlnPlot(sc.sub.ct, g, cols = cell.cols) + 
    ggtitle(NULL) +
    NoLegend() + theme(axis.title.x.bottom = element_blank())
  p2 = FeaturePlot(sc.sub.ct, g, cols = exp.cols(50), slot = "data", order = T) + 
    ggtitle(NULL) +
    theme_dr() +
    theme(axis.title = element_blank(),panel.grid = element_blank())
  p3 = VlnPlot(sc.sub.uc, g, cols = cell.cols) + 
    ggtitle(NULL) +
    NoLegend() + theme(axis.title.x.bottom = element_blank())
  p4 = FeaturePlot(sc.sub.uc, g, cols = exp.cols(50), slot = "data", order = T) + 
    ggtitle(NULL) +
    theme_dr() +
    theme(axis.title = element_blank(),panel.grid = element_blank())
  
  pl1 = wrap_plots(p1, p2, widths = c(2,1))
  pl2 = wrap_plots(p3, p4, widths = c(2,1))
  
  text1 = wrap_elements(grid::textGrob('Control', gp = gpar(fontface = "bold")))
  text2 = wrap_elements(grid::textGrob('Disease', gp = gpar(fontface = "bold")))
  
  p = wrap_plots(text1, pl1, text2, pl2, heights = c(1,3,1,3)) + 
    plot_annotation(title = g, 
                    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  fn1 = paste0("05.Expression/", g, ".png")
  fn2 = paste0("05.Expression/", g, ".pdf")
  ggsave(fn1, p, width = 7, height = 7, units = "in", dpi = 300, bg = "white")
  ggsave(fn2, p, width = 7, height = 7, units = "in", dpi = 300, bg = "white")
}
