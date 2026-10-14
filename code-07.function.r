rm(list = ls())
library(tidyverse)
library(ggpubr)
library(ggsci)
library(clusterProfiler)
library(enrichplot)
library(msigdbr)

gsea.plot = function(res.kegg, top.hall, gene){
  gsdata <- do.call(rbind, lapply(top.hall, enrichplot:::gsInfo, object = res.kegg))
  gsdata$Description = factor(gsdata$Description, levels = top.hall)
  p1 = ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + theme_classic(14) + 
    theme(panel.grid.major = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(colour = "grey92"), 
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + 
    scale_x_continuous(expand = c(0, 0)) +
    scale_color_manual(values = c(pal_aaas()(10),"gold4")) +
    #scale_color_aaas() +
    ggtitle(gene) +
    geom_hline(yintercept = 0, color = "black", size = 0.8) +
    geom_line(aes_(y = ~runningScore, color = ~Description), size = 1) +
    theme(legend.position = "right", legend.title = element_blank(), legend.background = element_rect(fill = "transparent")) +
    ylab("Running Enrichment Score") + 
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.line.x = element_blank())
  i = 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1 }
  p2 = ggplot(gsdata, aes_(x = ~x)) + geom_linerange(aes_(ymin = ~ymin, ymax = ~ymax, color = ~Description)) + xlab(NULL) + ylab(NULL) + 
    theme_classic(14) + theme(legend.position = "none", 
                              axis.ticks = element_blank(), 
                              axis.text = element_blank(), 
                              axis.line.x = element_blank()) + 
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    # scale_color_aaas() +
    scale_color_manual(values = c(pal_aaas()(10),"gold4"))
  p = aplot::insert_bottom(p1, p2, height = 0.15)
  return(p)
}

df.m = msigdbr()
df.kegg = subset(df.m, gs_subcat %in% c("CP:KEGG"))[c(3,4)]

df.exp = read.delim2("cleandata/GSE107499.exp.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> as.data.frame()
df.pheno = read.delim2("cleandata/GSE107499.pheno.xls")
df.pheno = df.pheno[df.pheno$group == "Disease",]
df.exp = df.exp[df.pheno$sample]
df.exp = t(df.exp) |> as.data.frame()

gs = read.delim2("03.Model/08.Selected.xls")
gs = gs$gene

for(g in gs){
  df.sub = df.exp[-which(colnames(df.exp)==g)]
  exp.g = df.exp[[g]]
  res.cor = cor(df.sub, exp.g, method = "spearman") |> as.vector()
  names(res.cor) = colnames(df.sub)
  res.cor = res.cor[order(res.cor, decreasing = T)]
  res.cor = na.omit(res.cor)
  
  res.gsea = GSEA(res.cor, TERM2GENE = df.kegg, pvalueCutoff = 0.05, 
                  seed = 1, pAdjustMethod = "BH", eps = 0)
  top.kegg = res.gsea@result
  write.table(top.kegg, paste0("06.Function/01.GSEA/", g, ".xls"), quote = F, sep = "\t", row.names = F)
  top.kegg = subset(top.kegg, p.adjust<0.05)
  top.kegg = top.kegg[order(top.kegg$p.adjust, decreasing = F),]
  top.kegg = rbind(subset(top.kegg, NES > 0)[1:5,],
                   subset(top.kegg, NES < 0)[1:5,])
  top.kegg = top.kegg[order(top.kegg$NES, decreasing = T),]
  top.kegg = na.omit(top.kegg)
  top.kegg = top.kegg$Description
  
  p = gsea.plot(res.gsea, top.kegg, g)
  fn1 = paste0("06.Function/01.GSEA/", g, ".png")
  fn2 = paste0("06.Function/01.GSEA/", g, ".pdf")
  ggsave(fn1, p, width = 12, height = 6, units = "in", dpi = 300, bg = "white")
  ggsave(fn2, p, width = 12, height = 6, units = "in", dpi = 300, bg = "white")
}


rm(list = ls())
df.im = read.delim2("interm/cibersort.GSE107499.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> as.data.frame()
df.grp = read.delim2("cleandata/GSE107499.pheno.xls")
df.im = cbind(sample = rownames(df.im), df.im)
df.im = reshape2::melt(df.im, id.var = "sample", variable.name = "cell")
df.im$group = df.grp$group[match(df.im$sample, df.grp$sample)]
df.im = na.omit(df.im)

cell.cols = palettes("random", show_col = F)[1:22]
ggplot(df.im, aes(x = sample, y = value, fill = cell)) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_manual(values = cell.cols, name = "Cell",
                    guide = guide_legend(ncol = 1)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_grid(cols = vars(group), scale = "free", space = "free") + 
  ylab("Proportion") + xlab(NULL) +
  theme(axis.text.x.bottom = element_blank(),
        axis.ticks.x.bottom = element_blank())
ggsave("06.Function/02.CIBERSORT.Stack.png", width = 12, height = 6, units = "in", dpi = 300, bg = "white")
ggsave("06.Function/02.CIBERSORT.Stack.pdf", width = 12, height = 6, units = "in", dpi = 300, bg = "white")


rm(list = ls())
df.im = read.delim2("interm/cibersort.GSE107499.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> 
  as.data.frame()
df.grp = read.delim2("cleandata/GSE107499.pheno.xls")

df.im = df.im[df.grp$sample,]
df.diff = data.frame(
  cell = colnames(df.im),
  p.val = apply(df.im, 2, function(x){
    r = kruskal.test(x~df.grp$group)
    r$p.value
  }))
df.diff$p.adj = p.adjust(df.diff$p.val, method = "BH")
df.diff = df.diff[order(df.diff$p.adj, decreasing = F),]
write.table(df.diff, "06.Function/03.CIBERSORT.Diff.xls", quote = F, sep = "\t", row.names = F)

df.diff = subset(df.diff, p.adj < 0.05)
df.diff$p = signif(df.diff$p.adj, 3)
df.diff$plab = paste0("p.adj = ", df.diff$p)
df.diff$lab = paste0(df.diff$cell, "\n", df.diff$plab)
df.im = df.im[df.diff$cell]

df.im = aggregate(df.im, by = list(group = df.grp$group), FUN = "mean")
rownames(df.im) = df.im[[1]]
df.im = t(df.im[-1]) |> as.data.frame()
df.im$cell = rownames(df.im)
df.im$diff = df.im$Disease - df.im$Control

df.full = merge(df.im, df.diff, by = "cell")
df.full = df.full[order(df.full$diff, decreasing = F),]
df.full$cell = factor(df.full$cell, levels = df.full$cell)

ggplot(df.full, aes(y = cell, x = diff)) + 
  geom_vline(xintercept = 0) +
  geom_segment(aes(xend = 0)) +
  geom_point(aes(size = p.adj, color = diff)) + 
  scale_size_continuous(range = c(7,2), breaks = c(1e-9, 1e-7, 1e-4,1e-2), transform = "log10") +
  scale_color_gradient2(low = "darkblue", high = "darkred", midpoint = 0, 
                        limits = c(-0.12, 0.12), guide = guide_none()) +
  geom_text(data = subset(df.full, diff > 0), 
            mapping = aes(x = 0, label = cell), hjust = 1, 
            position = position_nudge(x = -0.001), 
            color = "darkred") + 
  geom_text(data = subset(df.full, diff < 0), 
            mapping = aes(x = 0, label = cell), hjust = 0, 
            position = position_nudge(x = 0.001),
            color = "darkblue") + 
  scale_x_continuous(limits = c(-0.12, 0.12)) +
  xlab("Mean difference between UC and NC samples") + 
  ylab(NULL) +
  theme_minimal() + 
  theme(axis.text.y = element_blank(),
        panel.background = element_rect(colour = "black", linewidth = 1))
ggsave("06.Function/04.CIBERSORT.Diff.png", width = 8, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("06.Function/04.CIBERSORT.Diff.pdf", width = 8, height = 5, units = "in", dpi = 300, bg = "white")


rm(list = ls())
df.im = read.delim2("interm/cibersort.GSE107499.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> 
  as.data.frame()
df.diff = read.delim2("06.Function/03.CIBERSORT.Diff.xls")
df.diff$p.adj = as.numeric(df.diff$p.adj)
df.diff = na.omit(df.diff)
df.diff = df.diff[df.diff$p.adj<0.05,]
df.im = df.im[df.diff$cell]
df.exp = read.delim2("cleandata/GSE107499.exp.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> as.data.frame()
gs = read.delim2("03.Model/08.Selected.xls")
df.exp = df.exp[gs$gene,rownames(df.im)] |> t() |> as.data.frame()
res.cor = linkET::correlate(df.exp, df.im)
df.cor = res.cor$r |> as.data.frame()
df.p = res.cor$p |> as.data.frame()
df.cor$gene = rownames(df.cor)
df.p$gene = rownames(df.p)
df.cor = reshape2::melt(df.cor, id.var = "gene", variable.name = "cell", value.name = "r")
df.p = reshape2::melt(df.p, id.var = "gene", variable.name = "cell", value.name = "p")
df.full = merge(df.cor, df.p)
write.table(df.full, "06.Function/05.Gene.Cell.Cor.xls", quote = F, sep = "\t", row.names = F)

df.full = na.omit(df.full)
df.full = df.full[df.full$p < 0.05,]
ggplot(df.full, aes(x = gene, y = cell)) + 
  geom_point(aes(size = p, color = r)) + 
  geom_text(aes(label = signif(r, 2))) +
  scale_color_gradient2(low = "blue", high = "red", limits = c(-0.7,0.7)) +
  scale_size_continuous(range = c(7,3)) +
  theme_minimal() + 
  xlab("Gene") + ylab("Cell") +
  theme(axis.text.x = element_text(hjust = 1, angle = 30))
ggsave("06.Function/06.Gene.Cell.Cor.png", width = 5, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("06.Function/06.Gene.Cell.Cor.pdf", width = 5, height = 7, units = "in", dpi = 300, bg = "white")


rm(list = ls())
df.score = read.delim2("interm/GSE107499.glycolysis.xls")
df.score$Glycolysis = as.numeric(df.score$Glycolysis)

df.im = read.delim2("interm/cibersort.GSE107499.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> 
  as.data.frame()
df.diff = read.delim2("06.Function/03.CIBERSORT.Diff.xls")
df.diff$p.adj = as.numeric(df.diff$p.adj)
df.diff = na.omit(df.diff)
df.diff = df.diff[df.diff$p.adj<0.05,]
df.im = df.im[df.diff$cell]
df.im = df.im[df.score$sample,]

df.exp = read.delim2("cleandata/GSE107499.exp.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> as.data.frame()
gs = read.delim2("03.Model/08.Selected.xls")
df.exp = df.exp[gs$gene,df.score$sample] |> t() |> as.data.frame()
df.im = cbind(df.exp, df.im)

res.cor = linkET::correlate(df.im, df.score['Glycolysis'])
df.cor = res.cor$r |> as.data.frame()
df.p = res.cor$p |> as.data.frame()
df.full = data.frame(feature = rownames(df.cor), r = df.cor$Glycolysis, p = df.p$Glycolysis)
write.table(df.full, "06.Function/07.Cell.Glycolysis.Cor.xls", quote = F, sep = "\t", row.names = F)

df.full = df.full[df.full$p < 0.05,]
df.full = df.full[order(df.full$r, decreasing = F),]
df.full$feature = factor(df.full$feature, levels = df.full$feature)
df.full$type = ifelse(df.full$feature %in% gs$gene, "gene", "cell")

ggplot(df.full, aes(y = feature, x = r)) + 
  geom_vline(xintercept = 0) +
  geom_segment(aes(xend = 0)) + 
  geom_point(aes(size = p, color = r, shape = type), alpha = 1) + 
  scale_size_continuous(range = c(6,3), breaks = c(1e-5, 1e-4, 1e-3), transform = "log10") + 
  scale_color_gradient2(low = "blue", high = "red", mid = 0) + 
  geom_text(aes(label = feature, x = 0), hjust = 1, 
            data = subset(df.full, r > 0), color = "red",
            position = position_nudge(x = -0.01)) +
  geom_text(aes(label = feature, x = 0), hjust = 0, 
            data = subset(df.full, r < 0), color = "blue",
            position = position_nudge(x = +0.01)) +
  geom_text(aes(label = signif(r,2)), hjust = 1.4, data = subset(df.full, r < 0)) +
  geom_text(aes(label = signif(r,2)), hjust = -0.4, data = subset(df.full, r > 0)) +
  scale_x_continuous(limits = c(-0.55, 0.55)) +
  xlab("Correlation with glycolysis") +
  ylab("Feature") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        panel.background = element_rect(colour = "black", linewidth = 1))
ggsave("06.Function/08.Gene.Cell.Glycolysis.Cor.png", width = 8, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("06.Function/08.Gene.Cell.Glycolysis.Cor.pdf", width = 8, height = 5, units = "in", dpi = 300, bg = "white")



