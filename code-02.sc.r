rm(list = ls())
library(Seurat)
library(ggsci)
library(tidyverse)
library(AUCell)
library(gghalves)
library(ggpubr)

sc.obj.merge = readRDS("interm/integrated.rds")
p1 = VlnPlot(sc.obj.merge, "nFeature_RNA", group.by = "orig.ident") + 
  NoLegend() + 
  scale_fill_d3() +
  xlab(NULL) +
  theme(axis.text.x.bottom = element_blank())
p2 = VlnPlot(sc.obj.merge, "nCount_RNA", group.by = "orig.ident") + 
  NoLegend() + 
  scale_fill_d3() +
  xlab(NULL) +
  theme(axis.text.x.bottom = element_blank())
p3 = VlnPlot(sc.obj.merge, "percent.mt", group.by = "orig.ident") + 
  scale_fill_d3() +
  NoLegend()
p = p1+p2+p3
ggsave("01.SC/01.QC.png", p, width = 9, height = 9, units = "in", dpi = 300, bg = "white")
ggsave("01.SC/01.QC.pdf", p, width = 9, height = 9, units = "in", dpi = 300, bg = "white")


dims = 1:50
sc.obj.merge = FindNeighbors(sc.obj.merge, reduction = "integrated.cca", dims = dims)
sc.obj.merge = FindClusters(sc.obj.merge, resolution = 2)

sc.obj.merge = RunUMAP(sc.obj.merge, dims = 1:50, reduction = "integrated.cca",
                       n.neighbors = 8, min.dist = 0.5)
DimPlot(sc.obj.merge, label = F, repel = T, shuffle = T)
ggsave("01.SC/02.Cluster.png", width = 6.7, height = 6, units = "in", dpi = 300, bg = "white")
ggsave("01.SC/02.Cluster.pdf", width = 6.7, height = 6, units = "in", dpi = 300, bg = "white")

DotPlot(sc.obj.merge, cols = c("grey100","darkred"), scale.max = 80, scale.min = 0,
        scale.by = "size", dot.scale = 6, col.min = 0,
        features = c("TOP2A", "PCNA", #TA
                     "SLC26A2", "FABP1", # colonocyte
                     "TRPM5", "SH2D6", # Tuft
                     "MUC2", "TFF3", # goblet
                     "DCN","COL1A2",
                     "PECAM1","CLDN5", # endothelial
                     "S100B", "CDH19" , #glial
                     "CD3D","CD3E",
                     "NKG7","KLRD1",
                     "FOXP3","IL2RA",
                     "CD79A","MS4A1",
                     "CTSW", "KLRB1", # mait
                     "CD38", "MZB1",
                     "CD68","C1QA", # macro
                     "FCGR3B", "CSF3R", # neutrophil
                     "MS4A2","TPSAB1")) +
  coord_flip()
ggsave("01.SC/03.Cluster.Gene.png", width = 11, height = 8, units = "in", dpi = 300, bg = "white")
ggsave("01.SC/03.Cluster.Gene.pdf", width = 11, height = 8, units = "in", dpi = 300, bg = "white")

clst.type = rep("", length = nlevels(sc.obj.merge$seurat_clusters))
names(clst.type) = seq_len(nlevels(sc.obj.merge$seurat_clusters))-1
clst.type[c("0","1","2","4","13","17","20","24")] = "Plasma"
clst.type[c("5")] = "B"
clst.type[c("3","8","12","28","25")] = "T"
clst.type[c("9")] = "Treg"
clst.type[c("6","7","22")] = "NKT"
clst.type[c("33")] = "MAIT"
clst.type[c("18","23")] = "Macrophage"
clst.type[c("29")] = "Neutrophil"
clst.type[c("27")] = "Mast"
clst.type[c("16","19","26")] = "Fibroblast"
clst.type[c("31")] = "Endothelial"
clst.type[c("14")] = "TA"
clst.type[c("21","32")] = "Goblet"
clst.type[c("10","11","15")] = "Colonocyte"
clst.type[c("30")] = "Tuft"
clst.type[c("34")] = "Glial"
which(clst.type == "")

cell.type = sc.obj.merge@meta.data['seurat_clusters']
cell.type$type = clst.type[cell.type$seurat_clusters]
cells = c("Colonocyte","Tuft","Glial","TA","Goblet",
          "Endothelial","Fibroblast","Macrophage",
          "Neutrophil","Mast",
          "Plasma","B","T","Treg","NKT","MAIT")
cell.type$type = factor(cell.type$type, levels = cells)

sc.obj.merge = AddMetaData(sc.obj.merge, metadata = cell.type$type, 
                           col.name = "cell.type")
Idents(sc.obj.merge) = sc.obj.merge$cell.type

cell.cols = rev(pal_d3("category20c")(16))
DimPlot(sc.obj.merge, cols = cell.cols, reduction = "umap", label = F) + 
  theme(aspect.ratio = 1)
ggsave("01.SC/04.Cell.Type.png", width = 6.5, height = 6, units = "in", dpi = 300, bg = "white")
ggsave("01.SC/04.Cell.Type.pdf", width = 6.5, height = 6, units = "in", dpi = 300, bg = "white")

DotPlot(sc.obj.merge, cols = c("grey100","darkred"), scale.max = 80, scale.min = 0,
        scale.by = "size", dot.scale = 6, col.min = 0,
        features = c("SLC26A2", "FABP1",
                     "TRPM5", "SH2D6",
                     "S100B", "CDH19",
                     "TOP2A", "PCNA", 
                     "MUC2", "TFF3",
                     "PECAM1","CLDN5",
                     "DCN","COL1A2",
                     "CD68","C1QA", 
                     "FCGR3B","CSF3R", 
                     "MS4A2","TPSAB1",
                     "CD38", "MZB1",
                     "CD79A","MS4A1",
                     "CD3D","CD3E",
                     "FOXP3","IL2RA",
                     "NKG7","KLRD1",
                     "CTSW", "KLRB1"
                     )) +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1))
ggsave("01.SC/05.Annot.Gene.png", width = 12, height = 6, units = "in", dpi = 300, bg = "white")
ggsave("01.SC/05.Annot.Gene.pdf", width = 12, height = 6, units = "in", dpi = 300, bg = "white")

saveRDS(sc.obj.merge, "interm/annotated.rds")


rm(list = ls())
sc.obj.merge = readRDS("interm/annotated.rds")
df.meta = sc.obj.merge@meta.data
df.meta = df.meta[c("Sample","cell.type")]
df.meta = group_by(df.meta, Sample, cell.type) |> 
  summarise(cnt = n()) |> 
  mutate(perc = cnt/sum(cnt))

df.grp = sc.obj.merge@meta.data[c("Sample","Group")]
df.grp = df.grp[!duplicated(df.grp),]
df.meta$Group = df.grp$Group[match(df.meta$Sample, df.grp$Sample)]

cell.cols = rev(pal_d3("category20c")(16))
ggplot(df.meta, aes(x = Sample, y = perc)) + 
  geom_bar(aes(fill = cell.type), stat = "identity", position = "fill") + 
  facet_wrap(.~Group, scale = "free") +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = cell.cols, name = "Cell type") + 
  ylab("Proportion") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("01.SC/06.Cell.Proportion.png", width = 8, height = 6, units = "in", dpi = 300, bg = "white")
ggsave("01.SC/06.Cell.Proportion.pdf", width = 8, height = 6, units = "in", dpi = 300, bg = "white")


rm(list = ls())
library(msigdbr)
df.m = msigdbr(category = "H")
df.m = df.m[df.m$gs_name == "HALLMARK_GLYCOLYSIS",]
gs = df.m$gene_symbol |> unique()
write.table(data.frame(gene = gs), "01.SC/08.Glycolysis.Gene.xls", quote = F, sep = "\t", row.names = F)


rm(list = ls())
sc.obj.merge = readRDS("interm/annotated.rds")

set.seed(2025)
cells_rankings = AUCell_buildRankings(sc.obj.merge@assays$RNA$counts)

gs = read.delim2("01.SC/07.Glycolysis.Gene.xls")
gs = gs$gene
cells_AUC = AUCell_calcAUC(gs, cells_rankings)
cells_assignment = AUCell_exploreThresholds(cells_AUC)

thr = cells_assignment$geneSet$aucThr$selected
sc.obj.merge$AUC = as.numeric(getAUC(cells_AUC)[1,])
df.auc = sc.obj.merge@meta.data[c('cell.type','Sample','Group','AUC')]
df.auc$cell.type = as.character(df.auc$cell.type)
df.auc$cell = rownames(df.auc)
write.table(df.auc, "interm/sc.auc.tsv", quote = F, sep = "\t", row.names = F)

df.test = sapply(unique(df.auc$cell.type), \(cell){
  df.sub = df.auc[df.auc$cell.type == cell,]
  wilcox.test(df.sub$AUC~df.sub$Group)$p.value
}) |> data.frame()
df.test = cbind(rownames(df.test), df.test)
colnames(df.test) = c("cell.type","pval")
df.test$padj = p.adjust(df.test$pval, method = "BH")
df.test = df.test[order(df.test$padj, decreasing = F),]
df.test = df.test[df.test$cell.type != "Neutrophil",]
write.table(df.test,"01.SC/08.Glycolysis.Diff.xls", quote = F, sep = "\t", row.names = F)

df.auc = df.auc[df.auc$cell.type!="Neutrophil",]
df.mean = group_by(df.auc, Group, cell.type) |> 
  summarise(mean = mean(AUC))
df.auc$cell.type = factor(df.auc$cell.type, levels = df.test$cell.type)
df.mean$cell.type = factor(df.mean$cell.type, levels = df.test$cell.type)

df.test$lab = ifelse(df.test$padj > 0.05, "n.s", 
                     ifelse(df.test$padj > 0.01, "*",
                            ifelse(df.test$padj > 0.001, "**", "***")))

ggplot(df.auc, aes(x = cell.type, y = AUC, fill = Group)) + 
  geom_half_violin(data = subset(df.auc, Group == "NC"), 
                   side = "l", color = "transparent", 
                   position = position_nudge(x=-0.005)) +
  geom_half_violin(data = subset(df.auc, Group == "UC"), 
                   side = "r", color = "transparent",
                   position = position_nudge(x=0.005)) + 
  geom_point(data = subset(df.mean, Group == "NC"), 
             mapping = aes(y = mean),
             position = position_nudge(x=-0.05), 
             show.legend = F) +
  geom_point(data = subset(df.mean, Group == "UC"), 
             mapping = aes(y = mean),
             position = position_nudge(x=+0.05), 
             show.legend = F) + 
  geom_text(data = df.test,
            inherit.aes = F,
            mapping = aes(x = cell.type, y = 0.16, label = lab)) +
  scale_fill_aaas() +
  scale_y_continuous(limits = c(0,0.17)) +
  theme_bw() +
  xlab("Cell type") + 
  theme(axis.text.x.bottom = element_text(angle = 30, hjust = 1))
ggsave("01.SC/09.Glycolysis.Diff.png", width = 9, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("01.SC/09.Glycolysis.Diff.pdf", width = 9, height = 5, units = "in", dpi = 300, bg = "white")


rm(list = ls())
sc.obj.merge = readRDS("interm/annotated.rds")
df.meta = sc.obj.merge@meta.data
df.meta = df.meta[c("Sample","cell.type")]
df.meta = group_by(df.meta, Sample, cell.type) |> 
  summarise(cnt = n()) |> 
  mutate(perc = cnt/sum(cnt))
df.out = reshape2::dcast(df.meta, cell.type~Sample, value.var = "perc")
df.out[is.na(df.out)] = 0
rownames(df.out) = df.out[[1]]
df.out = df.out[-1]

df.grp = sc.obj.merge@meta.data[c("Sample","Group")]
df.grp = df.grp[!duplicated(df.grp),]
df.out = df.out[df.grp$Sample]

apply(df.out, 1, \(x){
  wilcox.test(x~df.grp$Group)$p.value
})

df.meta$Group = df.grp$Group[match(df.meta$Sample, df.grp$Sample)]
df.meta = df.meta[df.meta$cell.type == "Colonocyte",]

ggplot(df.meta, aes(x = Group, y = perc)) + 
  geom_violin(aes(fill = Group), alpha = 0.5) + 
  geom_boxplot(aes(fill = Group), width = 0.1) + 
  scale_fill_aaas(guide = guide_none()) + 
  stat_compare_means(label.x = 1.5, hjust = 0.5) +
  ylab("Colonocyte proportion") +
  theme_classic2()
ggsave("01.SC/10.Colonocyte.png", width = 4.5, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("01.SC/10.Colonocyte.pdf", width = 4.5, height = 5, units = "in", dpi = 300, bg = "white")


rm(list = ls())
sc.obj.merge = readRDS("interm/annotated.rds")
res.marker = FindMarkers(sc.obj.merge, ident.1 = "Colonocyte", only.pos = T)
res.marker = cbind(gene = rownames(res.marker), res.marker)
res.marker = res.marker[c(1,4,5,3,2,6)]
res.marker = res.marker[res.marker$pct.1 > 0.1, ]
write.table(res.marker, "01.SC/11.Colonocyte.Marker.xls", quote = F, sep = "\t", row.names = F)


rm(list = ls())
library(msigdbr)
library(GSVA)
sc.obj.merge = readRDS("interm/annotated.rds")
sc.obj.merge = subset(sc.obj.merge, cell.type == "Colonocyte")
df.exp = sc.obj.merge@assays$RNA$data |> as.data.frame()

df.m = msigdbr()
df.kegg = subset(df.m, gs_subcat %in% c("CP:KEGG"))[c(3,4)]
gs.kegg = split(df.kegg$gene_symbol, df.kegg$gs_name)

gsvaPar = ssgseaParam(as.matrix(df.exp), gs.kegg)
res.gsva = gsva(gsvaPar)
res.gsva = as.data.frame(res.gsva)
rownames(res.gsva) = str_remove(rownames(res.gsva), "KEGG_")

df.grp = sc.obj.merge@meta.data[c("Sample","Group")]
df.grp$cell = rownames(df.grp)
df.grp = df.grp[!duplicated(df.grp),]
res.gsva = res.gsva[df.grp$cell]

res.gsva = t(res.gsva) |> as.data.frame()
df.median = aggregate(res.gsva, by = list(group = df.grp$Group), FUN = median)
rownames(df.median) = df.median[[1]]
df.median = t(df.median[-1]) |> as.data.frame()

p.diff = apply(res.gsva, 2, \(x){
  wilcox.test(x~df.grp$Group)$p.value
})
df.diff = data.frame(pathway = colnames(res.gsva),
                     pval = p.diff,
                     padj = p.adjust(p.diff, "BH"))
df.diff = merge(df.median, df.diff, by = 0)
df.diff = df.diff[c(4,2,3,5,6)]
df.diff = df.diff[order(df.diff$padj, decreasing = F),]
write.table(df.diff, "01.SC/12.Pathway.Diff.xls", quote = F, sep = "\t", row.names = F)

df.diff = rbind(subset(df.diff, NC > UC)[1:5,],
                subset(df.diff, NC < UC)[1:5,])
df.diff$lab = ifelse(df.diff$padj < 0.001, "***",
                     ifelse(df.diff$padj < 0.01, "**", 
                            ifelse(df.diff$padj < 0.05, "*", "n.s")))

df.plot = res.gsva[df.diff$pathway]
df.plot$cell = rownames(df.plot)
df.plot = reshape2::melt(df.plot, id.var = "cell", variable.name = "pathway")
df.plot$lab = df.diff$lab[match(df.plot$pathway, df.diff$pathway)]
df.plot$group = df.grp$Group[match(df.plot$cell, df.grp$cell)]
df.plot$pathway = factor(df.plot$pathway, df.diff$pathway)

ggplot(df.plot, aes(x = value, y = pathway)) + 
  geom_boxplot(aes(color = group)) + 
  scale_color_aaas(name = "Group") + 
  xlab("Score") + ylab("Pathway") +
  theme_classic2()
ggsave("01.SC/13.Pathway.Diff.png", width = 10, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("01.SC/13.Pathway.Diff.pdf", width = 10, height = 7, units = "in", dpi = 300, bg = "white")


rm(list = ls())
library(CellChat)
sc.obj.merge = readRDS("interm/annotated.rds")
sc.obj.merge = subset(sc.obj.merge, cell.type != "Neutrophil")
sc.obj.merge$cell.type = droplevels(sc.obj.merge$cell.type, exclude = setdiff(levels(sc.obj.merge$cell.type),unique(sc.obj.merge$cell.type)))
sc.obj.merge$samples = sc.obj.merge$Sample
cellchat.high = createCellChat(object = subset(sc.obj.merge, Group == "UC"), 
                               group.by = "cell.type", assay = "RNA")
cellchat.low = createCellChat(object = subset(sc.obj.merge, Group == "NC"), 
                              group.by = "cell.type", assay = "RNA")

CellChatDB = CellChatDB.human
CellChatDB.use = subsetDB(CellChatDB)
cellchat.high@DB = cellchat.low@DB = CellChatDB.use

cellchat.high = subsetData(cellchat.high)
cellchat.high = identifyOverExpressedGenes(cellchat.high)
cellchat.high = identifyOverExpressedInteractions(cellchat.high)
cellchat.high = smoothData(cellchat.high, adj = PPI.human)
cellchat.high = computeCommunProb(cellchat.high, type = "triMean")
cellchat.high = filterCommunication(cellchat.high, min.cells = 10)
cellchat.high = computeCommunProbPathway(cellchat.high)
cellchat.high = aggregateNet(cellchat.high)

cellchat.low = subsetData(cellchat.low)
cellchat.low = identifyOverExpressedGenes(cellchat.low)
cellchat.low = identifyOverExpressedInteractions(cellchat.low)
cellchat.low = smoothData(cellchat.low, adj = PPI.human)
cellchat.low = computeCommunProb(cellchat.low, type = "triMean")
cellchat.low = filterCommunication(cellchat.low, min.cells = 10)
cellchat.low = computeCommunProbPathway(cellchat.low)
cellchat.low = aggregateNet(cellchat.low)

p1 = netVisual_heatmap(cellchat.high, measure = "count")
p2 = netVisual_heatmap(cellchat.high, measure = "weight")
ht_opt$legend_gap = unit(10, "mm")
png("01.SC/14.CellChat.UC.png", width = 9, height = 5, units = "in", bg = "white", res = 300)
draw(p1 + p2, newpage = F, column_title = "UC", 
     column_title_gp = gpar(fontsize = 14),
     auto_adjust = T, heatmap_legend_side = "right")
dev.off()
pdf("01.SC/14.CellChat.UC.pdf", width = 9, height = 5)
draw(p1 + p2, newpage = F, column_title = "UC", 
     column_title_gp = gpar(fontsize = 14),
     auto_adjust = T, heatmap_legend_side = "right")
dev.off()

p1 = netVisual_heatmap(cellchat.low, measure = "count")
p2 = netVisual_heatmap(cellchat.low, measure = "weight")
ht_opt$legend_gap = unit(10, "mm")
png("01.SC/15.CellChat.NC.png", width = 9, height = 5, units = "in", bg = "white", res = 300)
draw(p1 + p2, newpage = F, column_title = "NC", 
     column_title_gp = gpar(fontsize = 14),
     auto_adjust = T, heatmap_legend_side = "right")
dev.off()
pdf("01.SC/15.CellChat.NC.pdf", width = 9, height = 5)
draw(p1 + p2, newpage = F, column_title = "NC", 
     column_title_gp = gpar(fontsize = 14),
     auto_adjust = T, heatmap_legend_side = "right")
dev.off()

cellchat.high = netAnalysis_computeCentrality(cellchat.high)
cellchat.low = netAnalysis_computeCentrality(cellchat.low)
object.list <- list(NC = cellchat.low, UC = cellchat.high)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2)) + 
  scale_fill_aaas()
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), 
                           measure = "weight") + 
  scale_fill_aaas()
gg1 + gg2
ggsave("01.SC/16.Interaction.Diff.png", width = 6, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("01.SC/16.Interaction.Diff.pdf", width = 6, height = 4, units = "in", dpi = 300, bg = "white")


num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link))
gg <- list()
for (i in 1:length(object.list)) {
  object.list[[i]] = netAnalysis_computeCentrality(object.list[[i]])
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], 
                                               title = names(object.list)[i], 
                                               weight.MinMax = weight.MinMax)
  gg[[i]] = gg[[i]] + xlim(0,10) + ylim(0,10)
}
patchwork::wrap_plots(plots = gg) + plot_layout(guides = "collect")
ggsave("01.SC/17.Interaction.Diff.Dir.png", width = 9, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("01.SC/17.Interaction.Diff.Dir.pdf", width = 9, height = 5, units = "in", dpi = 300, bg = "white")


gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

png("01.SC/18.Interaction.Diff.Cell.png", width = 9, height = 5, units = "in", bg = "white", res = 300)
draw(gg1 + gg2, newpage = F, column_title = "UC vs NC", 
     column_title_gp = gpar(fontsize = 14),
     auto_adjust = T, heatmap_legend_side = "right")
dev.off()
pdf("01.SC/18.Interaction.Diff.Cell.pdf", width = 9, height = 5)
draw(gg1 + gg2, newpage = F, column_title = "UC vs NC", 
     column_title_gp = gpar(fontsize = 14),
     auto_adjust = T, heatmap_legend_side = "right")
dev.off()


rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = 1, 
        targets.use = NULL, stacked = T, do.stat = TRUE, 
        color.use = c("NC" = "#3B4992FF","UC" = "#EE0000FF")) + 
  ggtitle("Colonocyte as source") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("01.SC/19.Pathway.Colonocyte.Source.png", width = 5, height = 6, units = "in", dpi = 300, bg = "white")
ggsave("01.SC/19.Pathway.Colonocyte.Source.pdf", width = 5, height = 6, units = "in", dpi = 300, bg = "white")

rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, 
        targets.use = 1, stacked = T, do.stat = TRUE, 
        color.use = c("NC" = "#3B4992FF","UC" = "#EE0000FF")) + 
  ggtitle("Colonocyte as target") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("01.SC/20.Pathway.Colonocyte.Target.png", width = 5, height = 6, units = "in", dpi = 300, bg = "white")
ggsave("01.SC/20.Pathway.Colonocyte.Target.pdf", width = 5, height = 6, units = "in", dpi = 300, bg = "white")

png("01.SC/21.LR.Colonocyte.Source.png", width = 9, height = 8, units = "in", bg = "white", res = 300)
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:15), 
                 comparison = c(1,2), angle.x = 45,
                 color.text = c("NC" = "#3B4992FF","UC" = "#EE0000FF"))
dev.off()
pdf("01.SC/21.LR.Colonocyte.Source.pdf", width = 9, height = 8)
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:15), 
                 comparison = c(1,2), angle.x = 45,
                 color.text = c("NC" = "#3B4992FF","UC" = "#EE0000FF"))
dev.off()


png("01.SC/22.LR.Colonocyte.Target.png", width = 9, height = 8, units = "in", bg = "white", res = 300)
netVisual_bubble(cellchat, sources.use = 1:15, targets.use = 1, 
                 comparison = c(1,2), angle.x = 45,
                 color.text = c("NC" = "#3B4992FF","UC" = "#EE0000FF"))
dev.off()
pdf("01.SC/22.LR.Colonocyte.Target.pdf", width = 9, height = 8)
netVisual_bubble(cellchat, sources.use = 1:15, targets.use = 1, 
                 comparison = c(1,2), angle.x = 45,
                 color.text = c("NC" = "#3B4992FF","UC" = "#EE0000FF"))
dev.off()



