rm(list = ls())
library(tidyverse)
library(ggsci)
library(ConsensusClusterPlus)
library(IOBR)
library(GSVA)
library(msigdbr)

df.exp = read.delim2("cleandata/GSE107499.exp.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> as.data.frame()
gs = read.delim2("01.SC/07.Glycolysis.Gene.xls")
gs = gs$gene
gsvaPar = ssgseaParam(as.matrix(df.exp), list(Glycolysis = gs))
res.gsva = gsva(gsvaPar)
res.gsva = as.data.frame(res.gsva)
res.gsva = t(res.gsva) |> as.data.frame()
res.gsva$sample = rownames(res.gsva)
write.table(res.gsva, "interm/GSE107499.glycolysis.xls", quote = F, sep = "\t", row.names = F)


rm(list = ls())
df.exp = read.delim2("cleandata/GSE107499.exp.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> as.data.frame()
df.pheno = read.delim2("cleandata/GSE107499.pheno.xls")
df.pheno = df.pheno[df.pheno$group == "Disease",]
gs = read.delim2("03.Model/08.Selected.xls")
gs = gs$gene
df.exp = df.exp[gs,df.pheno$geo_accession]

res = ConsensusClusterPlus(as.matrix(df.exp), maxK = 10, 
                           clusterAlg = "pam", distance = "pearson",
                           seed = 2025, title = "cluster", plot = "pdf")

grp = res[[4]]$consensusClass
df.grp = data.frame(sample = names(grp), cluster = grp)
df.grp$cluster = paste0("Cluster", df.grp$cluster)
write.table(df.grp, "interm/consensus.cluster.tsv", quote = F, sep = "\t", row.names = F)


rm(list = ls())
df.score = read.delim2("interm/GSE107499.glycolysis.xls")
df.grp = read.delim2("interm/consensus.cluster.tsv")
df.full = merge(df.grp, df.score)
df.full$Glycolysis = as.numeric(df.full$Glycolysis)

ggplot(df.full, aes(x = cluster, y = Glycolysis)) + 
  geom_violin(aes(fill = cluster), alpha = 0.5) + 
  geom_boxplot(aes(fill = cluster), width = 0.2) + 
  stat_compare_means(label.x = 2.5, hjust = 0.5) + 
  scale_fill_d3(guide = guide_none()) +
  ylab('Glycolysis score') +
  xlab("Cluster") +
  theme_classic2()
ggsave("04.Cluster/03.Glycolysis.png", width = 5, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("04.Cluster/03.Glycolysis.pdf", width = 5, height = 5, units = "in", dpi = 300, bg = "white")


rm(list = ls())
df.exp = read.delim2("cleandata/GSE107499.exp.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> as.data.frame()
df.grp = read.delim2("interm/consensus.cluster.tsv")
gs = read.delim2("01.SC/07.Glycolysis.Gene.xls")
df.exp = df.exp[gs$gene, df.grp$sample]
df.exp = na.omit(df.exp)

p.diff = apply(df.exp, 1, \(x){
  kruskal.test(x~df.grp$cluster)$p.value
})

df.diff = data.frame(gene = rownames(df.exp),
                     pval = p.diff,
                     padj = p.adjust(p.diff, "BH"))
df.diff = df.diff[order(df.diff$padj, decreasing = F),]
df.diff = df.diff[df.diff$padj<0.05,]
write.table(df.diff, "04.Cluster/04.Glycolysis.Gene.xls", quote = F, sep = "\t", row.names = F)

df.exp = df.exp[df.diff$gene,]
df.diff$lab = ifelse(df.diff$padj > 0.01, "*", 
                     ifelse(df.diff$padj > 0.001, "**", "***"))
rownames(df.exp) = paste0(df.diff$gene, " (", df.diff$lab, ")")
df.grp = df.grp[order(df.grp$cluster, decreasing = F),]
df.exp = df.exp[df.grp$sample]
rownames(df.grp) = df.grp$sample
gaps = cumsum(table(df.grp$cluster))
gaps = gaps[-length(gaps)]

cluster.cols = pal_d3()(4)
names(cluster.cols) = paste0("Cluster", 1:4)

png("04.Cluster/05.Glycolysis.Gene.png", width = 9, height = 9, units = "in", bg = "white", res = 300)
pheatmap::pheatmap(df.exp, scale = "row", cluster_cols = F, 
                   color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100),
                   show_colnames = F, border_color = "white",
                   gaps_col = gaps,
                   annotation_colors = list(cluster = cluster.cols),
                   annotation_col = df.grp['cluster'])
dev.off()
pdf("04.Cluster/05.Glycolysis.Gene.pdf", width = 9, height = 9)
pheatmap::pheatmap(df.exp, scale = "row", cluster_cols = F, 
                   color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100),
                   show_colnames = F, border_color = "white",
                   gaps_col = gaps,
                   annotation_colors = list(cluster = cluster.cols),
                   annotation_col = df.grp['cluster'])
dev.off()


rm(list = ls())
df.exp = read.delim2("cleandata/GSE107499.exp.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> as.data.frame()
df.im = deconvo_cibersort(df.exp, arrays = T, perm = 1)
df.im = df.im[1:23]
colnames(df.im) = str_remove(colnames(df.im), "_CIBERSORT")
write.table(df.im, "interm/cibersort.GSE107499.xls", sep = "\t", quote = F, row.names = F)


rm(list = ls())
df.im = read.delim2("interm/cibersort.GSE107499.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> as.data.frame()
df.grp = read.delim2("interm/consensus.cluster.tsv")
df.im = df.im[df.grp$sample,]

p.diff = apply(df.im, 2, \(x){
  kruskal.test(x~df.grp$cluster)$p.value
})

df.diff = data.frame(cell = colnames(df.im),
                     pval = p.diff,
                     padj = p.adjust(p.diff, "BH"))
df.diff = df.diff[order(df.diff$padj, decreasing = F),]


rm(list = ls())
df.exp = read.delim2("cleandata/GSE107499.exp.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> as.data.frame()
df.m = msigdbr()
df.m = df.m[df.m$gs_name %in% c("BIOCARTA_IL6_PATHWAY",
                                "HALLMARK_INFLAMMATORY_RESPONSE",
                                "GOBP_TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY"),]
gs = split(df.m$gene_symbol, df.m$gs_name)

gs.2 = read.gmt("rawdata/GOBP_INTERFERON_MEDIATED_SIGNALING_PATHWAY.v2025.1.Hs.gmt")
gs$IFN = gs.2$gene

gsvaPar = ssgseaParam(as.matrix(df.exp), gs)
res.gsva = gsva(gsvaPar)
res.gsva = as.data.frame(res.gsva)
res.gsva = t(res.gsva) |> as.data.frame()

df.grp = read.delim("interm/consensus.cluster.tsv")
res.gsva = res.gsva[df.grp$sample,]

colnames(res.gsva) = c("IL6", "TNFa", "Inflam", "IFN")
p.diff = apply(res.gsva, 2, \(x){
  kruskal.test(x~df.grp$cluster)$p.value
})

df.diff = data.frame(pathway = colnames(res.gsva),
                     pval = p.diff,
                     padj = p.adjust(p.diff, "BH"))
df.diff = df.diff[order(df.diff$padj, decreasing = F),]

res.gsva$sample = rownames(res.gsva)
res.gsva = reshape2::melt(res.gsva, id.var = "sample", variable.name = "path")
res.gsva$padj = df.diff$padj[match(res.gsva$path, df.diff$pathway)]
res.gsva$path = c("IL6"="IL6 Signaling",
                  "TNFa"="TNFα signaling",
                  "Inflam"="Inflammatory response",
                  "IFN"="IFN signaling")[res.gsva$path]
res.gsva$lab = paste0(res.gsva$path, "\np.adj = ", signif(res.gsva$padj,3))
res.gsva$cluster = df.grp$cluster[match(res.gsva$sample, df.grp$sample)]

ggplot(res.gsva, aes(x = cluster, y = value)) + 
  geom_boxplot(aes(color = cluster)) + 
  facet_wrap(~lab, scale = "free") +
  scale_color_d3(name = NULL) +
  theme_bw() + 
  xlab(NULL) + ylab("Score") +
  theme(axis.text.x.bottom = element_blank(), 
        axis.ticks.length.x.bottom = unit(0,"mm"))
ggsave("04.Cluster/06.Pathway.png", width = 6, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("04.Cluster/06.Pathway.pdf", width = 6, height = 7, units = "in", dpi = 300, bg = "white")

