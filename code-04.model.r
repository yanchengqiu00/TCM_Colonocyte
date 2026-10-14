rm(list = ls())
library(tidyverse)
library(glmnet)
library(ggsci)
library(ggpubr)
library(ggvenn)
library(pROC)
library(Seurat)

df.exp = read.delim2("cleandata/GSE107499.exp.xls", row.names = 1)
gs = read.delim2("02.Target/05.Intersection.xls")
gs = gs$gene
df.pheno = read.delim2("cleandata/GSE107499.pheno.xls")

df.exp = df.exp[gs, df.pheno$geo_accession]
df.exp = apply(df.exp, c(1,2), as.numeric) |> t() |> as.data.frame()

set.seed(2025)
res.lasso = cv.glmnet(as.matrix(df.exp), df.pheno$group, family = "binomial", nfolds = 10, alpha = 1)
plot(res.lasso)
ggsave("03.Model/01.GSE107499.Lasso.CV.png", plot(res.lasso), width = 8, height = 7, dpi = 300, units = "in", bg = "white")
ggsave("03.Model/01.GSE107499.Lasso.CV.pdf", plot(res.lasso), width = 8, height = 7, dpi = 300, units = "in", bg = "white")
ggsave("03.Model/02.GSE107499.Lasso.Coef.png", plot(res.lasso$glmnet.fit, xvar = 'lambda'), width = 8, height = 7, dpi = 300, units = "in", bg = "white")
ggsave("03.Model/02.GSE107499.Lasso.Coef.pdf", plot(res.lasso$glmnet.fit, xvar = 'lambda'), width = 8, height = 7, dpi = 300, units = "in", bg = "white")
res.lasso$lambda.min # 0.007481599

df.coef = coef(res.lasso, s = res.lasso$lambda.min)
df.coef = cbind(gene = rownames(df.coef), coefficient = df.coef[,1]) |> as.data.frame()
df.coef$coefficient = as.numeric(df.coef$coefficient)
df.coef = subset(df.coef, coefficient != 0)
df.coef = df.coef[-1,]
df.coef = df.coef[order(df.coef$coefficient, decreasing = F),]
df.coef$gene = factor(df.coef$gene, levels = df.coef$gene)

ggplot(df.coef, aes(y = gene, x = coefficient)) + 
  geom_bar(aes(fill = gene), stat = "identity") + 
  geom_text(aes(label = signif(coefficient,3)), hjust = 0, data = subset(df.coef, coefficient>0)) +
  geom_text(aes(label = signif(coefficient,3)), hjust = 1, data = subset(df.coef, coefficient<0)) +
  scale_fill_d3("category20") +
  ggtitle("GSE107499") +
  xlim(c(-5.5,4)) + 
  ylab("Gene") + xlab("Coefficient") +
  guides(fill = guide_none()) +
  theme_classic2() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("03.Model/03.GSE107499.Coefficient.png", width = 6, height = 4, dpi = 300, units = "in", bg = "white")
ggsave("03.Model/03.GSE107499.Coefficient.pdf", width = 6, height = 4, dpi = 300, units = "in", bg = "white")
write.table(df.coef, "interm/GSE107499.coef.txt", quote = F, sep = "\t", row.names = F)


rm(list = ls())
df.exp = read.delim2("cleandata/GSE87466.exp.xls", row.names = 1)
gs = read.delim2("02.Target/05.Intersection.xls")
gs = gs$gene
df.pheno = read.delim2("cleandata/GSE87466.pheno.xls")

df.exp = df.exp[gs, df.pheno$geo_accession]
df.exp = apply(df.exp, c(1,2), as.numeric) |> t() |> as.data.frame()

set.seed(2025)
res.lasso = cv.glmnet(as.matrix(df.exp), df.pheno$group, family = "binomial", nfolds = 10, alpha = 1)
plot(res.lasso)
ggsave("03.Model/04.GSE87466.Lasso.CV.png", plot(res.lasso), width = 8, height = 7, dpi = 300, units = "in", bg = "white")
ggsave("03.Model/04.GSE87466.Lasso.CV.pdf", plot(res.lasso), width = 8, height = 7, dpi = 300, units = "in", bg = "white")
ggsave("03.Model/05.GSE87466.Lasso.Coef.png", plot(res.lasso$glmnet.fit, xvar = 'lambda'), width = 8, height = 7, dpi = 300, units = "in", bg = "white")
ggsave("03.Model/05.GSE87466.Lasso.Coef.pdf", plot(res.lasso$glmnet.fit, xvar = 'lambda'), width = 8, height = 7, dpi = 300, units = "in", bg = "white")
res.lasso$lambda.min # 0.01531288

df.coef = coef(res.lasso, s = res.lasso$lambda.min)
df.coef = cbind(gene = rownames(df.coef), coefficient = df.coef[,1]) |> as.data.frame()
df.coef$coefficient = as.numeric(df.coef$coefficient)
df.coef = subset(df.coef, coefficient != 0)
df.coef = df.coef[-1,]
df.coef = df.coef[order(df.coef$coefficient, decreasing = F),]
df.coef$gene = factor(df.coef$gene, levels = df.coef$gene)

ggplot(df.coef, aes(y = gene, x = coefficient)) + 
  geom_bar(aes(fill = gene), stat = "identity") + 
  geom_text(aes(label = signif(coefficient,3)), hjust = 0, data = subset(df.coef, coefficient>0)) +
  geom_text(aes(label = signif(coefficient,3)), hjust = 1, data = subset(df.coef, coefficient<0)) +
  scale_fill_d3("category20") +
  xlim(c(-2.5,3.5)) + 
  ggtitle("GSE87466") +
  ylab("Gene") + xlab("Coefficient") +
  guides(fill = guide_none()) +
  theme_classic2() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("03.Model/06.GSE87466.Coefficient.png", width = 6, height = 4, dpi = 300, units = "in", bg = "white")
ggsave("03.Model/06.GSE87466.Coefficient.pdf", width = 6, height = 4, dpi = 300, units = "in", bg = "white")
write.table(df.coef, "interm/GSE87466.coef.txt", quote = F, sep = "\t", row.names = F)


rm(list = ls())
gs.1 = read.delim2("interm/GSE107499.coef.txt")
gs.2 = read.delim2("interm/GSE87466.coef.txt")

ggvenn(list(GSE107499 = gs.1$gene, GSE87466 = gs.2$gene), 
       stroke_size = 0.3, set_name_size = 6, text_size = 5, show_percentage = F) + 
  ylim(c(-2,2)) + 
  scale_fill_bmj()
ggsave("03.Model/07.Lasso.Venn.png", width = 6, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("03.Model/07.Lasso.Venn.pdf", width = 6, height = 4, units = "in", dpi = 300, bg = "white")

inter = intersect(gs.1$gene, gs.2$gene)
write.table(data.frame(gene = inter), "03.Model/08.Selected.xls", quote = F, sep = "\t", row.names = F)


rm(list = ls())
df.exp = read.delim2("cleandata/GSE107499.exp.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> as.data.frame()
gs = read.delim2("03.Model/08.Selected.xls")
gs = gs$gene
df.pheno = read.delim2("cleandata/GSE107499.pheno.xls")
df.exp = df.exp[gs,df.pheno$geo_accession] |> t() |> as.data.frame()

res.roc = apply(df.exp, 2, function(x){
  roc(df.pheno$group, x)
})

res.auc = sapply(res.roc, pROC::auc)
res.auc = data.frame(Model = names(res.auc), AUC = res.auc)
res.auc$AUC = round(res.auc$AUC, 3)
res.auc$lab = paste0(res.auc$Model, " (AUC = ", res.auc$AUC, ")")
res.auc = res.auc[order(res.auc$AUC, decreasing = T),]
res.roc = res.roc[res.auc$Model]
names(res.roc) = res.auc$lab

ggroc(res.roc, legacy.axes = T, size = 0.8) + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "grey50", linetype = 2) +
  theme_bw() + 
  scale_color_manual(values = IOBR::palette1, 
                     guide = guide_legend(ncol = 1)) +
  ggtitle("GSE107499") +
  theme(panel.border = element_rect(color = "black", fill = "transparent"), 
        aspect.ratio = 1, legend.title = element_blank(), 
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"), 
        plot.title = element_text(hjust = 0.5))
ggsave("03.Model/09.GSE107499.ROC.png", width = 6, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("03.Model/09.GSE107499.ROC.pdf", width = 6, height = 5, units = "in", dpi = 300, bg = "white")


rm(list = ls())
df.exp = read.delim2("cleandata/GSE87466.exp.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> as.data.frame()
gs = read.delim2("03.Model/08.Selected.xls")
gs = gs$gene
df.pheno = read.delim2("cleandata/GSE87466.pheno.xls")
df.exp = df.exp[gs,df.pheno$geo_accession] |> t() |> as.data.frame()

res.roc = apply(df.exp, 2, function(x){
  roc(df.pheno$group, x)
})

res.auc = sapply(res.roc, pROC::auc)
res.auc = data.frame(Model = names(res.auc), AUC = res.auc)
res.auc$AUC = round(res.auc$AUC, 3)
res.auc$lab = paste0(res.auc$Model, " (AUC = ", res.auc$AUC, ")")
res.auc = res.auc[order(res.auc$AUC, decreasing = T),]
res.roc = res.roc[res.auc$Model]
names(res.roc) = res.auc$lab

ggroc(res.roc, legacy.axes = T, size = 0.8) + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "grey50", linetype = 2) +
  theme_bw() + 
  scale_color_manual(values = IOBR::palette1, 
                     guide = guide_legend(ncol = 1)) +
  ggtitle("GSE87466") +
  theme(panel.border = element_rect(color = "black", fill = "transparent"), 
        aspect.ratio = 1, legend.title = element_blank(), 
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"), 
        plot.title = element_text(hjust = 0.5))
ggsave("03.Model/10.GSE87466.ROC.png", width = 6, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("03.Model/10.GSE87466.ROC.pdf", width = 6, height = 5, units = "in", dpi = 300, bg = "white")


rm(list = ls())
df.exp = read.delim2("cleandata/GSE107499.exp.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> as.data.frame()
df.pheno = read.delim2("cleandata/GSE107499.pheno.xls")
gs = read.delim2("03.Model/08.Selected.xls")
gs = gs$gene
df.exp = df.exp[gs,df.pheno$geo_accession]

df.diff = read.delim2("03.DEG/01.DEG.xls")
df.exp = cbind(gene = rownames(df.exp), df.exp)
df.exp = reshape2::melt(df.exp, id.var = "gene", variable.name = "sample")
df.exp$group = df.pheno$group[match(df.exp$sample, df.pheno$geo_accession)]
df.exp$padj = df.diff$adj.P.Val[match(df.exp$gene, df.diff$Gene)] |> as.numeric()
df.exp$lab = paste0(df.exp$gene, "\np.adj = ", signif(df.exp$padj,3))

ggplot(df.exp, aes(x = group, y = value)) + 
  geom_boxplot(aes(color = group)) + 
  facet_wrap(~lab, ncol = 5, scale = "free") + 
  scale_color_d3(name = NULL) +
  xlab(NULL) + 
  ylab("Expression") +
  ggtitle("GSE107499") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "top",
        axis.text.x.bottom = element_blank(), 
        axis.ticks.length.x.bottom = unit(0,"mm"))
ggsave("03.Model/11.GSE107499.Boxplot.png", width = 6, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("03.Model/11.GSE107499.Boxplot.pdf", width = 6, height = 4, units = "in", dpi = 300, bg = "white")


rm(list = ls())
df.exp = read.delim2("cleandata/GSE87466.exp.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> as.data.frame()
df.pheno = read.delim2("cleandata/GSE87466.pheno.xls")
gs = read.delim2("03.Model/08.Selected.xls")
gs = gs$gene
df.exp = df.exp[gs,df.pheno$geo_accession]

p.diff = apply(df.exp, 1, \(x){
  kruskal.test(x~df.pheno$group)$p.value
})
df.diff = data.frame(gene = rownames(df.exp), 
                     p.val = p.diff,
                     p.adj = p.adjust(p.diff, method = "BH"))

df.exp = cbind(gene = rownames(df.exp), df.exp)
df.exp = reshape2::melt(df.exp, id.var = "gene", variable.name = "sample")
df.exp$group = df.pheno$group[match(df.exp$sample, df.pheno$geo_accession)]
df.exp$padj = df.diff$p.adj[match(df.exp$gene, df.diff$gene)]
df.exp$lab = paste0(df.exp$gene, "\np.adj = ", signif(df.exp$padj,3))

ggplot(df.exp, aes(x = group, y = value)) + 
  geom_boxplot(aes(color = group)) + 
  facet_wrap(~lab, ncol = 5, scale = "free") + 
  scale_color_d3(name = NULL) +
  xlab(NULL) + 
  ylab("Expression") +
  ggtitle("GSE87466") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "top",
        axis.text.x.bottom = element_blank(), 
        axis.ticks.length.x.bottom = unit(0,"mm"))
ggsave("03.Model/12.GSE87466.Boxplot.png", width = 6, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("03.Model/12.GSE87466.Boxplot.pdf", width = 6, height = 4, units = "in", dpi = 300, bg = "white")


rm(list = ls())
df.tar = read.delim2("02.Target/01.Drug.Targets.xls")
gs = read.delim2("03.Model/08.Selected.xls")
df.tar = df.tar[df.tar$Gene %in% gs$gene,]

df.drug = df.tar[1:2]
df.gene = df.tar[2:3]
df.drug = df.drug[!duplicated(df.drug),]
df.gene = df.gene[!duplicated(df.gene),]
colnames(df.drug) = colnames(df.gene) = c("source","target")
df.out = rbind(df.drug, df.gene)
write.table(df.out, "interm/hub.mol.drug.txt", quote = F, sep = "\t", row.names = F)


rm(list = ls())
library(rms)
library(regplot)
df.exp = read.delim2("cleandata/GSE107499.exp.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> as.data.frame()
gs = read.delim2("03.Model/08.Selected.xls")
gs = gs$gene
df.pheno = read.delim2("cleandata/GSE107499.pheno.xls")
df.exp = df.exp[gs, df.pheno$geo_accession]
df.exp = t(df.exp) |> as.data.frame()
df.exp = apply(df.exp, 2, scale) |> as.data.frame()
rownames(df.exp) = df.pheno$sample
df.exp$group = df.pheno$group

dd = datadist(df.exp)
options(datadist="dd")
fit = lrm(group~TOP1+SLPI+AKR1C3, data = df.exp, x = T, y = T)
df.exp$pred = predict(fit, df.exp)
roc.train = roc(response = df.exp$group, predictor = df.exp$pred)

df.test = read.delim2("cleandata/GSE87466.exp.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> as.data.frame()
df.test.pheno = read.delim2("cleandata/GSE87466.pheno.xls")
df.test = df.test[gs, df.test.pheno$geo_accession] |> t() |> as.data.frame()
df.test = apply(df.test, 2, scale) |> as.data.frame()
rownames(df.test) = df.test.pheno$sample
df.test$group = df.test.pheno$group
df.test$pred = predict(fit, df.test)
roc.test = roc(response = df.test$group, predictor = df.test$pred)
res.roc = list(`Train set` = roc.train, `Test set` = roc.test)

res.auc = sapply(res.roc, pROC::auc)
res.auc = data.frame(Model = names(res.auc), AUC = res.auc)
res.auc$AUC = round(res.auc$AUC, 3)
res.auc$lab = paste0(res.auc$Model, " (AUC = ", res.auc$AUC, ")")
res.auc = res.auc[order(res.auc$AUC, decreasing = T),]
res.roc = res.roc[res.auc$Model]
names(res.roc) = res.auc$lab

ggroc(res.roc, legacy.axes = T, size = 0.8) + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "grey50", linetype = 2) +
  theme_bw() + 
  scale_color_manual(values = IOBR::palette1, 
                     guide = guide_legend(ncol = 1)) +
  theme(panel.border = element_rect(color = "black", fill = "transparent"), 
        aspect.ratio = 1, legend.title = element_blank(), 
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"), 
        plot.title = element_text(hjust = 0.5),
        legend.position = "inside", 
        legend.position.inside = c(0.75, 0.25))
ggsave("03.Model/14.Model.ROC.png", width = 5, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("03.Model/14.Model.ROC.pdf", width = 5, height = 5, units = "in", dpi = 300, bg = "white")


set.seed(2025)
cal.train = calibrate(fit, method='boot', B=100)
png("03.Model/15.Calibration.Train.png", width = 5, height = 5, units = "in", bg = "white", res = 300)
plot(cal.train,
     xlim = c(0,1),
     ylim = c(0,1),
     xlab = "Prediced Probability",
     ylab = "Observed Probability"
) 
dev.off()
pdf("03.Model/15.Calibration.Train.pdf", width = 5, height = 5)
plot(cal.train,
     xlim = c(0,1),
     ylim = c(0,1),
     xlab = "Prediced Probability",
     ylab = "Observed Probability"
) 
dev.off()

library(ggDCA)
ggplot(dca(fit, model.names = "Logistic Regression"), 
       linetype = 1, color = pal_d3()(3), lwd = 0.7) + 
  scale_y_continuous(limits = c(-0.52, 0.75)) +
  theme_classic2(base_size = 12) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.25,0.2),
        legend.title = element_blank(),
        legend.background = element_blank())
ggsave("03.Model/16.DCA.png", width = 6, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("03.Model/16.DCA.pdf", width = 6, height = 4, units = "in", dpi = 300, bg = "white")



