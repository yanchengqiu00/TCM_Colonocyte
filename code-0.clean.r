rm(list = ls())
library(org.Hs.eg.db)
library(tidyverse)
library(GEOquery)
library(AnnoProbe)

gs = getGEO(filename = "rawdata/GSE107499_series_matrix.txt.gz", getGPL = F)
df.exp = Biobase::exprs(gs) |> as.data.frame()
df.exp = limma::normalizeBetweenArrays(df.exp) |> as.data.frame()
df.pheno = gs |> pData() |>  as.data.frame()
df.annot = read.delim2("rawdata/GPL15207-17536.txt.gz", comment.char = "#")
df.annot = subset(df.annot, !str_detect(Gene.Symbol, "//"))
df.exp$gene = df.annot$Gene.Symbol[match(rownames(df.exp), df.annot$ID)]
df.exp = subset(df.exp, gene != "")
df.exp = na.omit(df.exp)
df.exp = aggregate(.~gene, FUN = max, data = df.exp)
df.pheno$sample = df.pheno$geo_accession
df.pheno$group = ifelse(df.pheno$`tissue state:ch1`=="Lesional", "Disease","Control")
write.table(df.exp, "cleandata/GSE107499.exp.xls", sep = "\t", quote = F, row.names = F)
write.table(df.pheno, "cleandata/GSE107499.pheno.xls", sep = "\t", quote = F, row.names = F)


rm(list = ls())
gs = getGEO(filename = "rawdata/GSE75214_series_matrix.txt.gz", getGPL = F)
df.exp = Biobase::exprs(gs) |> as.data.frame()
df.exp = limma::normalizeBetweenArrays(df.exp) |> as.data.frame()
df.pheno = gs |> pData() |>  as.data.frame()
df.pheno$sample = df.pheno$geo_accession
df.pheno = df.pheno[df.pheno$`disease:ch1` != "CD",]
df.pheno$group = ifelse(df.pheno$`disease:ch1`=="ulcerative colitis", "Disease","Control")
df.exp = df.exp[df.pheno$sample]
df.annot = getGEO(filename = "rawdata/GPL6244.annot.gz") |> Table()
df.annot = subset(df.annot, !str_detect(`Gene symbol`, "//"))
df.exp$gene = df.annot$`Gene symbol`[match(rownames(df.exp), df.annot$ID)]
df.exp = subset(df.exp, gene != "")
df.exp = na.omit(df.exp)
df.exp = aggregate(.~gene, FUN = max, data = df.exp)
write.table(df.exp, "cleandata/GSE75214.exp.xls", sep = "\t", quote = F, row.names = F)
write.table(df.pheno, "cleandata/GSE75214.pheno.xls", sep = "\t", quote = F, row.names = F)


rm(list = ls())
gs = getGEO(filename = "rawdata/GSE87466_series_matrix.txt.gz", getGPL = F)
df.exp = Biobase::exprs(gs) |> as.data.frame()
df.exp = limma::normalizeBetweenArrays(df.exp) |> as.data.frame()
df.pheno = gs |> pData() |>  as.data.frame()
df.pheno$sample = df.pheno$geo_accession
df.pheno$group = ifelse(df.pheno$`disease:ch1`=="Ulcerative Colitis (UC)", "Disease","Control")
df.annot = getGEO(filename = "rawdata/GPL13158.annot.gz") |> Table()
df.annot = subset(df.annot, !str_detect(`Gene symbol`, "//"))
df.exp$gene = df.annot$`Gene symbol`[match(rownames(df.exp), df.annot$ID)]
df.exp = subset(df.exp, gene != "")
df.exp = na.omit(df.exp)
df.exp = aggregate(.~gene, FUN = max, data = df.exp)
write.table(df.exp, "cleandata/GSE87466.exp.xls", sep = "\t", quote = F, row.names = F)
write.table(df.pheno, "cleandata/GSE87466.pheno.xls", sep = "\t", quote = F, row.names = F)

