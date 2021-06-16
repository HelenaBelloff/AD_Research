#!/usr/bin/env Rscript
library(GOtest) ##from Minghui Wang
library(msigdb) ##from Minghui Wang

options(stringsAsFactors=FALSE)

GOsets = c('C5.BP','C5.CC', 'C5.MF')
gosets_genes = msigdb.genesets(sets=GOsets, type='symbols', species='human',return.data.frame=T)

genes <- read.table("Gene_Names.txt", "\t", header = TRUE)
universe = genes$gene_names


ad_deg <- read.delim("ad_spearman.txt", "\t", header = TRUE)


result = GOtest(x=ad_deg[,c('Gene', 'group')], go=gosets_genes, query.population=universe, background='query', name.x='significant_corr', name.go='GOsets', method='hypergeometric')

write.table(result,"/sc/arion/projects/zhangb03a/belloh02/FULL_AD_ANALYSIS/Full_AD_Spearman/Full_AD_Spearman_Hypergeom.tsv", sep="\t", quote=F)