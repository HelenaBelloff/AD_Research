#!/usr/bin/env Rscript
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
options(stringsAsFactors=F)
library(DESeq2)
library(limma)
library(edgeR)
library(qvalue)

wkdir = "/sc/arion/projects/zhangb03a/belloh02/Regulation";
setwd(wkdir)

deg_out = "DEGS_Full_AD/"

data.file <- "genes.tsv"
annot.col <- 1

print("loading data.df")
data.Df <- read.delim(file = data.file,sep = "\t",header = T)
datExpr <- as.matrix(data.Df[,-(annot.col)]);
rownames(datExpr) <- data.Df[,1] # should be gene names

datExpr = datExpr[complete.cases(datExpr),]
print("finished loading data.")


classification = read.delim2(file="new_meta_BM36.tsv", header=T, sep="\t") # meta file containing info on clinical covariates
classification$cluster = classification$ADsubtypeclass
rownames(classification) = classification$SynapseId
classification = classification[match(colnames(datExpr),rownames(classification)),]
classification = classification[complete.cases(classification),]


datExpr = datExpr[,rownames(classification)]
classification = classification[colnames(datExpr),]


clusters = unique(classification$cluster)
clusters = clusters[clusters!=""]


logCPM = datExpr
design = model.matrix(~0+cluster,data=classification)

fit = lmFit(datExpr,design) # fitting groups together

for(cluster in clusters){
        if((cluster=="control")){next}
        print(cluster)
        contrast.matrix <- makeContrasts(paste0("cluster",cluster,"-clustercontrol"),levels=design)
        fit2 = contrasts.fit(fit,contrast.matrix)
        fit2 <- eBayes(fit2)
        diffexp_genes = topTable(fit2,coef=1,adjust.method = "BH",number=Inf)
        diffexp_genes$adj.P.Val = qvalue(p=diffexp_genes$P.Value,pi0.method="smoother")$qvalues
        print(nrow(diffexp_genes[diffexp_genes[,"adj.P.Val"]<0.05,]))
        write.table(diffexp_genes, file=paste(deg_out,"test_",cluster,"_new_deg_limma.tsv",sep=""),sep="\t",quote=F,row.names = T)
}


