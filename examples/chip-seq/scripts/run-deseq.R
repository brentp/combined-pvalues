library(DESeq)
fname = commandArgs(TRUE)[1]
tbl = read.delim(fname, header=T, sep="\t")

cds = newCountDataSet(tbl,  c("A3", "A4"))
cds = estimateSizeFactors(cds)

# choice of these depends on de-seq version.
#cds = estimateVarianceFunctions(cds, method="blind", sharingMode="fit-only")
cds = estimateVarianceFunctions(cds, method="blind")

res = nbinomTest(cds, "A3", "A4")
write.table(res, row.names=F, sep="\t", quote=F)

#resvalid = res[!is.na(res$padj),]
#resSig = resvalid[ resvalid$padj < .05, ]

#write.table(resSig, row.names=F, sep="\t", quote=F)
