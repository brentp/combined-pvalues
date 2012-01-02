library(DESeq)
fname = commandArgs(TRUE)[1]
tbl = read.delim(fname, header=T, sep="\t")

tbl = newCountDataSet(tbl,  c("A3", "A4"))
tbl = estimateSizeFactors(tbl)

# choice of these depends on de-seq version.
#cds = estimateVarianceFunctions(cds, method="blind", sharingMode="fit-only")
tbl = estimateVarianceFunctions(tbl, method="blind")

tbl = nbinomTest(tbl, "A3", "A4")
tbl = tbl[!is.na(tbl$pval), colnames(tbl) %in% c("id", "pval")]
write.table(tbl, row.names=F, sep="\t", quote=F)
