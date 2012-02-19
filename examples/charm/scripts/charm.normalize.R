library(charm)
require(BSgenome.Hsapiens.UCSC.hg18)

sample_description = commandArgs(TRUE)[1] # e.g. sample-description.txt
out_methp = commandArgs(TRUE)[2]  # e.g. methp.txt
colname = commandArgs(TRUE)[3] # e.g. Barcode
method = commandArgs(TRUE)[4] # quantile or sqn

pd = read.delim(sample_description, sep=",")
dat = readCharm(files=pd$Filename, sampleKey=pd,
                path="/vol2/home/brentp/src/combined-pvalues/examples/charm/data")

# need this since pmQuality is not exported.
#attach(getNamespace("charm"))
#qual = format(pmQuality(dat), digits=2)
#colnames(qual) = pd[seq(1, nrow(pd), 2), colname]
#rownames(qual) = paste(pmChr(dat), pmPosition(dat), sep="_")
#message(paste(out_methp, ".qual.txt", sep=""))
#write.table(qual, sep="\t", quote=F, file=paste(out_methp, ".qual.txt", sep=""), row.names=T)

ctrlIdx = getControlIndex(dat, subject=Hsapiens)

p  = methp(dat, controlIndex=ctrlIdx, betweenSampleNorm=method)
rm(ctrlIdx)
rownames(p) = paste(pmChr(dat), pmPosition(dat), sep="_")
message("before COLS")
colnames(p) = pd[seq(1, nrow(pd), 2), colname]
message("COLS")
rm(dat)

tp = format(t(p), nsmall=3, digits=1)
rm(p)

message("writing...")
write.table(tp, sep="\t", quote=FALSE, file=out_methp)

