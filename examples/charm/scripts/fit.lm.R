options(scipen=14) # stop it from printing 1e6 instead of 1000000

fname = commandArgs(TRUE)[1]

df = read.delim("data/natgen2009.csv", header=T, sep=",")
df = df[seq(2, nrow(df), 2),]
df$TissueType = factor(df$TissueType)
df$DiseaseState = factor(df$DiseaseState)

set.seed(42)
shuff_df = df[sample(1:nrow(df), nrow(df)),]


for (tiss in c("colon", "frontalcortex", "liver", "spleen")) {
    df[,tiss] = df$TissueType == tiss
}


data = read.delim(fname, header=T, sep="\t")
stopifnot(all(df$SampleID == data$ID))
names = colnames(data)

write(paste("#chrom", "start", "end", "p.disease", "p.tissue", 
                                      "p.shuff.disease", "p.shuff.tissue",
 #"p.colon", "p.frontalcortex", "p.liver", "p.spleen", 
        sep="\t"), stdout())
for (i in 2:ncol(data)) {

    methp = data[,i]
    methp = log(methp / (1 - methp))

    m = lm(methp ~ DiseaseState + TissueType, data=df)
    r = drop1(m, ~ ., test="F")
    p.tissue = format(r["TissueType", "Pr(>F)"], digits=4)
    p.disease = format(r["DiseaseState", "Pr(>F)"], digits=4)

    # shuffled
    m = lm(methp ~ DiseaseState + TissueType, data=shuff_df)
    r = drop1(m, ~ ., test="F")
    p.shuff.tissue = format(r["TissueType", "Pr(>F)"], digits=4)
    p.shuff.disease = format(r["DiseaseState", "Pr(>F)"], digits=4)


    chrom_start = unlist(strsplit(names[i], "_", fixed=TRUE))
    start = as.numeric(chrom_start[2])

    write(paste(chrom_start[1], start - 1, start + 50, p.disease, p.tissue, p.shuff.disease, p.shuff.tissue, sep="\t"), stdout())
    #stopifnot(p.asthma == p.drop.asthma)

}

