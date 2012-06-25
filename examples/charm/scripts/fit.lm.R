options(scipen=14) # stop it from printing 1e6 instead of 1000000

fname = commandArgs(TRUE)[1]

df = read.delim("data/natgen2009.csv", header=T, sep=",")
df = df[seq(2, nrow(df), 2),]
df$TissueType = factor(df$TissueType)
df$DiseaseState = factor(df$DiseaseState)

set.seed(42)
for(i in seq(1, 20)){
    order1 = sample(1:nrow(df), nrow(df))
    df[,paste("shuff_tissue", i, sep="_")] = df$TissueType[order1]
    df[,paste("shuff_disease", i, sep="_")] = df$DiseaseState[order1]

    #shuff_df = df[sample(1:nrow(df), nrow(df)),]
}

for (tiss in c("colon", "frontalcortex", "liver", "spleen")) {
    df[,tiss] = df$TissueType == tiss
}


data = read.delim(fname, header=T, sep="\t")
stopifnot(all(df$SampleID == data$ID))
names = colnames(data)


write(paste("#chrom", "start", "end", "p.disease", "p.tissue", 
    paste("p-shuff-tiss", 1:20, sep="-", collapse="\t"),
    paste("p-shuff-disease", 1:20, sep="-", collapse="\t"),
 #"p.colon", "p.frontalcortex", "p.liver", "p.spleen", 
        sep="\t"), stdout())
for (i in 2:ncol(data)) {

    methp = data[,i]
    methp = log(methp / (1 - methp))

    m = lm(methp ~ DiseaseState + TissueType, data=df)
    r = drop1(m, ~ ., test="F")
    p.tissue = format(r["TissueType", "Pr(>F)"], digits=4)
    p.disease = format(r["DiseaseState", "Pr(>F)"], digits=4)

    p.shuff.tissue = rep(NA, 20)
    p.shuff.disease = rep(NA, 20)
    # shuffled
    for(i in 1:20){
        m = lm(as.formula(paste("methp ~ shuff_disease_", i, " + shuff_tissue_", i, sep="")), data=df)
        r = drop1(m, ~ ., test="F")
        p.shuff.tissue[i] = format(r[paste("shuff_tissue_", i, sep=""), "Pr(>F)"], digits=4)
        p.shuff.disease[i] = format(r[paste("shuff_disease_", i, sep=""), "Pr(>F)"], digits=4)
    }


    chrom_start = unlist(strsplit(names[i], "_", fixed=TRUE))
    start = as.numeric(chrom_start[2])

    write(paste(chrom_start[1], start - 1, start + 50, p.disease, p.tissue,
        paste(p.shuff.disease, collapse="\t"), 
        paste(p.shuff.tissue, collapse="\t")
        , sep="\t"), stdout())
    #stopifnot(p.asthma == p.drop.asthma)

}

