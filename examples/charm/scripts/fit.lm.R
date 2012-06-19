options(scipen=14) # stop it from printing 1e6 instead of 1000000

fname = commandArgs(TRUE)[1]

df = read.delim("data/natgen2009.csv", header=T, sep=",")
df = df[seq(2, nrow(df), 2),]
df$TissueType = factor(df$TissueType)
df$DiseaseState = factor(df$DiseaseState)


for (tiss in c("colon", "frontalcortex", "liver", "spleen")) {
    df[,tiss] = df$TissueType == tiss
}


data = read.delim(fname, header=T, sep="\t")
stopifnot(all(df$SampleID == data$ID))
names = colnames(data)

write(paste("#chrom", "start", "end", "p.disease", "p.tissue", "p.colon",
        "p.frontalcortex", "p.liver", "p.spleen", sep="\t"), stdout())
for (i in 2:ncol(data)) {

    methp = data[,i]
    methp = log(methp / (1 - methp))

    m = lm(methp ~ DiseaseState + TissueType, data=df)
    s = summary(m)
    p.disease = format(s$coefficients["DiseaseStatetumor", "Pr(>|t|)"], digits=4)
    m = lm(methp ~ TissueType + DiseaseState, data=df)
    s = summary(m)
    r = drop1(m, ~ TissueType, test="F")
    p.tissue = format(r["TissueType", "Pr(>F)"], digits=4)

    s = summary(lm(methp ~ colon + DiseaseState + frontalcortex + liver + spleen, data=df))
    p.colon = format(s$coefficients["colonTRUE", "Pr(>|t|)"], digits=4)

    s = summary(lm(methp ~ frontalcortex + DiseaseState + colon + liver + spleen, data=df))
    p.fc = format(s$coefficients["frontalcortexTRUE", "Pr(>|t|)"], digits=4)

    s = summary(lm(methp ~ liver + frontalcortex + DiseaseState + colon + spleen, data=df))
    p.liver = format(s$coefficients["liverTRUE", "Pr(>|t|)"], digits=4)

    s = summary(lm(methp ~ spleen + liver + frontalcortex + DiseaseState + colon, data=df))
    p.spleen = format(s$coefficients["spleenTRUE", "Pr(>|t|)"], digits=4)

    chrom_start = unlist(strsplit(names[i], "_", fixed=TRUE))
    start = as.numeric(chrom_start[2])

    write(paste(chrom_start[1], start - 1, start + 50, p.disease, p.tissue, p.colon, p.fc, p.liver, p.spleen, sep="\t"), stdout())
    #stopifnot(p.asthma == p.drop.asthma)

}

