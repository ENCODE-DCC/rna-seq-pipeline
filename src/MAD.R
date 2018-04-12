## input two table, each is from one rep
## version: 1.0

##table1 and table2 are the filenames with data
organizeExp <- function(table1,table2,col1,gene1=NA,tx1=NA,align=c("gene","tx"),col2=col1,gene2=gene1,tx2=tx1){
    align <- match.arg(align)
    t1 <- read.delim(table1,header=T,stringsAsFactors=F)
    t2 <- read.delim(table2,header=T,stringsAsFactors=F)
    exp1 <- t1[,col1]
    exp2 <- t2[,col2]
    if(!is.na(gene1)){
        geneid1 <- t1[,gene1]
        geneid2 <- t2[,gene2]
    }
    if(!is.na(tx1)){
        txid1 <- t1[,tx1]
        txid2 <- t2[,tx2]
    }
    if(align=="gene"){
        if(length(unique(geneid1))!=length(geneid1)) stop("not unique")
        exp2 <- exp2[match(geneid1,geneid2)]
        geneid2 <- geneid2[match(geneid1,geneid2)]
        if(!identical(geneid1,geneid2)) stop("not identical")
        data.frame(gene=geneid1,rep1=exp1,rep2=exp2)
    }else{
        if(length(unique(txid1))!=length(txid1)) stop("not unique")
        exp2 <- exp2[match(txid1,txid2)]
        txid2 <- txid2[match(txid1,txid2)]
        if(!identical(txid1,txid2)) stop("not identical")
        data.frame(transcript=txid1,rep1=exp1,rep2=exp2)
    }
}

### This cutoff depends on the units of what we are computing

Acutoff <-0  ##this means we ignore FPKM < 1 when computing mean

###para will have filenames supplied by
para <- commandArgs(trailingOnly = TRUE)
reps <- organizeExp(para[1],para[2],7,1)
nozero <- which(reps$rep1!=0 | reps$rep2!=0)
reps_part <- reps[nozero,]
logrep1 <- log2(reps$rep1[nozero])
logrep2 <- log2(reps$rep2[nozero])
A <- (logrep1+logrep2)/2
M <- logrep1-logrep2

cat("{\n")
cat("\"MAD of log ratios\":", round(median(abs(M)[A>Acutoff])*1.4826,3),",","\n")

##if you want to compute pearson and spearman on same data it's easy:
cat("\"Pearson correlation\":",cor(logrep1[A>Acutoff],logrep2[A>Acutoff]),",","\n")
cat("\"Spearman correlation\":",cor(logrep1[A>Acutoff],logrep2[A>Acutoff],method="spearman"),",","\n")
cat("\"SD of log ratios\":", round(sqrt(mean(M[A>Acutoff]^2)),3),"\n")
cat("}\n")

###if you want to make plot
bitmap("MAplot.png")
plot(A,M)
##
