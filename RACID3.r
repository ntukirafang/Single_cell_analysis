knitr::opts_chunk$set(echo = TRUE)
options(repos=structure(c(CRAN="http://cran.mirrors.hoobly.com/")))
install.packages("devtools")
library(devtools)
install_github("dgrun/RaceID3_StemID2_package")
library(RaceID)
require(Matrix)
require(RaceID)
setwd("C:/Users/jimmy.fang/Documents")
x <- readMM("matrix.mtx")
f <- read.csv("features.tsv",sep="\t",header=FALSE)
b <- read.csv("barcodes.tsv",sep="\t",header=FALSE)
rownames(x) <- f[,1]
colnames(x) <- b[,1]
sc <- SCseq(x)
sc <- filterdata(sc,mintotal=1000, minexpr=2)
sc <- compdist(sc,metric="pearson")
sc <- clustexp(sc)
plotsaturation(sc,disp=FALSE)
plotsaturation(sc,disp=TRUE)
plotjaccard(sc)
sc <- clustexp(sc,cln=15,sat=FALSE)
sc <- findoutliers(sc)
plotbackground(sc)
plotsensitivity(sc)
plotoutlierprobs(sc)
clustheatmap(sc)
sc <- comptsne(sc)
sc <- compfr(sc,knn=10)
sc <- compumap(sc)
plotmap(sc)
plotmap(sc,um=TRUE)
types <- sub("(\\_\\d+)$","", colnames(sc@ndata))
subset <- types[grep("AAACCTGCAGTATCTG-1",types)]
plotsymbolsmap(sc,types,subset=subset,fr=TRUE)
d  <- clustdiffgenes(sc,4,pvalue=.01)
dg <- d$dg
head(dg,35)
types <- sub("(\\_\\d+)$","", colnames(sc@ndata))
genes <- head(rownames(dg)[dg$fc>1],35)
plotmarkergenes(sc,genes,samples=types)
fractDotPlot(sc, genes, cluster=c(2,6,7,8,10), zsc=TRUE)
A <- names(sc@cpart)[sc@cpart %in% c(2,4)]
B <- names(sc@cpart)[sc@cpart %in% c(3)]
x <- diffexpnb(getfdata(sc,n=c(A,B)), A=A, B=B )
plotdiffgenesnb(x,pthr=.05,lthr=.5,mthr=-1,Aname="Cl.2",Bname="Cl.3,5",show_names=TRUE,padj=TRUE)
ltr <- Ltree(sc)
ltr <- compentropy(ltr)
ltr <- projcells(ltr,cthr=5,nmode=TRUE,fr=TRUE,knn=10)
ltr <- lineagegraph(ltr)
ltr <- comppvalue(ltr,pthr=0.01)
plotgraph(ltr,showCells=FALSE,showMap=TRUE)