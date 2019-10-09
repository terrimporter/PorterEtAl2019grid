# Teresita M. Porter, Oct. 9, 2019

library(vegan)
library(ggplot2)
library(ggpubr)

###################################################################
# Edit rarecurve function to remove the horizontal lines
###################################################################

rarecurve2 <- function (x, step = 1, sample, xlab = "Sample Size", ylab = "Species", 
                        label = TRUE, col, lty, ...) 
{
  x <- as.matrix(x)
  if (!identical(all.equal(x, round(x)), TRUE)) 
    stop("function accepts only integers (counts)")
  if (missing(col)) 
    col <- par("col")
  if (missing(lty)) 
    lty <- par("lty")
  tot <- rowSums(x)
  S <- specnumber(x)
  if (any(S <= 0)) {
    message("empty rows removed")
    x <- x[S > 0, , drop = FALSE]
    tot <- tot[S > 0]
    S <- S[S > 0]
  }
  nr <- nrow(x)
  col <- rep(col, length.out = nr)
  lty <- rep(lty, length.out = nr)
  out <- lapply(seq_len(nr), function(i) {
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) 
      n <- c(n, tot[i])
    drop(rarefy(x[i, ], n))
  })
  Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(out, max)
  plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab, 
       type = "n", ...)
  if (!missing(sample)) {
    abline(v = sample) 
    rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"), 
                                           y = z, xout = sample, rule = 1)$y)
    #    abline(h = rare, lwd = 0.5) #turn off horizontal lines
  }
  for (ln in seq_along(out)) {
    N <- attr(out[[ln]], "Subsample")
    lines(N, out[[ln]], col = col[ln], lty = lty[ln], ...)
  }
  if (label) { 
    ordilabel(cbind(tot, S), labels = rownames(x), ...)
  }
  invisible(out)
}
###################################################################

# Read infile prepared by python script
## LV2016_2.csv is the clean taxonomic assignment table
A<-read.table(file="LV2016_2.csv", head=TRUE, sep=",")

# Select phylum Arthropoda only
Arth_df<-A[A$phylum=="Arthropoda",]

#split out ILC and NZC datasets
# ESV 1C1E
ILC_A<-Arth_df[Arth_df$site=="ILC45" & Arth_df$extractions==1,]
NZC_A<-Arth_df[Arth_df$site=="NZC85" & Arth_df$extractions==1,]
# ESV 1C3E
ILC_D<-Arth_df[Arth_df$site=="ILC45" & Arth_df$cores==1 & Arth_df$extractions==3,]
NZC_D<-Arth_df[Arth_df$site=="NZC85" & Arth_df$cores==1 & Arth_df$extractions==3,]
# ESV XC3E
ILC_G<-Arth_df[Arth_df$site=="ILC45" & Arth_df$cores!=1 & Arth_df$extractions==3,]
NZC_G<-Arth_df[Arth_df$site=="NZC85" & Arth_df$cores!=1 & Arth_df$extractions==3,]

# Pivot to make matrices for vegan, ESVs in rows, sites in columns, reads in cells
ILC_A.matrix<-dcast(ILC_A, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)
NZC_A.matrix<-dcast(NZC_A, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)

ILC_D.matrix<-dcast(ILC_D, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)
NZC_D.matrix<-dcast(NZC_D, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)

ILC_G.matrix<-dcast(ILC_G, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)
NZC_G.matrix<-dcast(NZC_G, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)

# Move marker_OTU to row names
rownames(ILC_A.matrix)<-ILC_A.matrix$marker_OTU
ILC_A.matrix<-ILC_A.matrix[,-1]
rownames(NZC_A.matrix)<-NZC_A.matrix$marker_OTU
NZC_A.matrix<-NZC_A.matrix[,-1]

rownames(ILC_D.matrix)<-ILC_D.matrix$marker_OTU
ILC_D.matrix<-ILC_D.matrix[,-1]
rownames(NZC_D.matrix)<-NZC_D.matrix$marker_OTU
NZC_D.matrix<-NZC_D.matrix[,-1]

rownames(ILC_G.matrix)<-ILC_G.matrix$marker_OTU
ILC_G.matrix<-ILC_G.matrix[,-1]
rownames(NZC_G.matrix)<-NZC_G.matrix$marker_OTU
NZC_G.matrix<-NZC_G.matrix[,-1]

# Transpose to get sites in rows, ESVs in columns, reads in cells
ILC_A.matrix.t<-t(ILC_A.matrix)
NZC_A.matrix.t<-t(NZC_A.matrix)

ILC_D.matrix.t<-t(ILC_D.matrix)
NZC_D.matrix.t<-t(NZC_D.matrix)

ILC_G.matrix.t<-t(ILC_G.matrix)
NZC_G.matrix.t<-t(NZC_G.matrix)

#remove columns with only zeros
ILC_A_notnull<-ILC_A.matrix.t[,colSums(ILC_A.matrix.t) !=0]
NZC_A_notnull<-NZC_A.matrix.t[,colSums(NZC_A.matrix.t) !=0]

ILC_D_notnull<-ILC_D.matrix.t[,colSums(ILC_D.matrix.t) !=0]
NZC_D_notnull<-NZC_D.matrix.t[,colSums(NZC_D.matrix.t) !=0]

ILC_G_notnull<-ILC_G.matrix.t[,colSums(ILC_G.matrix.t) !=0]
NZC_G_notnull<-NZC_G.matrix.t[,colSums(NZC_G.matrix.t) !=0]

#remove rows with only zeros
ILC_A_notnull2<-ILC_A_notnull[rowSums(ILC_A_notnull) !=0,]
NZC_A_notnull2<-NZC_A_notnull[rowSums(NZC_A_notnull) !=0,]

ILC_D_notnull2<-ILC_D_notnull[rowSums(ILC_D_notnull) !=0,]
NZC_D_notnull2<-NZC_D_notnull[rowSums(NZC_D_notnull) !=0,]

ILC_G_notnull2<-ILC_G_notnull[rowSums(ILC_G_notnull) !=0,]
NZC_G_notnull2<-NZC_G_notnull[rowSums(NZC_G_notnull) !=0,]

#calculate 15th percentile for rrarefy function
ILC_A_15percentile<-quantile(rowSums(ILC_A_notnull2), prob=0.15)
# 2843.8
NZC_A_15percentile<-quantile(rowSums(NZC_A_notnull2), prob=0.15)
# 3170.9

ILC_D_15percentile<-quantile(rowSums(ILC_D_notnull2), prob=0.15)
# 2641.85
NZC_D_15percentile<-quantile(rowSums(NZC_D_notnull2), prob=0.15)
# 2302.3

ILC_G_15percentile<-quantile(rowSums(ILC_G_notnull2), prob=0.15)
# 2540.6
NZC_G_15percentile<-quantile(rowSums(NZC_G_notnull2), prob=0.15)
# 4317.4

###################################################################
##### Plot rarefaction curves
###################################################################

set.seed(1234)

#print rarefaction curves for each sample to assess read coverage per sample, indicate 15th percentile
pdf("SF3_rarecurves.pdf")
par(mfrow=c(3,3))

rarecurve2(ILC_A_notnull2, sample=ILC_A_15percentile, step=100, main="ILC - 1C1E", xlab="Reads", ylab="ESVs", col=c("#4DAF4A","#FF7F00","#377EB8",
                                                                                                                    "#4DAF4A","#FF7F00",
                                                                                                                    rep(c("#4DAF4A","#FF7F00","#377EB8"),6)),
                                                                                                                    cex=1, cex.main=1, cex.lab=1, cex.axis=1, label=FALSE)
rarecurve2(ILC_D_notnull2, sample=ILC_D_15percentile, step=100, main="ILC - 1C3E", xlab="Reads", ylab="ESVs", col=rep(c("#4DAF4A","#FF7F00","#377EB8"),36), 
                                                                                                                    cex=1, cex.main=1, cex.lab=1, cex.axis=1, label=FALSE)
rarecurve2(ILC_G_notnull2, sample=ILC_G_15percentile, step=100, main="ILC - XC3E", xlab="Reads", ylab="ESVs", col=rep(c("#4DAF4A","#FF7F00","#377EB8"),20),
                                                                                                                    cex=1, cex.main=1, cex.lab=1, cex.axis=1, label=FALSE)

rarecurve2(NZC_A_notnull2, sample=NZC_A_15percentile, step=100, main="NZC - 1C1E", xlab="Reads", ylab="ESVs", col=rep(c("#4DAF4A","#FF7F00","#377EB8"),8),
                                                                                                                    cex=1, cex.main=1, cex.lab=1, cex.axis=1, label=FALSE)
rarecurve2(NZC_D_notnull2, sample=NZC_D_15percentile, step=100, main="NZC - 1C3E", xlab="Reads", ylab="ESVs", col=c(rep(c("#4DAF4A","#FF7F00","#377EB8"),3),
                                                                                                                    "#4DAF4A","#FF7F00",
                                                                                                                    rep(c("#4DAF4A","#FF7F00","#377EB8"),2),
                                                                                                                    "#FF7F00","#377EB8",
                                                                                                                    rep(c("#4DAF4A","#FF7F00","#377EB8"),3),
                                                                                                                    "#4DAF4A","#FF7F00",
                                                                                                                    rep(c("#4DAF4A","#FF7F00","#377EB8"),5),
                                                                                                                    "#4DAF4A","#FF7F00",
                                                                                                                    rep(c("#4DAF4A","#FF7F00","#377EB8"),16),
                                                                                                                    "#4DAF4A","#FF7F00",
                                                                                                                    rep(c("#4DAF4A","#FF7F00","#377EB8"),2)),
                                                                                                                    cex=1, cex.main=1, cex.lab=1, cex.axis=1, label=FALSE)
rarecurve2(NZC_G_notnull2, sample=NZC_G_15percentile, step=100, main="NZC - XC3E", xlab="Reads", ylab="ESVs", col=c("#4DAF4A",
                                                                                                                    "#4DAF4A",
                                                                                                                    "#4DAF4A","#FF7F00","#377EB8",
                                                                                                                    "#FF7F00","#377EB8",
                                                                                                                    "#FF7F00","#377EB8",
                                                                                                                    "#FF7F00","#377EB8",
                                                                                                                    rep(c("#4DAF4A","#FF7F00","#377EB8"),16),
                                                                                                                    "#4DAF4A"),
                                                                                                                    cex=1, cex.main=1, cex.lab=1, cex.axis=1, label=FALSE)

plot.new()
frame()
plot.new()
legend("topleft", 
       legend=c("Bryophyte","Organic","Mineral"), 
       lty=1,
       col=c("#4DAF4A","#377EB8","#FF7F00"),
       box.lty=0)

dev.off()
