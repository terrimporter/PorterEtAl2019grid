# Teresita M. Porter, Feb. 13, 2019

library(vegan)
library(reshape2)
library(gridExtra)
library(grid)
library(ggplot2)
library(plyr)
library(data.table)
library(pairwiseAdonis)
library("car")
library(cowplot)

###################################################################
##### NEW FUNCTION FOR SCREE PLOT #################################
###################################################################

#new function edited from http://www.pmc.ucsc.edu/~mclapham/Rtips/ordination.htm
nmds.scree<-function(x) {
  plot(rep(1,10),replicate(10,metaMDS(x, autotransform=F,k=1)$stress/100),xlim=c(1,nrow(x)),ylim=c(0,0.5),xlab="Number dimensions",ylab="Stress",main="NMDS stress plot")
  # originally for(i in 1:(nrow(x)-2)) { but change to just 10 dim
  for(i in 1:10) {
    points(rep(i+1,10),replicate(10,metaMDS(x, autotransform=F,k=i+1)$stress/100))
  }
}

########################################################
##### NEW FUNCTIUON TO GET TARGET GRIDCOORD DF's
########################################################

# New function to grab target gridcoord rows (use for 2 + pooled cores)

get_target_gridcoord<- function (matrix, df) {
  
  out=NULL
  df.list<-c()
  for (i in 1:(ncol(matrix))) {
    gridcoord.list<-matrix[,i]
    
    for (j in 1:length(gridcoord.list)) {
      targetgridcoord<-gridcoord.list[[j]] 
      df.j<-df[df$gridcoord==targetgridcoord,]
      out=rbind(out, df.j)
    }
    df.list[[i]]<-out
  }
  return(df.list)
}

########################################################
##### NEW FUNCTION TO GET DISTANCE MATRICES
########################################################

get_distmat <- function (list) {
  
  distmat.list<-c()
  for (i in 1:length(list)) {
    df.i<-list[[i]]
    df.i.pivot<-dcast(df.i, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)
    rownames(df.i.pivot)<-df.i.pivot$marker_OTU
    df.i.pivot<-df.i.pivot[,-1]
    df.i.pivot.t<-t(df.i.pivot)
    notnull<-df.i.pivot.t[,colSums(df.i.pivot.t) !=0]
    notnull2<-notnull[rowSums(notnull) !=0,]
    percentile<-quantile(rowSums(notnull2), prob=0.15)
    set.seed(1234)
    rarefied<-rrarefy(notnull2,sample=percentile)
    distmat<-vegdist(rarefied, "bray", binary=TRUE)
    distmat.list[[i]]<-distmat
  }
  return(distmat.list)
}

########################################################
##### NEW FUNCTION TO ADD SITES TO METADATA AUTOMATICALLY
########################################################

# for each dist matrix, create groupings for mrpp, start by adding site column

add_sites<- function (distmat.list) {
  
  meta.df.list<-c()
  for (i in 1:length(distmat.list)) {
    distmat<-as.matrix(distmat.list[[i]])
    meta.df<-data.frame(names=rownames(distmat))
    
    for (j in 1:length(rownames(meta.df))) {
      
      if (grepl("^ILC", meta.df$names[j])) {
        meta.df$site[j]<-"ILC45"
      }
      else  {
        meta.df$site[j]<-"NZC85"
      }
    }
    rownames(meta.df)<-meta.df$names
    meta.df.list[[i]]<-meta.df
  }
  return(meta.df.list)
}

########################################################
##### NEW FUNCTION TO ADD LAYERS TO METADATA AUTOMATICALLY
########################################################

# for each dist matrix, create groupings for mrpp, add layer column

add_layers<- function (distmat.list, meta.list) {
  
  meta.df.list<-c()
  for (i in 1:length(distmat.list)) {
    distmat<-as.matrix(distmat.list[[i]])
    dist.df<-data.frame(names=rownames(distmat))
    metamat<-meta.list[[i]]
    
    for (j in 1:length(rownames(dist.df))) {
      
      if (grepl("B$", dist.df$names[j])) {
        metamat$layer[j]<-"Bryophyte"
      }
      else if (grepl("O$", dist.df$names[j])) {
        metamat$layer[j]<-"Organic"
      }
      else {
        metamat$layer[j]<-"Mineral"
      }
    }
    meta.df.list[[i]]<-metamat
  }
  return(meta.df.list)
}

########################################################
##### NEW FUNCTION TO GET WITHIN AND BETWEEN GROUP MEANS
########################################################

# Write function to get within and between group mean Sorensen distances for each rep

get_meandist <- function (distmat.list, meta.list) {
  
  group.df<-as.data.frame(matrix(nrow=length(distmat.list), ncol=6))
  x<-c("group1","within1","between1",
       "group2","within2","between2")
  colnames(group.df)<-x
  
  for (i in 1:length(distmat.list)) {
    distmat<-as.dist(distmat.list[[i]])
    metadf<-meta.list[[i]]

    # calc meandist for sites
    group.md <- meandist(distmat, metadf$site)
    group.md.summary<-summary(group.md)
    
    # get within and between group means, add to table
    group.df$group1[i]<-"site"
    group.df$within1[i]<-group.md.summary$W
    group.df$between1[i]<-group.md.summary$B
    
    # calc meandist for layers
    group.md <- meandist(distmat, metadf$layer)
    group.md.summary<-summary(group.md)
    
    # get within and between group means, add to table
    group.df$group2[i]<-"layer"
    group.df$within2[i]<-group.md.summary$W
    group.df$between2[i]<-group.md.summary$B
    
  }
  
  return(group.df)
  
}

###################################################################
# Read in metadata, edit GRDI format to remove excess columns because they mess up nrows
###################################################################

ENV<-read.table(file="metadata.csv", sep=",",head=TRUE)

# Move grdiname to rownames
rownames(ENV)<-ENV$grdiname

###################################################################

# Read infile prepared by python script
## LV2016_2.csv is the clean taxonomic assignment table
A<-read.table(file="LV2016_2.csv", head=TRUE, sep=",")

# Select phylum Arthropoda only
Arth_df<-A[A$phylum=="Arthropoda",]

# Create matrix of just 1C3E samples
B<-Arth_df[(Arth_df$cores==1 & Arth_df$extractions==3),]

# Create matrix of just 1C3E and XC3E samples
# Grab 1C3E
C1<-Arth_df[Arth_df$cores==1 & Arth_df$extractions==3,]

# Randomly subsample 8 unique gridcoord without replacement so 1C3E can be subsampled to create balanced design
C1.gridcoord<-sample(unique(C1$gridcoord), 4, replace=FALSE)

# Retrieve gridcoord
C1<-C1[C1$gridcoord=="11" |
         C1$gridcoord=="16" |
         C1$gridcoord=="33" |
         C1$gridcoord=="56",]

# Grab XC3E
C2<-Arth_df[Arth_df$cores!=1 & Arth_df$extractions==3,]

# combine C1 with C2
C3<-rbind(C1,C2)

# Create matrix of just 1C1E and 1C3E samples
# Grab 1C3E
D1<-Arth_df[Arth_df$cores==1 & Arth_df$extractions==3,]

# Match 8 gridcoord from 1C3E to match those in 1C1E for balanced design
D1<-D1[D1$gridcoord=="11" |
         D1$gridcoord=="15" |
         D1$gridcoord=="23" |
         D1$gridcoord=="31" |
         D1$gridcoord=="35" |
         D1$gridcoord=="43" |
         D1$gridcoord=="51" |
         D1$gridcoord=="55",]

# Grab 1C1E
D2<-Arth_df[Arth_df$cores==1 & Arth_df$extractions==1,]

# combine D1 with D2
D3<-rbind(D1,D2)

############################################
#### Create Fig 4A plot, just 1C3E samples
############################################

# Pivot to make matrix for vegan
matrix<-dcast(B, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)

# Move marker_OTU to row names
rownames(matrix)<-matrix$marker_OTU
matrix<-matrix[,-1]

# Transpose to get sites in rows, ESVs in columns
matrix2<-t(matrix)

# Remove columns with only zeros
notnull<-matrix2[,colSums(matrix2) !=0]

# Remove rows with only zeros
notnull2<-notnull[rowSums(notnull) !=0,]

# Calculate 15th percentile for rrarefy function
percentile<-quantile(rowSums(notnull2), prob=0.15)

# Set random seed for rarefaction
set.seed(1234)

# Rarefy the dataset down to the 15th percentile
df<-rrarefy(notnull2,sample=percentile)

# Convert to presence-absence matrix
df[df>0] <-1

# Scree plots to determine number of dimensions to use for NMDS

#pdf("scree.pdf")
#nmds.scree(df)
#dev.off()

# Do 2 dimensional NMDS, no environment file for now 
nmds2<-metaMDS(df, k=2,trymax=100)
# stress = 0.146

# Create grouping matrix for samples by grabbing row names from above matrix
sample_df<-data.frame(row.names(df))

# Rename the column
names(sample_df)<-"GRDIname"

# Copy column to row names
row.names(sample_df)<-sample_df$GRDIname

# Split first column into their own fields
sample_df[,2:5]<-do.call('rbind', strsplit(as.character(sample_df$GRDIname),'_',fixed=TRUE))

# Remove first column
sample_df<-sample_df[,-1]

# Rename columns
names(sample_df)<-c("site","date","expt","rep_layer")

# Grab just layer and save to own column
sample_df$layer <- substr(sample_df$rep_layer, 3, 3)

# Grab sites & scores from NMDS output
B_site.sc <- data.frame(scores(nmds2, display = "sites"))

# Put it all in one df for ggplot
mergedB <- merge(B_site.sc,sample_df,by="row.names")

colnames(mergedB)[colnames(mergedB)=="Row.names"] <- "GRDIname"

# Create factors
mergedB$layer<-factor(mergedB$layer, levels=c("B","O","M"))
mergedB$expt<-factor(mergedB$expt, levels=c("1C3E"), labels=c("1"))
mergedB$site<-factor(mergedB$site, levels=c("ILC45","NZC85"))

# Compile coord for convex hulls
chulls <- ddply(mergedB, .(site), function(mergedB) mergedB[chull(mergedB$NMDS1, mergedB$NMDS2), ])

# Create scatter plot
# p1 legends off
p1<-ggplot(data=mergedB, aes(x=NMDS1, y=NMDS2)) + 
  ggtitle("A)") +
  geom_polygon(data=chulls, aes(x=NMDS1, y=NMDS2, fill=site), alpha=0.5) +
  geom_point(data=mergedB, aes(x=NMDS1, y=NMDS2, color=layer), pch=0) +
  scale_fill_manual(values=c("black","grey")) +
  scale_color_manual(values=c("#4DAF4A","#377EB8","#FF7F00")) +
  theme_bw() +
  theme(
         plot.title = element_text(size = 10),
         axis.text.x = element_text(hjust=1),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         legend.key=element_blank(),
        axis.title = element_text(size=10),
        axis.text = element_text(size=10),
        legend.title = element_blank(),
        legend.text = element_text(size=8))+
  guides(fill=FALSE,
         color=FALSE)

#p1b legends on
p1b<-ggplot(data=mergedB, aes(x=NMDS1, y=NMDS2)) + 
  ggtitle("A)") +
  geom_polygon(data=chulls, aes(x=NMDS1, y=NMDS2, fill=site), alpha=0.5) +
  geom_point(data=mergedB, aes(x=NMDS1, y=NMDS2, color=layer), pch=0) +
  scale_fill_manual(values=c("black","grey")) +
  scale_color_manual(values=c("#4DAF4A","#377EB8","#FF7F00")) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10),
    axis.text.x = element_text(hjust=1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(),
    axis.title = element_text(size=10),
    axis.text = element_text(size=10),
    legend.title = element_blank(),
    legend.text = element_text(size=8))

# Process metadata
# Only keep the samples that are in df
ENV2<-subset(ENV,rownames(ENV) %in% rownames(df) )

# Create factors
ENV2$layer<-factor(ENV2$layer, levels=c("Moss","Organic","Mineral"), labels=c("Bryophyte","Organic","Mineral"))

# Assess dispersion (variance) using ANOVA
# Create distance matrix based on P-A data using jaccard dissimilarity
sor<-vegdist(df, "bray", binary=TRUE)

# Calculate beta dispersion (homogeneity needed for adonis)
#site
bd_site<-betadisper(sor, factor(ENV2[,2]))
#layer
bd_layer<-betadisper(sor, factor(ENV2[,4]))

# check for heterogeneity of beta dispersions within groups
anova(bd_site) # n/s
anova(bd_layer) # n/s

pdf("BetaDispersion_boxplot_A.pdf")
par(mfrow=c(2,2))
boxplot(bd_site)
mtext("A) Site", side=3, adj=0, line=1.2)
boxplot(bd_layer)
mtext("B) Layer", side=3, adj=0, line=1.2)
dev.off()

# Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
pdf("stressplot_A.pdf")
stressplot(nmds2)
gof <-goodness(nmds2)
gof
plot(nmds2, display = "sites", type="n")
points(nmds2, display="sites",cex=2*gof/mean(gof))
dev.off()

# Use ADONIS to test significance of groupings (esp. sites, layers, don't expect any diff for expt)

# 1. Test the interaction among all groups first to see how to proceed (use strata so randomizations occur within each site)
adonis(sor~site*layer, data=ENV2, permutations=999, strata=ENV2$site)
#site R2 = 0.08243, pvalue = 0.038 *
#layer n/a
#site:layer n/s
# No significant interactions so go ahead and look at main effects (only site is significant)

# make a community table by merging df with ENV2
comm<-merge(df, ENV2, by="row.names")

# cite pairwiseAdonis
# Martinez Arbizu, P. (2017). pairwiseAdonis: Pairwise multilevel comparison using adonis. R package version 0.0.1.

# do pairwise tests using custom function above with strata
pairwise.adonis(comm[,2:2536], paste(comm$site))
#        pairs      F.Model     R2     p.value   p.adjusted sig
# 1 ILC45 vs NZC85 13.54945 0.06088288   0.001      0.001  **

pairwise.adonis(comm[,2:2536], paste(comm$layer))
#            pairs        F.Model     R2       p.value  p.adjusted sig
# 1 Bryophyte vs Mineral 18.169333 0.11415096   0.001      0.003   *
#   2 Bryophyte vs Organic  9.613643 0.06557127   0.001      0.003   *
#   3   Mineral vs Organic 12.722034 0.08440726   0.001      0.003   *

###################################################
# Look at average sorensen dissimilarity within and among groups for 1, 2, 4, 6, 8, 15 randomly sampled cores
## bioinformatic pooling
###################################################

# randomly sample 1 gridcoord, repeat this sampling 4 times for 1C3E
B.gridcoord<-t(as.matrix(sapply(rep(1,4), function (x) sample(unique(B$gridcoord), x, replace=FALSE))))

# Retrieve target gridcoord rows for each of 4 samples
#B.samples<-as.data.frame(sapply(B.gridcoord, function(x) B[B$gridcoord==x,]))
B.samples<-get_target_gridcoord(B.gridcoord,B)

# get distance matrices
distmat.list.1<-get_distmat(B.samples)

# add sites to metadata for each sample
meta.list.1<-add_sites(distmat.list.1)

# add layers to metadata for each sample
meta.list.1<-add_layers(distmat.list.1,meta.list.1)

# get within and between group mean Sorensen dissimilarities
group.1<-get_meandist(distmat.list.1, meta.list.1)

#move rownames to new col sample
group.1<-setDT(group.1, keep.rownames = TRUE)[]

# reshape df for melt
top<-group.1[,1:4]
names<-c("rn","group","within","between")
names(top)<-names
bottom<-group.1[,c(1,5:7)]
names(bottom)<-names
group.1<-rbind(top,bottom)

# melt for ggplot
group.1<-melt(group.1, id.var=c("rn","group"), 
     variable.name=c("disttype"),
     value.name="meandist")

# add column to indicate pooled cores
group.1$cores<-"1"

# randomly sample 2 gridcoord, repeat this sampling 4 times for 1C3E
B.gridcoord<-sapply(rep(2,4), function (x) sample(unique(B$gridcoord), x, replace=FALSE))

# Retrieve target gridcoord rows for each of 4 samples
B.samples<-get_target_gridcoord(B.gridcoord,B)

# get distance matrices
distmat.list.2<-get_distmat(B.samples)

# add sites to metadata for each sample
meta.list.2<-add_sites(distmat.list.2)

# add layers to metadata for each sample
meta.list.2<-add_layers(distmat.list.2,meta.list.2)

# get within and between group mean Sorensen dissimilarities
group.2<-get_meandist(distmat.list.2, meta.list.2)

#move rownames to new col sample
group.2<-setDT(group.2, keep.rownames = TRUE)[]

# reshape df for melt
top<-group.2[,1:4]
names<-c("rn","group","within","between")
names(top)<-names
bottom<-group.2[,c(1,5:7)]
names(bottom)<-names
group.2<-rbind(top,bottom)

# melt for ggplot
group.2<-melt(group.2, id.var=c("rn","group"), 
              variable.name=c("disttype"),
              value.name="meandist")

# add column to indicate pooled cores
group.2$cores<-"2"

# randomly sample 4 gridcoord, repeat this sampling 4 times for 1C3E
B.gridcoord<-sapply(rep(4,4), function (x) sample(unique(B$gridcoord), x, replace=FALSE))

# Retrieve target gridcoord rows for each of 4 samples
B.samples<-get_target_gridcoord(B.gridcoord,B)

# get distance matrices
distmat.list.4<-get_distmat(B.samples)

# add sites to metadata for each sample
meta.list.4<-add_sites(distmat.list.4)

# add layers to metadata for each sample
meta.list.4<-add_layers(distmat.list.4,meta.list.4)

# get within and between group mean Sorensen dissimilarities
group.4<-get_meandist(distmat.list.4, meta.list.4)

#move rownames to new col sample
group.4<-setDT(group.4, keep.rownames = TRUE)[]

# reshape df for melt
top<-group.4[,1:4]
names<-c("rn","group","within","between")
names(top)<-names
bottom<-group.4[,c(1,5:7)]
names(bottom)<-names
group.4<-rbind(top,bottom)

# melt for ggplot
group.4<-melt(group.4, id.var=c("rn","group"), 
              variable.name=c("disttype"),
              value.name="meandist")

# add column to indicate pooled cores
group.4$cores<-"4"

# randomly sample 6 gridcoord, repeat this sampling 4 times for 1C3E
B.gridcoord<-sapply(rep(6,4), function (x) sample(unique(B$gridcoord), x, replace=FALSE))

# Retrieve target gridcoord rows for each of 4 samples
B.samples<-get_target_gridcoord(B.gridcoord,B)

# get distance matrices
distmat.list.6<-get_distmat(B.samples)

# add sites to metadata for each sample
meta.list.6<-add_sites(distmat.list.6)

# add layers to metadata for each sample
meta.list.6<-add_layers(distmat.list.6,meta.list.6)

# get within and between group mean Sorensen dissimilarities
group.6<-get_meandist(distmat.list.6, meta.list.6)

#move rownames to new col sample
group.6<-setDT(group.6, keep.rownames = TRUE)[]

# reshape df for melt
top<-group.6[,1:4]
names<-c("rn","group","within","between")
names(top)<-names
bottom<-group.6[,c(1,5:7)]
names(bottom)<-names
group.6<-rbind(top,bottom)

# melt for ggplot
group.6<-melt(group.6, id.var=c("rn","group"), 
              variable.name=c("disttype"),
              value.name="meandist")

# add column to indicate pooled cores
group.6$cores<-"6"

# randomly sample 8 gridcoord, repeat this sampling 4 times for 1C3E
B.gridcoord<-sapply(rep(8,4), function (x) sample(unique(B$gridcoord), x, replace=FALSE))

# Retrieve target gridcoord rows for each of 4 samples
B.samples<-get_target_gridcoord(B.gridcoord,B)

# get distance matrices
distmat.list.8<-get_distmat(B.samples)

# add sites to metadata for each sample
meta.list.8<-add_sites(distmat.list.8)

# add layers to metadata for each sample
meta.list.8<-add_layers(distmat.list.8,meta.list.8)

# get within and between group mean Sorensen dissimilarities
group.8<-get_meandist(distmat.list.8, meta.list.8)

#move rownames to new col sample
group.8<-setDT(group.8, keep.rownames = TRUE)[]

# reshape df for melt
top<-group.8[,1:4]
names<-c("rn","group","within","between")
names(top)<-names
bottom<-group.8[,c(1,5:7)]
names(bottom)<-names
group.8<-rbind(top,bottom)

# melt for ggplot
group.8<-melt(group.8, id.var=c("rn","group"), 
              variable.name=c("disttype"),
              value.name="meandist")

# add column to indicate pooled cores
group.8$cores<-"8"

# randomly sample 15 gridcoord, repeat this sampling 4 times for 1C3E
B.gridcoord<-sapply(rep(15,4), function (x) sample(unique(B$gridcoord), x, replace=FALSE))

# Retrieve target gridcoord rows for each of 4 samples
B.samples<-get_target_gridcoord(B.gridcoord,B)

# get distance matrices
distmat.list.15<-get_distmat(B.samples)

# add sites to metadata for each sample
meta.list.15<-add_sites(distmat.list.15)

# add layers to metadata for each sample
meta.list.15<-add_layers(distmat.list.15,meta.list.15)

# get within and between group mean Sorensen dissimilarities
group.15<-get_meandist(distmat.list.15, meta.list.15)

#move rownames to new col sample
group.15<-setDT(group.15, keep.rownames = TRUE)[]

# reshape df for melt
top<-group.15[,1:4]
names<-c("rn","group","within","between")
names(top)<-names
bottom<-group.15[,c(1,5:7)]
names(bottom)<-names
group.15<-rbind(top,bottom)

# melt for ggplot
group.15<-melt(group.15, id.var=c("rn","group"), 
              variable.name=c("disttype"),
              value.name="meandist")

# add column to indicate pooled cores
group.15$cores<-"15"

# combine meandist for 1-15 pooled cores
group.combined<-rbind(group.1, group.2, group.4, group.6, group.8, group.15)

# create factors
group.combined$cores<-factor(group.combined$cores, levels=c("1","2","4","6","8","15"))
group.combined$group<-factor(group.combined$group, levels=c("site","layer"),
                             labels=c("Sites","Strata"))
group.combined$disttype<-factor(group.combined$disttype, levels=c("within","between"),
                                labels=c("Within","Between"))

# boxplot cores vs meandist for within & between series
# box legends off
box<-ggplot(group.combined, aes(x=cores, y=meandist, fill=disttype)) +
  ggtitle("A)") +
  facet_wrap(~group) +
  geom_boxplot() +
  labs(x="Bioinformatically Pooled cores",y="Average Sorensen dissimilarity") +
  theme_bw() +
  theme(
        plot.title = element_text(size = 10),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none")

# boxb legends on
boxb<-ggplot(group.combined, aes(x=cores, y=meandist, fill=disttype)) +
  ggtitle("A)") +
  facet_wrap(~group) +
  geom_boxplot() +
  labs(x="Bioinformatically pooled samples",y="Average Sorensen dissimilarity") +
  theme_bw() +
  theme(
    text=element_text(size=12),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.position = "right") +
  guides(fill=guide_legend(title="Groups"))

# Check for normality with Shapiro-Wilk’s test, sig result means not normal
shapiro.test(group.combined$meandist)
# W = 0.82419, p-value = 2.517e-09

#visual inspection, sometimes small sample sizes can pass normality tests
qqPlot(group.combined$meandist)

# Use Kruskal-Wallis test to check for any significant differences among groups
# Do just for Between Sites
group.new<-group.combined[group.combined$disttype=="Between" & group.combined$group=="Sites",]
kruskal.test(meandist ~ cores, data = group.new)
# pval = 0.04383
# Do just for Between Strata
group.new2<-group.combined[group.combined$disttype=="Between" & group.combined$group=="Strata",]
kruskal.test(meandist ~ cores, data = group.new2)
# pval = 0.01619

# Use multiple pairwise-comparison bewteen groups to check for specific diffs across cores 
#(for Between, separately for Sites & Strata)
# p.adjust method Benjamini & Hochberg (1995)
pairwise.wilcox.test(group.new$meandist, group.new$cores,
                     p.adjust.method = "BH")
# all n/s after p.adjust for multiple comparisons
pairwise.wilcox.test(group.new2$meandist, group.new2$cores,
                     p.adjust.method = "BH")
# all n/s after p.adjust for multiple comparisons

################################################
#### Create Fig 4B plot, 1C3E and XC3E samples
################################################

# Pivot to make matrix for vegan
matrix<-dcast(C3, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)

# Move marker_OTU to row names
rownames(matrix)<-matrix$marker_OTU
matrix<-matrix[,-1]

# Transpose to get sites in rows, ESVs in columns
matrix2<-t(matrix)

# Remove columns with only zeros
notnull<-matrix2[,colSums(matrix2) !=0]

# Remove rows with only zeros
notnull2<-notnull[rowSums(notnull) !=0,]

# Calculate 15th percentile for rrarefy function
percentile<-quantile(rowSums(notnull2), prob=0.15)

# Set random seed
set.seed(1234)

# Rarefy the dataset down to the 15th percentile
df<-rrarefy(notnull2,sample=percentile)

# Convert to presence-absence matrix
df[df>0] <-1

##### Do scree plots
#Scree plots to determine number of dimensions to use for NMDS

#pdf("scree.pdf")
#nmds.scree(df)
#dev.off()

### Do NMDS
#Do 2 dimensional NMDS, no environment file for now 
nmds2<-metaMDS(df, k=2,trymax=100)
# stress = 0.108

# Create grouping matrix for samples
# Grab row names from matrix above
sample_df<-data.frame(row.names(df))

# Rename the column
names(sample_df)<-"GRDIname"

# Copy column to row names
row.names(sample_df)<-sample_df$GRDIname

# Split first column into their own fields
sample_df[,2:5]<-do.call('rbind', strsplit(as.character(sample_df$GRDIname),'_',fixed=TRUE))

# Remove first column
sample_df<-sample_df[,-1]

# Rename columns
names(sample_df)<-c("site","date","expt","rep_layer")

# Grab just layer and save to own column
sample_df$layer <- substr(sample_df$rep_layer, 3, 3)

# Recode 12/13/14/15C3E with 1215C3E for ILC and NZC to enable compariability
sample_df$expt <- gsub("15C3E", "915C3E", sample_df$expt)
sample_df$expt <- gsub("9C3E", "915C3E", sample_df$expt)
sample_df$expt <- gsub("12C3E", "915C3E", sample_df$expt)
sample_df$expt <- gsub("13C3E", "915C3E", sample_df$expt)
sample_df$expt <- gsub("14C3E", "915C3E", sample_df$expt)

# Grab sites & scores
C_site.sc <- data.frame(scores(nmds2, display = "sites"))

# Put it all in one df
mergedC <- merge(C_site.sc,sample_df,by="row.names")
colnames(mergedC)[colnames(mergedC)=="Row.names"] <- "GRDIname"

# Create factors
mergedC$layer<-factor(mergedC$layer, levels=c("B","O","M"), 
                                   labels=c("Bryophyte","Organic","Mineral"))
mergedC$expt<-factor(mergedC$expt, levels=c("1C3E","2C3E","4C3E","6C3E","8C3E","915C3E"), 
                                 labels=c("1","2","4","6","8","9-15"))
mergedC$site<-factor(mergedC$site, levels=c("ILC45","NZC85"), 
                                 labels=c("Island Lake", "Nimitz"))

# Create convex hulls
chulls <- ddply(mergedC, .(site), function(mergedC) mergedC[chull(mergedC$NMDS1, mergedC$NMDS2), ])

# Create scatter plot
# p2 legends off
p2<-ggplot(data=mergedC, aes(x=NMDS1, y=NMDS2)) + 
  ggtitle("B)") +
  geom_polygon(data=chulls, aes(x=NMDS1, y=NMDS2, fill=site), alpha=0.5) +
  geom_point(data=mergedC, aes(x=NMDS1, y=NMDS2, pch=expt, color=layer)) +
  scale_fill_manual(name="Site",values=c("black","grey")) +
  scale_color_manual(values=c("#4DAF4A","#377EB8","#FF7F00")) +
  scale_shape_manual(values=c(0,1,2,3,4,5)) +
  theme_bw() +
  theme(
        text = element_text(size=12),
        axis.text.x = element_text(hjust=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        legend.title = element_blank())+
  guides(
        fill = FALSE, 
        colour = FALSE,
        shape = FALSE)

#p2b legends on
p2b<-ggplot(data=mergedC, aes(x=NMDS1, y=NMDS2)) + 
  ggtitle("B)") +
  geom_polygon(data=chulls, aes(x=NMDS1, y=NMDS2, fill=site), alpha=0.5) +
  geom_point(data=mergedC, aes(x=NMDS1, y=NMDS2, pch=expt, color=layer)) +
  scale_fill_manual(name="Site", values=c("black","grey")) +
  scale_color_manual(name="Strata", values=c("#4DAF4A","#377EB8","#FF7F00")) +
  scale_shape_manual(name="Pooled Cores", values=c(0,1,2,3,4,5)) +
  theme_bw() +
  theme(
    text = element_text(size=12),
    axis.text.x = element_text(hjust=1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank()) +
  guides(fill = guide_legend(order=1),
         color = guide_legend(order=2),
         shape = guide_legend(order=3))

# Process metadata
# Only keep the samples that are in df
ENV2<-subset(ENV,rownames(ENV) %in% rownames(df) )

# Recode ENV2 expt 
ENV2$expt <- gsub("15C3E", "915C3E", ENV2$expt)
ENV2$expt <- gsub("9C3E", "915C3E", ENV2$expt)
ENV2$expt <- gsub("12C3E", "915C3E", ENV2$expt)
ENV2$expt <- gsub("13C3E", "915C3E", ENV2$expt)
ENV2$expt <- gsub("14C3E", "915C3E", ENV2$expt)

# Create factors
ENV2$expt<-factor(ENV2$expt, 
                 levels=c("1C3E","2C3E","4C3E","6C3E","8C3E","915C3E"), 
                 labels=c("1","2","4","6","8","9-15"))
ENV2$layer<-factor(ENV2$layer, levels=c("Moss","Organic","Mineral"), labels=c("Bryophyte","Organic","Mineral"))

# Create distance matrix based on P-A data using jaccard dissimilarity
sor<-vegdist(df, "bray", binary=TRUE)

# Calculate beta dispersion (homogeneity needed for adonis)
#expt
bd_expt<-betadisper(sor, factor(ENV2[,3]))
#site
bd_site<-betadisper(sor, factor(ENV2[,2]))
#layer
bd_layer<-betadisper(sor, factor(ENV2[,4]))

# check for homogeneity of beta dispersions within groups
anova(bd_expt) # n/s *removed heterogenous beta dispersions by creating balanced design
anova(bd_site) # n/s
anova(bd_layer) # 0.01167 * heterogeous beta dispersions across layers now

#PERMANOVA (which is basically adonis()) was found to be largely unaffected by heterogeneity in Anderson & Walsh's simulations but only for balanced designs.
#Anderson MJ, Walsh DCI. PERMANOVA, ANOSIM, and the Mantel test in the face of heterogeneous dispersions: What null hypothesis are you testing? Ecological monographs [Internet] 2013; 83: 557. Available from: http://doi.org/10.1890/12-2010.1

pdf("BetaDispersion_boxplot_B.pdf")
par(mfrow=c(2,2))
boxplot(bd_expt, las=2)
mtext("A) Experiment", side=3, adj=0, line=1.2)
boxplot(bd_site)
mtext("B) Site", side=3, adj=0, line=1.2)
boxplot(bd_layer)
mtext("C) Layer", side=3, adj=0, line=1.2)
dev.off()

##### Shephards curve and goodness of fit calcs
# Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
pdf("stressplot_B.pdf")
stressplot(nmds2)
gof <-goodness(nmds2)
gof
plot(nmds2, display = "sites", type="n")
points(nmds2, display="sites",cex=2*gof/mean(gof))
dev.off()

# Use ADONIS to test significance of groupings (esp. sites, layers, don't expect any diff for expt)
# Test significance of expt groupings (use strata so randomizations occur within each site)
# don't bother including expt because of significant hetergeneity of beta dispersions

# 1. Test the interaction first to see how to proceed
adonis(sor~site*layer*expt, data=ENV2, permutations=999, strata=ENV2$site)
#site n/s
#layer R2 = 0.06954, pvalue = 0.001 *** <- sig. beta disper detected within layers but balanced design so permanova largely unaffected by heterogeneity
#expt n/s
#site:layer n/s
#site:expt n/s
#layer:expt n/s
#site:layer:expt n/s
# No significant interactive effects, go ahead and go with main tests

#https://stats.stackexchange.com/questions/314184/betadisper-is-significant-can-i-still-do-adonis-and-anosim
#PERMANOVA (which is basically adonis()) was found to be largely unaffected by heterogeneity in Anderson & Walsh's simulations but only for balanced designs.
#Anderson MJ, Walsh DCI. PERMANOVA, ANOSIM, and the Mantel test in the face of heterogeneous dispersions: What null hypothesis are you testing? Ecological monographs [Internet] 2013; 83: 557. Available from: http://doi.org/10.1890/12-2010.1

#2. do pairwise comparisons, to check for sig differences across Bryophyte and Organic layers across sites
# make a community table by merging df with ENV2
comm<-merge(df, ENV2, by="row.names")

# cite pairwiseAdonis
# Martinez Arbizu, P. (2017). pairwiseAdonis: Pairwise multilevel comparison using adonis. R package version 0.0.1.

# do pairwise tests using custom function above with strata
pairwise.adonis(comm[,2:2548], paste(comm$site))
# pairs  F.Model         R2 p.value p.adjusted sig
# 1 ILC45 vs NZC85 11.81701 0.07682514   0.001      0.001  **
pairwise.adonis(comm[,2:2548], paste(comm$layer))
# pairs   F.Model         R2 p.value p.adjusted sig
# 1 Bryophyte vs Mineral 19.790419 0.17391990   0.001      0.003   *
#   2 Bryophyte vs Organic  9.469663 0.09152116   0.001      0.003   *
#   3   Mineral vs Organic 15.280058 0.13982476   0.001      0.003   *
pairwise.adonis(comm[,2:2548], paste(comm$expt))
# pairs   F.Model         R2 p.value p.adjusted sig
# 1  9-15 vs 1 1.1418094 0.02422074   0.223          1    
# 2  9-15 vs 2 0.9274855 0.01976423   0.507          1    
# 3  9-15 vs 4 0.7829594 0.01673600   0.781          1    
# 4  9-15 vs 6 1.0673009 0.02267606   0.320          1    
# 5  9-15 vs 8 0.5780683 0.01241074   0.988          1    
# 6     1 vs 2 1.0174486 0.02163981   0.383          1    
# 7     1 vs 4 1.1483435 0.02435597   0.216          1    
# 8     1 vs 6 1.0753451 0.02284306   0.334          1    
# 9     1 vs 8 1.0313413 0.02192881   0.367          1    
# 10    2 vs 4 0.9174287 0.01955411   0.542          1    
# 11    2 vs 6 1.1547728 0.02448899   0.246          1    
# 12    2 vs 8 0.8584848 0.01832080   0.631          1    
# 13    4 vs 6 0.8456669 0.01805219   0.697          1    
# 14    4 vs 8 0.7077132 0.01515196   0.915          1    
# 15    6 vs 8 0.7964822 0.01702013   0.786          1 

#getOption("max.print")
options(max.print=999999)

#print results to file
# sink("Fig4_B_allPairwise.txt")
# print(x)
# sink()

###########################################
#### Create Fig 4C plot, just 1C1E & 1C3E samples (subsampled down for balanced design)
###########################################

# Pivot to make matrix for vegan
matrix<-dcast(D3, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)

# Move marker_OTU to row names
rownames(matrix)<-matrix$marker_OTU
matrix<-matrix[,-1]

# Transpose to get sites in rows, ESVs in columns
matrix2<-t(matrix)

# Remove columns with only zeros
notnull<-matrix2[,colSums(matrix2) !=0]

# Remove rows with only zeros
notnull2<-notnull[rowSums(notnull) !=0,]

# Calculate 15th percentile for rrarefy function
percentile<-quantile(rowSums(notnull2), prob=0.15)

set.seed(1234)

# Rarefy the dataset down to the 15th percentile
df<-rrarefy(notnull2,sample=percentile)

# Convert to presence-absence matrix
df[df>0] <-1

##### Do scree plots

#Scree plots to determine number of dimensions to use for NMDS

#pdf("scree.pdf")
#nmds.scree(df)
#dev.off()

### Do NMDS

#Do 2 dimensional NMDS, no environment file for now 
nmds2<-metaMDS(df, k=2,trymax=100)
# stress = 0.123

# Create grouping matrix for samples
## Grab row names from matrix above
sample_df<-data.frame(row.names(df))
## Rename the column
names(sample_df)<-"GRDIname"
## Copy column to row names
row.names(sample_df)<-sample_df$GRDIname
## Split first column into their own fields
sample_df[,2:5]<-do.call('rbind', strsplit(as.character(sample_df$GRDIname),'_',fixed=TRUE))
## Remove first column
sample_df<-sample_df[,-1]
## Rename columns
names(sample_df)<-c("site","date","expt","rep_layer")
## Grab just layer and save to own column
sample_df$layer <- substr(sample_df$rep_layer, 3, 3)

# Grab sites & scores
D_site.sc <- data.frame(scores(nmds2, display = "sites"))

# Put it all in one df
mergedD <- merge(D_site.sc,sample_df,by="row.names")
colnames(mergedD)[colnames(mergedD)=="Row.names"] <- "GRDIname"

#create factors and levels for ILC
mergedD$layer<-factor(mergedD$layer, levels=c("B","O","M"), 
                     labels=c("Bryophyte","Organic","Mineral"))
mergedD$expt<-factor(mergedD$expt, levels=c("1C1E","1C3E"), 
                    labels=c("1","3"))
mergedD$site<-factor(mergedD$site, levels=c("ILC45","NZC85"),
                    labels=c("Island Lake","Nimitz"))

# Create convex hulls
chulls <- ddply(mergedD, .(site), function(mergedD) mergedD[chull(mergedD$NMDS1, mergedD$NMDS2), ])

# Create scatter plot
# p3 legends off
p3<-ggplot(data=mergedD, aes(x=NMDS1, y=NMDS2)) + 
  ggtitle("C)") +
  geom_polygon(data=chulls, aes(x=NMDS1, y=NMDS2, fill=site), alpha=0.5) +
  geom_point(data=mergedD, aes(x=NMDS1, y=NMDS2, pch=expt, color=layer)) +
  scale_fill_manual(values=c("black","grey")) +
  scale_color_manual(values=c("#4DAF4A","#377EB8","#FF7F00")) +
  scale_shape_manual(name="Extractions", values=c(6,8)) +
  theme_bw() +
  theme(
        text=element_text(size=12),
        axis.text.x = element_text(hjust=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        legend.title = element_blank(),
        legend.position="right",
        legend.spacing.y=unit(1,'cm')) +
  guides(fill = FALSE, 
         colour = FALSE,
         shape= FALSE)

# p3g legends on
p3b<-ggplot(data=mergedD, aes(x=NMDS1, y=NMDS2)) + 
  ggtitle("C)") +
  geom_polygon(data=chulls, aes(x=NMDS1, y=NMDS2, fill=site), alpha=0.5) +
  geom_point(data=mergedD, aes(x=NMDS1, y=NMDS2, pch=expt, color=layer)) +
  scale_fill_manual(name="Site", values=c("black","grey")) +
  scale_color_manual(name="Strata", values=c("#4DAF4A","#377EB8","#FF7F00")) +
  scale_shape_manual(name="DNA extractions", values=c(6,8)) +
  theme_bw() +
  theme(
    text = element_text(size=12),
    axis.text.x = element_text(hjust=1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank()) +
  guides(fill = guide_legend(order=1),
         color = guide_legend(order=2),
         shape = guide_legend(order=3))

lay<-rbind(c(1,1,1),
           c(2,3,4))

g<-grid.arrange(box, p2, p3, layout_matrix=lay)
ggsave("F3_NMDS.pdf", g, width = 8, height = 9, units = c("in"))

g2<-grid.arrange(boxb, p2b, p3b, layout_matrix=lay)
ggsave("F3_NMDSb.pdf", g2, width = 8, height = 9, units = c("in"))

# Prep metadata format
# Only keep the samples that are in df
ENV2<-subset(ENV,rownames(ENV) %in% rownames(df) )

# Create factors
ENV2$expt<-factor(ENV2$expt, 
                  levels=c("1C1E","1C3E"), labels=c("1C1E","1C3E"))
ENV2$layer<-factor(ENV2$layer, levels=c("Moss","Organic","Mineral"), labels=c("Bryophyte","Organic","Mineral"))

# Assess dispersion (variance) using ANOVA
# Create distance matrix based on P-A data using jaccard dissimilarity
sor<-vegdist(df, "bray", binary=TRUE)

# Calculate beta dispersion (homogeneity needed for adonis)
#expt
bd_expt<-betadisper(sor, factor(ENV2[,3]))
#site
bd_site<-betadisper(sor, factor(ENV2[,2]))
#layer
bd_layer<-betadisper(sor, factor(ENV2[,4]))

# check for homogeneity of beta dispersions within groups
anova(bd_expt) # n/s
anova(bd_site) # n/s
anova(bd_layer) # n/s

pdf("BetaDispersion_boxplot_C.pdf")
par(mfrow=c(2,2))
boxplot(bd_expt, las=2)
mtext("A) Experiment", side=3, adj=0, line=1.2)
boxplot(bd_site)
mtext("B) Site", side=3, adj=0, line=1.2)
boxplot(bd_layer)
mtext("C) Layer", side=3, adj=0, line=1.2)
dev.off()

##### Shephards curve and goodness of fit calcs
# Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
pdf("stressplot_C.pdf")
stressplot(nmds2)
gof <-goodness(nmds2)
gof
plot(nmds2, display = "sites", type="n")
points(nmds2, display="sites",cex=2*gof/mean(gof))
dev.off()

# Use ADONIS to test significance of groupings (esp. sites, layers, don't expect any diff for expt)
# Test significance of expt groupings (use strata so randomizations occur within each site)
# 1. Test the interaction first to see how to proceed
adonis(sor~site*layer*expt, data=ENV2, permutations=999, strata=ENV2$site)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# site             1     0.158 0.15756  0.4830 0.00457  0.001 ***
#   layer            2     1.809 0.90446  2.7723 0.05245  0.001 ***
#   expt             1     3.041 3.04057  9.3199 0.08817  0.001 ***
#   site:layer       2     0.418 0.20905  0.6408 0.01212  0.986    
# site:expt        1     0.145 0.14472  0.4436 0.00420  1.000    
# layer:expt       2     1.649 0.82435  2.5268 0.04781  0.002 ** 
#   site:layer:expt  2     0.516 0.25792  0.7906 0.01496  0.785    
# Residuals       82    26.752 0.32624         0.77573           
# Total           93    34.486                 1.00000           
# The significant interaction means that tests of main effects (individual factors alone) may not be meaningful

# 2. The next step is to do pair-wise comparisons for the factor of interest separately within each level of the other factor
# i.e. pair-wise comparisons of levels of factor A (sites) within each level of factor B (layers)

# make a community table by merging df with ENV2
comm<-merge(df, ENV2, by="row.names")

# cite pairwiseAdonis
# Martinez Arbizu, P. (2017). pairwiseAdonis: Pairwise multilevel comparison using adonis. R package version 0.0.1.

# do pairwise tests using custom function above with strata, just do layer and expt because no sig interactions
options(max.print=999999)
x<-pairwise.adonis(comm[,2:1627], paste(comm$site, comm$layer, comm$expt))
#all n/s
x2<-pairwise.adonis(comm[,2:1627], paste(comm$site))
# pairs  F.Model         R2 p.value p.adjusted sig
# 1 ILC45 vs NZC85 8.682022 0.08623209   0.001      0.001  **

w<-pairwise.adonis(comm[,2:1627], paste(comm$layer, comm$expt))
# pairs    F.Model         R2 p.value p.adjusted sig
# 1    Bryophyte 1C1E vs Mineral 1C1E 10.0505229 0.25094611   0.001      0.015   .
# 2    Bryophyte 1C1E vs Organic 1C1E  3.7530111 0.11458522   0.001      0.015   .
# 3  Bryophyte 1C1E vs Bryophyte 1C3E  0.4842159 0.01588415   0.988      1.000    
# 4    Bryophyte 1C1E vs Mineral 1C3E 12.9005674 0.30070855   0.001      0.015   .
# 5    Bryophyte 1C1E vs Organic 1C3E  3.1867902 0.09900926   0.002      0.030   .
# 6      Mineral 1C1E vs Organic 1C1E  5.7871269 0.16635829   0.001      0.015   .
# 7    Mineral 1C1E vs Bryophyte 1C3E  9.2322536 0.23532305   0.001      0.015   .
# 8      Mineral 1C1E vs Mineral 1C3E  0.7193820 0.02341785   0.931      1.000    
# 9      Mineral 1C1E vs Organic 1C3E  7.4086569 0.20348614   0.001      0.015   .
# 10   Organic 1C1E vs Bryophyte 1C3E  3.9413870 0.11964848   0.001      0.015   .
# 11     Organic 1C1E vs Mineral 1C3E  6.8269084 0.19055254   0.001      0.015   .
# 12     Organic 1C1E vs Organic 1C3E  0.5792374 0.02026777   0.963      1.000    
# 13   Bryophyte 1C3E vs Mineral 1C3E 11.9511528 0.28488258   0.001      0.015   .
# 14   Bryophyte 1C3E vs Organic 1C3E  3.3248578 0.10285762   0.001      0.015   .
# 15     Mineral 1C3E vs Organic 1C3E  8.8644077 0.23410924   0.001      0.015   .

# site by layer by expt ALL n/s

#sink("Fig4_C_alPairwise.txt")
#print(x)
#sink()
