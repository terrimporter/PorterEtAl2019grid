# Teresita M. Porter, July 8, 2019

library(vegan)
library(reshape2)
library(gridExtra)
library(grid)
library(ggplot2)
library(plyr)
library(data.table) # rowname to first col
library(goeveg) # scree

citation('goeveg')
citation('vegan')

########################################################
##### NEW FUNCTION TO GET TARGET GRIDCOORD DF's
########################################################

# write function to grab list of gridcoord from each column and retrieve appropriate rows

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
    out<-data.frame()
  }
  return(df.list)
}

###################################################################
# Read in metadata, edit GRDI format to remove excess columns because they mess up nrows
###################################################################

ENV<-read.table(file="metadata.csv", sep=",",head=TRUE)

# Move grdiname to rownames
rownames(ENV)<-ENV$grdiname

###################################################################

# Read infile
A<-read.table(file="LV2016_2.csv", head=TRUE, sep=",")

# Select phylum Arthropoda only
Arth_df<-A[A$phylum=="Arthropoda",]

# Create matrix of just 1C3E samples
B<-Arth_df[(Arth_df$cores==1 & Arth_df$extractions==3),]

# Create matrix of just 1C3E and XC3E samples
# Grab 1C3E
C1<-Arth_df[Arth_df$cores==1 & Arth_df$extractions==3,]
# Randomly subsample 4 unique gridcoord without replacement so can be compared with XC3E for balanced design
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
# Grab 1C1E
D2<-Arth_df[Arth_df$cores==1 & Arth_df$extractions==1,]

# Grab 1C3E
D1<-Arth_df[Arth_df$cores==1 & Arth_df$extractions==3,]
# Match 8 gridcoord from 1C1E from 1C3E also for direct comparison, balanced design
D1<-D1[D1$gridcoord=="11" |
         D1$gridcoord=="15" |
         D1$gridcoord=="23" |
         D1$gridcoord=="31" |
         D1$gridcoord=="35" |
         D1$gridcoord=="43" |
         D1$gridcoord=="51" |
         D1$gridcoord=="55",]

# combine D1 with D2
D3<-rbind(D1,D2)

############################################
#### Create Fig 4A plot, just 1C3E samples, all 36 samples
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
pdf("Scree_1C3E.pdf")
# check dims
dimcheckMDS(df)
dev.off()

# Do 2 dimensional NMDS, no environment file for now 
nmds2<-metaMDS(df, k=2,trymax=100)
# stress = 0.1291214
# linear fit R2 = 0.922, non-metric fit R2 = 0.983 (from stressplot below)

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
# Create distance matrix based on P-A data using binary bray curtis (Sorensen) dissimilarity
sor<-vegdist(df, "bray", binary=TRUE)

# Calculate beta dispersion (homogeneity needed for adonis)
# site
bd_site<-betadisper(sor, factor(ENV2[,2]))
# layer
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

# Shephards curve and goodness of fit calcs
# Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
pdf("stressplot_A.pdf")
stressplot(nmds2)
gof <-goodness(nmds2)
gof
plot(nmds2, display = "sites", type="n")
points(nmds2, display="sites",cex=2*gof/mean(gof))
dev.off()

# Use ADONIS to test for interactions among all groups first to see how to proceed (use strata so randomizations occur within each site)

# use strata so randomizations occurr within each site
set.seed(1234)
adonis(sor ~ site * layer, data=ENV2, permutations=999, strata = ENV2$site)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# site         1     6.539  6.5395 18.9193 0.08243  0.038 *
#   layer        2     0.975  0.4876  1.4108 0.01229  0.059 .
# site:layer   2     0.958  0.4792  1.3864 0.01208  0.055 .
# Residuals  205    70.858  0.3456         0.89319         
# Total      210    79.331                 1.00000         

################################################
#### Create Fig 4B plot, 1C3E and XC3E WITH OUTLIER REMOVED!!! THIS ONE, subsample 1C3E down to 4 for balanced design
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

# remove potential outlier, only NZC sample to cluster with ILC smaples
outlier<-"NZC85_20160701_6C3E_S3M"
notnull3<-notnull2[!rownames(notnull2) %in% outlier, ]

# Calculate 15th percentile for rrarefy function
percentile<-quantile(rowSums(notnull3), prob=0.15)

# Set random seed
set.seed(1234)

# Rarefy the dataset down to the 15th percentile
df<-rrarefy(notnull3, sample=percentile)

# Convert to presence-absence matrix
df[df>0] <-1

# Scree plots to determine number of dimensions to use for NMDS
pdf("Scree_XC3E_outlierremoved.pdf")
# check dims
dimcheckMDS(df)
dev.off()

### Do NMDS
# Do 2 dimensional NMDS 
nmds2<-metaMDS(df, k=2,trymax=100)
# stress = 0.106
# linear R2 = 0.989, non-metric fit = 0.95 (from stressplot below)
# NZC85_20160701_6C3E_S3M   1.74723771 -1.62779967 outlier

# Create grouping matrix for samples by grabbing row names from above matrix
sample_df<-data.frame(row.names(df))

# Copy column to row names
row.names(sample_df)<-sample_df[,1]

# rename column
names(sample_df)<-"GRDIname"

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
anova(bd_layer) # 0.007074 ** heterogeous beta dispersions across layers now

# NOTE: PERMANOVA (which is basically adonis()) was found to be largely unaffected by heterogeneity in Anderson & Walsh's simulations but only for balanced designs.
# Anderson MJ, Walsh DCI. PERMANOVA, ANOSIM, and the Mantel test in the face of heterogeneous dispersions: What null hypothesis are you testing? Ecological monographs [Internet] 2013; 83: 557. Available from: http://doi.org/10.1890/12-2010.1

pdf("BetaDispersion_boxplot_B_outlierremoved.pdf")
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
pdf("stressplot_B_outlierremoved.pdf")
stressplot(nmds2)
gof <-goodness(nmds2)
gof
plot(nmds2, display = "sites", type="n")
points(nmds2, display="sites",cex=2*gof/mean(gof))
dev.off()

# Use ADONIS to Test for interactions among groups

# Use strata so randomizations occur within each site)
set.seed(1234)
adonis(sor~site*layer*expt, data=ENV2, permutations=999, strata=ENV2$site)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# site              1     0.371 0.37100  1.1129 0.00751  0.020 *  
#   layer             2     3.509 1.75431  5.2626 0.07103  0.001 ***
#   expt              5     2.173 0.43450  1.3034 0.04398  0.056 .  
# site:layer        2     0.701 0.35051  1.0515 0.01419  0.325    
# site:expt         5     1.332 0.26632  0.7989 0.02696  0.908    
# layer:expt       10     3.272 0.32725  0.9817 0.06625  0.512    
# site:layer:expt  10     2.371 0.23712  0.7113 0.04800  0.997    
# Residuals       107    35.669 0.33335         0.72208           
# Total           142    49.397                 1.00000    

###########################################
#### Create Fig 4C plot, just 1C1E & 1C3E samples (1C3E subsampled down to 8 to compare with 1C1E for balanced design)
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

# remove outliers
# NZC85_20160701_1C1E_15O is the only NZC sample to cluster with ILC samples
# ILC45_20160701_1C1E_23M does not cluster with anything else
outlier<-c("NZC85_20160701_1C1E_15O","ILC45_20160701_1C1E_23M")
notnull3<-notnull2[!rownames(notnull2) %in% outlier, ]

# Calculate 15th percentile for rrarefy function
percentile<-quantile(rowSums(notnull3), prob=0.15)

set.seed(1234)

# Rarefy the dataset down to the 15th percentile
df<-rrarefy(notnull3, sample=percentile)

# Convert to presence-absence matrix
df[df>0] <-1

# Scree plots to determine number of dimensions to use for NMDS
pdf("Scree_1E3E_outliersRemoved.pdf")
# check dims
dimcheckMDS(df)
dev.off()

### Do NMDS

#Do 2 dimensional NMDS
nmds2<-metaMDS(df, k=2,trymax=100)
# stress = 0.123
# linear R2 = 0.932, non-linear R2 = 0.985 (from stressplot below)

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

lay<-rbind(c(1,2), 
           c(3,4))

g<-grid.arrange(p1, p2, p3, layout_matrix=lay)

ggsave("F2_NMDS.pdf", g, width = 8, height = 9, units = c("in"))

g2<-grid.arrange(p1b, p2b, p3b, layout_matrix=lay)
ggsave("F2_NMDSb.pdf", g2, width = 8, height = 9, units = c("in"))

# Prep metadata format
# Only keep the samples that are in df
ENV2<-subset(ENV,rownames(ENV) %in% rownames(df) )

# Create factors
ENV2$expt<-factor(ENV2$expt, 
                  levels=c("1C1E","1C3E"), labels=c("1C1E","1C3E"))
ENV2$layer<-factor(ENV2$layer, levels=c("Moss","Organic","Mineral"), labels=c("Bryophyte","Organic","Mineral"))

# Assess dispersion (variance) using ANOVA
# Create distance matrix based on P-A data using binary Bray Curtis (Sorensen) dissimilarity
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
anova(bd_layer) # p-value = 0.02692*
# See Anderson & Walsh, 2013 ref above

pdf("BetaDispersion_boxplot_C_outliersremoved.pdf")
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
pdf("stressplot_C_outliersremoved.pdf")
stressplot(nmds2)
gof <-goodness(nmds2)
gof
plot(nmds2, display = "sites", type="n")
points(nmds2, display="sites",cex=2*gof/mean(gof))
dev.off()

# Use ADONIS to test for interactions among groups

# use strata so randomizations occurr within site
set.seed(1234)
adonis(sor~site*layer*expt, data=ENV2, permutations=999, strata = ENV2$site)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# site             1     0.167 0.16653  0.5036 0.00498  0.001 ***
#   layer            2     0.992 0.49604  1.5000 0.02969  0.047 *  
#   expt             1     3.092 3.09203  9.3499 0.09254  0.001 ***
#   site:layer       2     0.699 0.34949  1.0568 0.02092  0.371    
# site:expt        1     0.153 0.15305  0.4628 0.00458  0.998    
# layer:expt       2     0.774 0.38715  1.1707 0.02317  0.206    
# site:layer:expt  2     1.079 0.53966  1.6319 0.03230  0.030 *  
#   Residuals       80    26.456 0.33070         0.79180           
# Total           91    33.412                 1.00000

# The significant interaction means that tests of main effects (individual factors alone) may not be meaningful