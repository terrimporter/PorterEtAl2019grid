# Teresita M. Porter, October 8, 2019

library(vegan)
library(reshape2)
library(gridExtra)
library(grid)
library(ggplot2)
library(plyr)
library(data.table) # rowname to first col
library(pairwiseAdonis)
library(goeveg) # scree
library(ggpubr)

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

# Read infile prepared by python script
## LV2016_2.csv is the clean taxonomic assignment table
A<-read.table(file="LV2016_2.csv", head=TRUE, sep=",")

# Select phylum Arthropoda only
Arth_df<-A[A$phylum=="Arthropoda",]

# Create matrix of just 1C3E samples
B<-Arth_df[(Arth_df$cores==1 & Arth_df$extractions==3),]

# Subsample down to various sample sizes to check for number needed to detect significant differences among sites
# Randomly subsample down to 2 (how low can we go?)
B.2.gridcoord<-sample(unique(B$gridcoord), 2, replace=FALSE)
# Retrieve gridcoord
B.2<-B[B$gridcoord=="21" |
         B$gridcoord=="66",]

# Randomly subsample down to 3 (how low can we go?) (contains same 2 grid coord as above)
B.3.gridcoord<-sample(unique(B$gridcoord), 3, replace=FALSE)
# Retrieve gridcoord
B.3<-B[B$gridcoord=="21" |
         B$gridcoord=="36" |
         B$gridcoord=="66",]

# Randomly subsample down to 4 (same sample size as XC3E expt) (contains same 3 grid coord as above)
B.4.gridcoord<-sample(unique(B$gridcoord), 4, replace=FALSE)
# Retrieve gridcoord
B.4<-B[B$gridcoord=="21" |
         B$gridcoord=="36" |
         B$gridcoord=="64" |
         B$gridcoord=="66",]

# Randomly subsample down to 8 (same sample size as 1E3E expt) (contains same 4 gridcoord as above, too)
B.8.gridcoord<-sample(unique(B$gridcoord), 8, replace=FALSE)
# Retrieve gridcoord
B.8<-B[B$gridcoord=="14" |
         B$gridcoord=="21" |
         B$gridcoord=="25" |
         B$gridcoord=="34" |
         B$gridcoord=="36" |
         B$gridcoord=="53" |
         B$gridcoord=="64" |
         B$gridcoord=="66",]

# Randomly subsample down to 20 (contains same 8 gridcoord as above, too)
set.seed(20)
B.20.gridcoord<-sample(unique(B$gridcoord), 20, replace=FALSE)
# Retrieve gridcoord, keep same 8 as above tough
B.20<-B[B$gridcoord=="14" |
          B$gridcoord=="15" |
          B$gridcoord=="16" |
         B$gridcoord=="21" |
          B$gridcoord=="22" |
         B$gridcoord=="25" |
          B$gridcoord=="32" |
          B$gridcoord=="33" |
         B$gridcoord=="34" |
        B$gridcoord=="35" |
         B$gridcoord=="36" |
          B$gridcoord=="42" |
          B$gridcoord=="46" |
          B$gridcoord=="51" |
         B$gridcoord=="53" |
          B$gridcoord=="54" |
          B$gridcoord=="63" |
         B$gridcoord=="64" |
          B$gridcoord=="65" |
         B$gridcoord=="66",]

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

# Do 2 dimensional NMDS, no environment file for now 
nmds2<-metaMDS(df, k=2,trymax=100)
# stress = 0.1291214
# linear fit R2 = 0.922, non-metric fit R2 = 0.983

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
mergedB$layer<-factor(mergedB$layer, levels=c("B","O","M"), labels=c("Bryophyte", "Organic", "Mineral"))
mergedB$expt<-factor(mergedB$expt, levels=c("1C3E"), labels=c("1"))
mergedB$site<-factor(mergedB$site, levels=c("NZC85", "ILC45"), labels=c("Nimitz", "Island Lake"))

# Compile coord for convex hulls
chulls <- ddply(mergedB, .(site), function(mergedB) mergedB[chull(mergedB$NMDS1, mergedB$NMDS2), ])

p1<-ggplot(data=mergedB, aes(x=NMDS1, y=NMDS2)) + 
  ggtitle("A) 1C3E") +
  geom_polygon(data=chulls, aes(x=NMDS1, y=NMDS2, fill=site), alpha=0.25) +
  geom_point(aes(color=layer), position="jitter") +
  scale_fill_manual(name="Site", values=c("black","grey")) +
  scale_color_manual(name="Strata", values=c("#4DAF4A","#377EB8","#FF7F00")) +
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
         color = guide_legend(order=2))

# get legend to plot in separate frame
l <- get_legend(p1)

p1b <- ggplot(data=mergedB, aes(x=NMDS1, y=NMDS2)) + 
  ggtitle("A) 1C3E") +
  geom_polygon(data=chulls, aes(x=NMDS1, y=NMDS2, fill=site), alpha=0.25) +
  geom_point(aes(color=layer), position="jitter") +
 # geom_text(aes(color=layer, label=expt), size = 3, position="jitter") +
  scale_fill_manual(name="Site", values=c("black","grey")) +
  scale_color_manual(name="Strata", values=c("#4DAF4A","#377EB8","#FF7F00")) +
  theme_bw() +
  theme(
    text = element_text(size=12),
    axis.text.x = element_text(hjust=1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(), 
    legend.position = "none")

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

# Shephards curve and goodness of fit calcs
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
set.seed(1234)
adonis(sor ~ site * layer, data=ENV2, permutations=999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# site         1     6.539  6.5395 18.9193 0.08243  0.001 ***
#   layer        2     0.975  0.4876  1.4108 0.01229  0.049 *  
#   site:layer   2     0.958  0.4792  1.3864 0.01208  0.055 .  
# Residuals  205    70.858  0.3456         0.89319           
# Total      210    79.331                 1.00000  

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
#percentile<-quantile(rowSums(notnull2), prob=0.15)

# Set random seed
set.seed(1234)

# Rarefy the dataset down to the 15th percentile
df<-rrarefy(notnull3, sample=percentile)
#df<-rrarefy(notnull2, sample=percentile)

# Convert to presence-absence matrix
df[df>0] <-1

# # # Scree plots to determine number of dimensions to use for NMDS
# pdf("Scree_XC3E_THISONE_outlierremoved.pdf")
# # check dims
#  dimcheckMDS(df)
#  dev.off()

### Do NMDS
#Do 2 dimensional NMDS, no environment file for now 
nmds2<-metaMDS(df, k=2,trymax=100)
# stress = 0.106
# linear R2 = 0.989, non-metric fit = 0.95
# NZC85_20160701_6C3E_S3M   1.74723771 -1.62779967 outlier?

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
                                 labels=c("1","2","4","6","8","15"))
mergedC$site<-factor(mergedC$site, levels=c("NZC85", "ILC45"), 
                                 labels=c("Nimitz", "Island Lake"))

# Create convex hulls
chulls <- ddply(mergedC, .(site), function(mergedC) mergedC[chull(mergedC$NMDS1, mergedC$NMDS2), ])

p2<-ggplot(data=mergedC, aes(x=NMDS1, y=NMDS2)) + 
  ggtitle("B) XC3E") +
  geom_polygon(data=chulls, aes(x=NMDS1, y=NMDS2, fill=site), alpha=0.25) +
  geom_text(aes(color=layer, label=expt), size=3, position="jitter") +
  scale_fill_manual(name="Site", values=c("black","grey")) +
  scale_color_manual(name="Strata", values=c("#4DAF4A","#377EB8","#FF7F00")) +
  theme_bw() +
  theme(
    text = element_text(size=12),
    axis.text.x = element_text(hjust=1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(),
    legend.position="none")

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
                 labels=c("1","2","4","6","8","15"))
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

#PERMANOVA (which is basically adonis()) was found to be largely unaffected by heterogeneity in Anderson & Walsh's simulations but only for balanced designs.
#Anderson MJ, Walsh DCI. PERMANOVA, ANOSIM, and the Mantel test in the face of heterogeneous dispersions: What null hypothesis are you testing? Ecological monographs [Internet] 2013; 83: 557. Available from: http://doi.org/10.1890/12-2010.1

pdf("BetaDispersion_boxplot_B_THISONE_outlierremoved.pdf")
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
pdf("stressplot_B_THISONE_outlierremoved.pdf")
stressplot(nmds2)
gof <-goodness(nmds2)
gof
plot(nmds2, display = "sites", type="n")
points(nmds2, display="sites",cex=2*gof/mean(gof))
dev.off()

# Use ADONIS to test significance of groupings

# 1. Test the interaction first to see how to proceed
set.seed(1234)
adonis(sor~site*layer*expt, data=ENV2, permutations=999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# site              1     0.371 0.37100  1.1129 0.00751  0.279    
# layer             2     3.509 1.75431  5.2626 0.07103  0.001 ***
#   expt              5     2.173 0.43450  1.3034 0.04398  0.050 *  
#   site:layer        2     0.701 0.35051  1.0515 0.01419  0.360    
# site:expt         5     1.332 0.26632  0.7989 0.02696  0.886    
# layer:expt       10     3.272 0.32725  0.9817 0.06625  0.527    
# site:layer:expt  10     2.371 0.23712  0.7113 0.04800  0.999    
# Residuals       107    35.669 0.33335         0.72208           
# Total           142    49.397                 1.00000  

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

# remove potential outlier, 
# NZC85_20160701_1C1E_15O is the only NZC sample to cluster with ILC samples
# ILC45_20160701_1C1E_23M does not cluster with anything else
outlier<-c("NZC85_20160701_1C1E_15O","ILC45_20160701_1C1E_23M")
notnull3<-notnull2[!rownames(notnull2) %in% outlier, ]

# Calculate 15th percentile for rrarefy function
#percentile<-quantile(rowSums(notnull2), prob=0.15)
percentile<-quantile(rowSums(notnull3), prob=0.15)

set.seed(1234)

# Rarefy the dataset down to the 15th percentile
#df<-rrarefy(notnull2, sample=percentile)
df<-rrarefy(notnull3, sample=percentile)

# Convert to presence-absence matrix
df[df>0] <-1

# # Scree plots to determine number of dimensions to use for NMDS
# pdf("Scree_1E3E.pdf")
# # check dims
# dimcheckMDS(df)
# dev.off()

### Do NMDS

#Do 2 dimensional NMDS, no environment file for now 
nmds2<-metaMDS(df, k=2,trymax=100)
# stress = 0.123
# linear R2 = 0.932, non-linear R2 = 0.985

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
mergedD$site<-factor(mergedD$site, levels=c("NZC85", "ILC45"),
                    labels=c("Nimitz", "Island Lake"))

# Create convex hulls
chulls <- ddply(mergedD, .(site), function(mergedD) mergedD[chull(mergedD$NMDS1, mergedD$NMDS2), ])

p3<-ggplot(data=mergedD, aes(x=NMDS1, y=NMDS2)) + 
  ggtitle("C) 1C1E + 1C3E") +
  geom_polygon(data=chulls, aes(x=NMDS1, y=NMDS2, fill=site), alpha=0.25) +
  geom_text(aes(label=expt, color=layer), size = 3, position="jitter") +
  scale_fill_manual(name="Site", values=c("black","grey")) +
  scale_color_manual(name="Strata", values=c("#4DAF4A","#377EB8","#FF7F00")) +
  theme_bw() +
  theme(
    text = element_text(size=12),
    axis.text.x = element_text(hjust=1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(), 
    legend.position="none") 

lay<-rbind(c(1,2), 
           c(3,4))

g<-grid.arrange(p1b, p2, p3, l, layout_matrix=lay)
ggsave("F3_NMDS.pdf", g, width = 8, height = 9, units = c("in"))

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
anova(bd_layer) # p-value = 0.02692*

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
set.seed(1234)
adonis(sor~site*layer*expt, data=ENV2, permutations=999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# site             1     0.167 0.16653  0.5036 0.00498  0.995    
# layer            2     0.992 0.49604  1.5000 0.02969  0.054 .  
# expt             1     3.092 3.09203  9.3499 0.09254  0.001 ***
#   site:layer       2     0.699 0.34949  1.0568 0.02092  0.344    
# site:expt        1     0.153 0.15305  0.4628 0.00458  0.997    
# layer:expt       2     0.774 0.38715  1.1707 0.02317  0.214    
# site:layer:expt  2     1.079 0.53966  1.6319 0.03230  0.025 *  
#   Residuals       80    26.456 0.33070         0.79180           
# Total           91    33.412                 1.00000   

# The significant interaction means that tests of main effects (individual factors alone) may not be meaningful

# 2. The next step is to do pair-wise comparisons for the factor of interest separately within each level of the other factor
# i.e. pair-wise comparisons of levels of factor A (sites) within each level of factor B (layers)

# cite pairwiseAdonis
# Martinez Arbizu, P. (2017). pairwiseAdonis: Pairwise multilevel comparison using adonis. R package version 0.0.1.

pairwise.adonis(sor, paste(ENV2$site, ENV2$layer, ENV2$expt))
# sites, layers, expt doesn't explain any sig differences
# pairs   F.Model         R2 p.value p.adjusted sig
# 1    NZC85 Bryophyte 1C3E vs NZC85 Organic 1C3E 0.8573559 0.06187009   0.434      1.000    
# 2    NZC85 Bryophyte 1C3E vs NZC85 Mineral 1C3E 0.8339700 0.05622028   0.631      1.000    
# 3  NZC85 Bryophyte 1C3E vs ILC45 Bryophyte 1C3E 0.6015637 0.04119858   0.837      1.000    
# 4    NZC85 Bryophyte 1C3E vs ILC45 Organic 1C3E 0.8501052 0.05724574   0.468      1.000    
# 5    NZC85 Bryophyte 1C3E vs ILC45 Mineral 1C3E 0.9884805 0.06594935   0.378      1.000    
# 6  NZC85 Bryophyte 1C3E vs NZC85 Bryophyte 1C1E 3.0496303 0.17886782   0.009      0.594    
# 7    NZC85 Bryophyte 1C3E vs NZC85 Organic 1C1E 2.4679666 0.15955339   0.004      0.264    
# 8    NZC85 Bryophyte 1C3E vs NZC85 Mineral 1C1E 4.0120672 0.22274329   0.001      0.066    
# 9  NZC85 Bryophyte 1C3E vs ILC45 Bryophyte 1C1E 2.2895426 0.14055291   0.014      0.924    
# 10   NZC85 Bryophyte 1C3E vs ILC45 Organic 1C1E 1.9575001 0.13087081   0.040      1.000    
# 11   NZC85 Bryophyte 1C3E vs ILC45 Mineral 1C1E 2.3139744 0.15110214   0.023      1.000    
# 12     NZC85 Organic 1C3E vs NZC85 Mineral 1C3E 0.8442254 0.06098032   0.490      1.000    
# 13   NZC85 Organic 1C3E vs ILC45 Bryophyte 1C3E 0.6519581 0.04775565   0.747      1.000    
# 14     NZC85 Organic 1C3E vs ILC45 Organic 1C3E 0.3491709 0.02615675   0.992      1.000    
# 15     NZC85 Organic 1C3E vs ILC45 Mineral 1C3E 0.9998958 0.07142166   0.393      1.000    
# 16   NZC85 Organic 1C3E vs NZC85 Bryophyte 1C1E 2.0534676 0.13641160   0.024      1.000    
# 17     NZC85 Organic 1C3E vs NZC85 Organic 1C1E 2.3251951 0.16231508   0.011      0.726    
# 18     NZC85 Organic 1C3E vs NZC85 Mineral 1C1E 4.4813623 0.25635086   0.001      0.066    
# 19   NZC85 Organic 1C3E vs ILC45 Bryophyte 1C1E 2.4100742 0.15639601   0.021      1.000    
# 20     NZC85 Organic 1C3E vs ILC45 Organic 1C1E 1.6474274 0.12071341   0.067      1.000    
# 21     NZC85 Organic 1C3E vs ILC45 Mineral 1C1E 2.1908140 0.15438255   0.020      1.000    
# 22   NZC85 Mineral 1C3E vs ILC45 Bryophyte 1C3E 0.9136122 0.06126029   0.472      1.000    
# 23     NZC85 Mineral 1C3E vs ILC45 Organic 1C3E 0.7197910 0.04889954   0.747      1.000    
# 24     NZC85 Mineral 1C3E vs ILC45 Mineral 1C3E 0.4402262 0.03048610   0.991      1.000    
# 25   NZC85 Mineral 1C3E vs NZC85 Bryophyte 1C1E 2.6790916 0.16062575   0.003      0.198    
# 26     NZC85 Mineral 1C3E vs NZC85 Organic 1C1E 2.3479058 0.15297891   0.002      0.132    
# 27     NZC85 Mineral 1C3E vs NZC85 Mineral 1C1E 4.3154479 0.23561793   0.001      0.066    
# 28   NZC85 Mineral 1C3E vs ILC45 Bryophyte 1C1E 2.4415619 0.14849939   0.006      0.396    
# 29     NZC85 Mineral 1C3E vs ILC45 Organic 1C1E 1.8307402 0.12344227   0.027      1.000    
# 30     NZC85 Mineral 1C3E vs ILC45 Mineral 1C1E 2.0381063 0.13552945   0.010      0.660    
# 31   ILC45 Bryophyte 1C3E vs ILC45 Organic 1C3E 0.7373031 0.05002972   0.692      1.000    
# 32   ILC45 Bryophyte 1C3E vs ILC45 Mineral 1C3E 0.9636132 0.06439709   0.448      1.000    
# 33 ILC45 Bryophyte 1C3E vs NZC85 Bryophyte 1C1E 2.5522658 0.15419435   0.010      0.660    
# 34   ILC45 Bryophyte 1C3E vs NZC85 Organic 1C1E 2.4104933 0.15641896   0.001      0.066    
# 35   ILC45 Bryophyte 1C3E vs NZC85 Mineral 1C1E 4.6260779 0.24836565   0.001      0.066    
# 36 ILC45 Bryophyte 1C3E vs ILC45 Bryophyte 1C1E 2.4637759 0.14964829   0.015      0.990    
# 37   ILC45 Bryophyte 1C3E vs ILC45 Organic 1C1E 1.8179415 0.12268516   0.031      1.000    
# 38   ILC45 Bryophyte 1C3E vs ILC45 Mineral 1C1E 2.2508744 0.14758986   0.010      0.660    
# 39     ILC45 Organic 1C3E vs ILC45 Mineral 1C3E 0.8661566 0.05826365   0.510      1.000    
# 40   ILC45 Organic 1C3E vs NZC85 Bryophyte 1C1E 2.2583074 0.13890176   0.010      0.660    
# 41     ILC45 Organic 1C3E vs NZC85 Organic 1C1E 2.2653976 0.14840083   0.015      0.990    
# 42     ILC45 Organic 1C3E vs NZC85 Mineral 1C1E 4.3602237 0.23748206   0.001      0.066    
# 43   ILC45 Organic 1C3E vs ILC45 Bryophyte 1C1E 2.3532603 0.14390160   0.022      1.000    
# 44     ILC45 Organic 1C3E vs ILC45 Organic 1C1E 1.6244340 0.11107671   0.059      1.000    
# 45     ILC45 Organic 1C3E vs ILC45 Mineral 1C1E 2.0680402 0.13724679   0.032      1.000    
# 46   ILC45 Mineral 1C3E vs NZC85 Bryophyte 1C1E 3.0070159 0.17681032   0.002      0.132    
# 47     ILC45 Mineral 1C3E vs NZC85 Organic 1C1E 2.0922566 0.13863113   0.009      0.594    
# 48     ILC45 Mineral 1C3E vs NZC85 Mineral 1C1E 3.8660109 0.21638915   0.002      0.132    
# 49   ILC45 Mineral 1C3E vs ILC45 Bryophyte 1C1E 2.2595561 0.13896788   0.012      0.792    
# 50     ILC45 Mineral 1C3E vs ILC45 Organic 1C1E 1.9176437 0.12854870   0.019      1.000    
# 51     ILC45 Mineral 1C3E vs ILC45 Mineral 1C1E 2.0322022 0.13518992   0.023      1.000    
# 52   NZC85 Bryophyte 1C1E vs NZC85 Organic 1C1E 2.2315289 0.14650722   0.012      0.792    
# 53   NZC85 Bryophyte 1C1E vs NZC85 Mineral 1C1E 5.2528687 0.27283564   0.001      0.066    
# 54 NZC85 Bryophyte 1C1E vs ILC45 Bryophyte 1C1E 2.0791411 0.12930673   0.034      1.000    
# 55   NZC85 Bryophyte 1C1E vs ILC45 Organic 1C1E 0.7283850 0.05305686   0.833      1.000    
# 56   NZC85 Bryophyte 1C1E vs ILC45 Mineral 1C1E 1.6547453 0.11291532   0.065      1.000    
# 57     NZC85 Organic 1C1E vs NZC85 Mineral 1C1E 2.0661846 0.13714053   0.001      0.066    
# 58   NZC85 Organic 1C1E vs ILC45 Bryophyte 1C1E 0.9686616 0.06934534   0.414      1.000    
# 59     NZC85 Organic 1C1E vs ILC45 Organic 1C1E 1.1669586 0.08862780   0.253      1.000    
# 60     NZC85 Organic 1C1E vs ILC45 Mineral 1C1E 0.6658312 0.05256909   0.823      1.000    
# 61   NZC85 Mineral 1C1E vs ILC45 Bryophyte 1C1E 1.0411350 0.06921918   0.457      1.000    
# 62     NZC85 Mineral 1C1E vs ILC45 Organic 1C1E 2.5403684 0.16346900   0.006      0.396    
# 63     NZC85 Mineral 1C1E vs ILC45 Mineral 1C1E 1.9064913 0.12789672   0.014      0.924    
# 64   ILC45 Bryophyte 1C1E vs ILC45 Organic 1C1E 0.8939279 0.06433947   0.492      1.000    
# 65   ILC45 Bryophyte 1C1E vs ILC45 Mineral 1C1E 0.7324459 0.05333688   0.753      1.000    
# 66     ILC45 Organic 1C1E vs ILC45 Mineral 1C1E 0.8447854 0.06576874   0.593      1.000 

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
