# Teresita M. Porter, July 8, 2019

library(vegan)
library(indicspecies)
library(reshape2)
library(corrplot)
library(metacoder)
library(grid)
library(tibble)
library(dplyr)
library(readr)
library(ggplot2)
library(gridExtra)
library(dichromat)
library("psych") #corr.test # remove
library(Hmisc)

###################################################################
# NEW FUNCTION TO EXTRACT SUBSTRING FROM THE RIGHT
###################################################################
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

###################################################################
# NEW FUNCTION TO MERGE order_family IF fBP <= 0.20
###################################################################
fix_family <- function (df) {
  for (row in 1:nrow(df)) {
    
    order<-df[row,"order"]
    familyBP<-df[row,"fBP"]
    
    if(familyBP<= 0.20) {
      df[row,"family"]<-"F_unk"
      df[row,"fBP"]<-"NA"
    }
    else {
      next
    }
  }
  return(df)
}

###################################################################
# NEW FUNCTION TO MERGE family_genus IF gBP <= 0.30
###################################################################
fix_genus <- function (df) {
  for (row in 1:nrow(df)) {
    
    family<-df[row,"family"]
    esv<-df[row,"marker_OTU"]
    genusBP<-df[row,"gBP"]
    
    if(genusBP<= 0.30) {
      df[row,"genus"]<-"G_unk"
      df[row,"gBP"]<-"NA"
    }
    else {
      next
    }
  }
  return(df)
}

###################################################################
# NEW FUNCTION TO MERGE genus_species IF sBP <= 0.70 (95% correct)
###################################################################
fix_species <- function (df) {
  for (row in 1:nrow(df)) {
    
    genus<-df[row,"genus"]
    esv<-df[row,"marker_OTU"]
    speciesBP<-df[row,"sBP"]
    
    if(speciesBP<= 0.70) {
      df[row,"species"]<-"S_unk"
      df[row,"sBP"]<-"NA"
    }
    else {
      next
    }
  }
  return(df)
}
###################################################################

# Read infile
A<-read.table(file="LV2016_2.csv", head=TRUE, sep=",")

# Select phylum Arthropoda only
Arth_df<-A[A$phylum=="Arthropoda",]

# Select samples from 1C1E (n=8 x 3 layers)
B1<-Arth_df[Arth_df$cores==1 & Arth_df$extractions==1,]

# Randomly subsample 4 unique gridcoord without replacement
B1.gridcoord<-sample(unique(B1$gridcoord), 4, replace=FALSE)

# Retrieve gridcoord
B1<-B1[B1$gridcoord=="11" |
         B1$gridcoord=="31" |
         B1$gridcoord=="51" |
         B1$gridcoord=="23",]

# Select samples from 1C3E (n=36 x 3 layers)
B2<-Arth_df[Arth_df$cores==1 & Arth_df$extractions==3,]

# Grab same gridcoord as 1C1E above for fair comparison
B2<-B2[B2$gridcoord=="11" |
         B2$gridcoord=="31" |
         B2$gridcoord=="51" |
         B2$gridcoord=="23",]

# Grab XC3E (n=4 x 3 layers x 5 pooled treatments)
B3<-Arth_df[Arth_df$cores!=1 & Arth_df$extractions==3,]

# combine B1, B2, and B3
B4<-rbind(B1,B2,B3)

# Pivot to make matrices for vegan
matrix.1<-dcast(B4, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)

# Move marker_OTU to row names
rownames(matrix.1)<-matrix.1$marker_OTU
matrix.1<-matrix.1[,-1]

# Transpose to get sites in rows, ESVs in columns
matrix2.1<-t(matrix.1)

# Remove columns that sum to zero
notnull.1<-matrix2.1[,colSums(matrix2.1) !=0]

#remove rows that sum tozero
notnull2.1<-notnull.1[rowSums(notnull.1) !=0,]

#calculate 15th percentile for rrarefy function
percentile.1<-quantile(rowSums(notnull2.1), prob=0.15)

set.seed(1234)

#Rarefy down to 15th percentile library size to normalize read depth across samples
df.1<-rrarefy(notnull2.1,sample=percentile.1)

# Convert to presence absence matrix
df.1[df.1>0] <-1

###################################################################
##### Do indicator analysis
###################################################################

# Turn data frame into tables for each experiment
# Used rarefied p-a matrices for indicator species analysis

# Break it down by experiment
X1C1E.1<-data.frame(df.1[grepl("_1C1E_", rownames(df.1)),])
X1C3E.1<-data.frame(df.1[grepl("_1C3E_", rownames(df.1)),])

# Break it down by cores
X2C3E.2<-data.frame(df.1[grepl("_2C3E_",rownames(df.1)),])
X4C3E.2<-data.frame(df.1[grepl("_4C3E_",rownames(df.1)),])
X6C3E.2<-data.frame(df.1[grepl("_6C3E_",rownames(df.1)),])
X8C3E.2<-data.frame(df.1[grepl("_8C3E_",rownames(df.1)),])
X915C3E.2<-data.frame(df.1[(grepl("_9C3E_",rownames(df.1)) |
                            grepl("_12C3E_",rownames(df.1)) |
                            grepl("_13C3E_",rownames(df.1)) |
                            grepl("_14C3E_",rownames(df.1)) |
                            grepl("_15C3E_",rownames(df.1)) ),])

# Sort out the groups for each experiment:  
# look at site indicators
groups_1C1E.1<-c(rep(1,12), 
            rep(2,12))
groups_1C3E.1<-c(rep(1,12), 
               rep(2,12))
groups_915C3E.2<-c(rep(1,12),
                  rep(2,12))
groups_2C3E.2<-c(rep(1,12),
               rep(2,12))
groups_4C3E.2<-c(rep(1,12),
               rep(2,12))
groups_6C3E.2<-c(rep(1,12),
               rep(2,12))
groups_8C3E.2<-c(rep(1,12),
               rep(2,12))

# Do indicator analysis
X1C1E_indval.1=multipatt(X1C1E.1, groups_1C1E.1, control=how(nperm=999))
X1C3E_indval.1=multipatt(X1C3E.1, groups_1C3E.1, control=how(nperm=999))

X915C3E_indval.2=multipatt(X915C3E.2, groups_915C3E.2, control=how(nperm=999))
X2C3E_indval.2=multipatt(X2C3E.2, groups_2C3E.2, control=how(nperm=999))
X4C3E_indval.2=multipatt(X4C3E.2, groups_4C3E.2, control=how(nperm=999))
X6C3E_indval.2=multipatt(X6C3E.2, groups_6C3E.2, control=how(nperm=999))
X8C3E_indval.2=multipatt(X8C3E.2, groups_8C3E.2, control=how(nperm=999))

# Grab the data frame containing the p-values from each analysis
X1C1E_indval_df.1<-data.frame(X1C1E_indval.1$sign)
X1C3E_indval_df.1<-data.frame(X1C3E_indval.1$sign)

X915C3E_indval_df.2<-data.frame(X915C3E_indval.2$sign)
X2C3E_indval_df.2<-data.frame(X2C3E_indval.2$sign)
X4C3E_indval_df.2<-data.frame(X4C3E_indval.2$sign)
X6C3E_indval_df.2<-data.frame(X6C3E_indval.2$sign)
X8C3E_indval_df.2<-data.frame(X8C3E_indval.2$sign)

# Grab ILC indicator ESVs with pvalues <= 0.05
ILC_1C1E_indic.1<-X1C1E_indval_df.1[which(X1C1E_indval_df.1$p.value<=0.05 & 
                                              X1C1E_indval_df.1$s.1==1 &
                                              X1C1E_indval_df.1$s.2==0),]
ILC_1C3E_indic.1<-X1C3E_indval_df.1[which(X1C3E_indval_df.1$p.value<=0.05 & 
                                              X1C3E_indval_df.1$s.1==1 &
                                              X1C3E_indval_df.1$s.2==0),]

ILC_915C3E_indic.2<-X915C3E_indval_df.2[which(X915C3E_indval_df.2$p.value<=0.05 & 
                                              X915C3E_indval_df.2$s.1==1 &
                                              X915C3E_indval_df.2$s.2==0),]
ILC_2C3E_indic.2<-X2C3E_indval_df.2[which(X2C3E_indval_df.2$p.value<=0.05 & 
                                                     X2C3E_indval_df.2$s.1==1 &
                                                     X2C3E_indval_df.2$s.2==0),]
ILC_4C3E_indic.2<-X4C3E_indval_df.2[which(X4C3E_indval_df.2$p.value<=0.05 & 
                                                     X4C3E_indval_df.2$s.1==1 &
                                                     X4C3E_indval_df.2$s.2==0),]
ILC_6C3E_indic.2<-X6C3E_indval_df.2[which(X6C3E_indval_df.2$p.value<=0.05 & 
                                                     X6C3E_indval_df.2$s.1==1 &
                                                     X6C3E_indval_df.2$s.2==0),]
ILC_8C3E_indic.2<-X8C3E_indval_df.2[which(X8C3E_indval_df.2$p.value<=0.05 & 
                                                     X8C3E_indval_df.2$s.1==1 &
                                                     X8C3E_indval_df.2$s.2==0),]

# Add site to the data frame
ILC_1C1E_indic.1$site <- "ILC"
ILC_1C3E_indic.1$site <- "ILC"

ILC_915C3E_indic.2$site <- "ILC"
ILC_2C3E_indic.2$site <- "ILC"
ILC_4C3E_indic.2$site <- "ILC"
ILC_6C3E_indic.2$site <- "ILC"
ILC_8C3E_indic.2$site <- "ILC"

# Add experiment to the data frame
ILC_1C1E_indic.1$expt <- "1C1E"
ILC_1C3E_indic.1$expt <- "1C3E"

ILC_915C3E_indic.2$expt <- "915C3E"
ILC_2C3E_indic.2$expt <- "2C3E"
ILC_4C3E_indic.2$expt <- "4C3E"
ILC_6C3E_indic.2$expt <- "6C3E"
ILC_8C3E_indic.2$expt <- "8C3E"

# Move rownames to column then reset row names
ILC_1C1E_indic.1$marker_OTU <- rownames(ILC_1C1E_indic.1)
ILC_1C3E_indic.1$marker_OTU <- rownames(ILC_1C3E_indic.1)

ILC_915C3E_indic.2$marker_OTU <- rownames(ILC_915C3E_indic.2)
ILC_2C3E_indic.2$marker_OTU <- rownames(ILC_2C3E_indic.2)
ILC_4C3E_indic.2$marker_OTU <- rownames(ILC_4C3E_indic.2)
ILC_6C3E_indic.2$marker_OTU <- rownames(ILC_6C3E_indic.2)
ILC_8C3E_indic.2$marker_OTU <- rownames(ILC_8C3E_indic.2)

rownames(ILC_1C1E_indic.1) <- NULL
rownames(ILC_1C3E_indic.1) <- NULL

rownames(ILC_915C3E_indic.2) <- NULL
rownames(ILC_2C3E_indic.2) <- NULL
rownames(ILC_4C3E_indic.2) <- NULL
rownames(ILC_6C3E_indic.2) <- NULL
rownames(ILC_8C3E_indic.2) <- NULL

# Grab NZC indicator ESVs with p-values <= 0.05
NZC_1C1E_indic.1<-X1C1E_indval_df.1[which(X1C1E_indval_df.1$p.value<=0.05 & 
                                        X1C1E_indval_df.1$s.1==0 &
                                        X1C1E_indval_df.1$s.2==1),]
NZC_1C3E_indic.1<-X1C3E_indval_df.1[which(X1C3E_indval_df.1$p.value<=0.05 & 
                                        X1C3E_indval_df.1$s.1==0 &
                                        X1C3E_indval_df.1$s.2==1),]

NZC_915C3E_indic.2<-X915C3E_indval_df.2[which(X915C3E_indval_df.2$p.value<=0.05 & 
                                            X915C3E_indval_df.2$s.1==0 &
                                            X915C3E_indval_df.2$s.2==1),]
NZC_2C3E_indic.2<-X2C3E_indval_df.2[which(X2C3E_indval_df.2$p.value<=0.05 & 
                                        X2C3E_indval_df.2$s.1==0 &
                                        X2C3E_indval_df.2$s.2==1),]
NZC_4C3E_indic.2<-X4C3E_indval_df.2[which(X4C3E_indval_df.2$p.value<=0.05 & 
                                        X4C3E_indval_df.2$s.1==0 &
                                        X4C3E_indval_df.2$s.2==1),]
NZC_6C3E_indic.2<-X6C3E_indval_df.2[which(X6C3E_indval_df.2$p.value<=0.05 & 
                                        X6C3E_indval_df.2$s.1==0 &
                                        X6C3E_indval_df.2$s.2==1),]
NZC_8C3E_indic.2<-X8C3E_indval_df.2[which(X8C3E_indval_df.2$p.value<=0.05 & 
                                        X8C3E_indval_df.2$s.1==0 &
                                        X8C3E_indval_df.2$s.2==1),]

# Add site to the data frame
NZC_1C1E_indic.1$site <- "NZC"
NZC_1C3E_indic.1$site <- "NZC"

NZC_915C3E_indic.2$site <- "NZC"
NZC_2C3E_indic.2$site <- "NZC"
NZC_4C3E_indic.2$site <- "NZC"
NZC_6C3E_indic.2$site <- "NZC"
NZC_8C3E_indic.2$site <- "NZC"

# Add experiment to the data frame
NZC_1C1E_indic.1$expt <- "1C1E"
NZC_1C3E_indic.1$expt <- "1C3E"

NZC_915C3E_indic.2$expt <- "915C3E"
NZC_2C3E_indic.2$expt <- "2C3E"
NZC_4C3E_indic.2$expt <- "4C3E"
NZC_6C3E_indic.2$expt <- "6C3E"
NZC_8C3E_indic.2$expt <- "8C3E"

# Move rownames to column then reset row names
NZC_1C1E_indic.1$marker_OTU <- rownames(NZC_1C1E_indic.1)
NZC_1C3E_indic.1$marker_OTU <- rownames(NZC_1C3E_indic.1)

NZC_915C3E_indic.2$marker_OTU <- rownames(NZC_915C3E_indic.2)
NZC_2C3E_indic.2$marker_OTU <- rownames(NZC_2C3E_indic.2)
NZC_4C3E_indic.2$marker_OTU <- rownames(NZC_4C3E_indic.2)
NZC_6C3E_indic.2$marker_OTU <- rownames(NZC_6C3E_indic.2)
NZC_8C3E_indic.2$marker_OTU <- rownames(NZC_8C3E_indic.2)

rownames(NZC_1C1E_indic.1) <- NULL
rownames(NZC_1C3E_indic.1) <- NULL

rownames(NZC_915C3E_indic.2) <- NULL
rownames(NZC_2C3E_indic.2) <- NULL
rownames(NZC_4C3E_indic.2) <- NULL
rownames(NZC_6C3E_indic.2) <- NULL
rownames(NZC_8C3E_indic.2) <- NULL

# Combine the data frames
indic.1 <- rbind(
  ILC_1C1E_indic.1, 
  ILC_1C3E_indic.1,
  NZC_1C1E_indic.1, 
  NZC_1C3E_indic.1)

indic.2 <- rbind(
  ILC_915C3E_indic.2,            
  ILC_2C3E_indic.2, 
  ILC_4C3E_indic.2, 
  ILC_6C3E_indic.2,
  ILC_8C3E_indic.2, 
  NZC_915C3E_indic.2, 
  NZC_2C3E_indic.2, 
  NZC_4C3E_indic.2, 
  NZC_6C3E_indic.2,
  NZC_8C3E_indic.2)

# combine into one matrix
indic <- rbind(indic.1, indic.2)

# Add taxonomy to indicator ESVs (keep all indic ESVs, just add the taxonomy) left outer join
merged <- merge(x=indic, y=B4, by="marker_OTU", all.x=TRUE)

# create pivot table populated with ESV counts
pivot <- dcast(merged, site.x + layer + marker_OTU ~ expt, value.var="marker_OTU", fun.aggregate = length)

# convert ESV counts to presence-absence
pivot[,4:10][pivot[,4:10] >0] <-1

#######################################################
### Fix lineage format
#######################################################

# select indicator ESVs based on marker_OTU from original balanced matrix
Z<-as.data.frame(as.matrix(B4), stringsAsFactors = FALSE)

# filter by bootstrap support
Zf<- fix_family(Z)
Zg<- fix_genus(Zf)
Zs<- fix_species(Zg)

# merge lineage into single column
Zs$colineage <- paste(Zs[,15], Zs[,14], sep ="__" )
Zs$sklineage <- paste(Zs[,18], Zs[,17], sep = "__" )
Zs$klineage <-  paste(Zs[,21], Zs[,20], sep = "__" )
Zs$plineage <-  paste(Zs[,24], Zs[,23], sep = "__" )
Zs$clineage <-  paste(Zs[,27], Zs[,26], sep = "__" )
Zs$olineage <-  paste(Zs[,30], Zs[,29], sep = "__" )
Zs$flineage <-  paste(Zs[,33], Zs[,32], sep = "__" )
Zs$glineage <-  paste(Zs[,36], Zs[,35], sep = "__" )
Zs$slineage <-  paste(Zs[,39], Zs[,38], sep = "__" )

#start lineage at class rank
Zs$lineage <- apply( Zs[,45:49] , 1 , paste , collapse = ";" )

# add ESV to end of lineage to go beyond species resolution
Zs$lineage <- paste(Zs$lineage, Zs$marker_OTU, sep=";")

# make sure reads are integers
Zs$reads<-as.numeric(factor(Zs$reads))

# make abundance matrix with ESVs in rows and sites in columns
pivot2<-dcast(Zs, marker_OTU + lineage ~ GRDIname, value.var="reads", fun.aggregate = sum)

# Narrow down this abundance matrix to just the indicator ESVs (don't need to be significant)
indic2<-merge(x=indic, y=pivot2, by="marker_OTU", all.x=TRUE)

# Remove unneeded columns
indic3<-indic2[,-c(2:8)]

# remove duplicate rows
indic4<-unique(indic3)

# make all numeric (double) columns integer instead
indic4[,3:length(names(indic4))] <- lapply(indic4[,3:length(names(indic4))], as.integer)

# set any NAs to zero
indic4[is.na(indic4)] <- 0

# reset rownames
rownames(indic4) <- NULL

# rename for readability
esvs<-indic4

## recode rownames in samples 
oldnames<-names(esvs)
newnames<-gsub("15C3E", "915C3E", oldnames)
newnames<-gsub("12C3E", "915C3E", newnames)
newnames<-gsub("13C3E", "915C3E", newnames)
newnames<-gsub("14C3E", "915C3E", newnames)
newnames<-gsub( "9C3E", "915C3E", newnames)
names(esvs)<-newnames

# grab first two columns
esvs1<-esvs[,1:2]

# sort columns alphabetically
esvs2<-esvs[ , sort(names(esvs[3:length(names(esvs))]))]
# ensure ESVs are all integers
esvs2[]<-lapply(esvs2, as.integer)
# Remove columns that sum to zero
esvs2<-esvs2[,colSums(esvs2) !=0]
# remove rows that sum to zero
esvs2<-esvs2[rowSums(esvs2) !=0,]

# merge
esvs<-merge(esvs1,esvs2,by="row.names")
esvs<-esvs[,-1]

############
# Plot just the global site ESVs

# create mapping file
map<-esvs[,1:2]

# melt for ggplot
combo <- melt(pivot, id=c("marker_OTU","site.x","layer"))

# merge map to result combo table
TABLE.1<-merge(combo, map, by=c("marker_OTU", "marker_OTU"))

# remove duplicate rows
TABLE.2<-unique(TABLE.1)

# exclude unknown taxa from lineage
#TABLE.2$lineage<-gsub("(.*);species__S_unk","\\1", TABLE.2$lineage)
#TABLE.2$lineage<-gsub("(.*);genus__G_unk","\\1", TABLE.2$lineage)
#TABLE.2$lineage<-gsub("(.*);family__F_unk","\\1", TABLE.2$lineage)

# exclude ranks from lineage
TABLE.2$lineage<-gsub("(.*)species__(.*)","\\1\\2", TABLE.2$lineage)
TABLE.2$lineage<-gsub("(.*)genus__(.*)","\\1\\2", TABLE.2$lineage)
TABLE.2$lineage<-gsub("(.*)family__(.*)","\\1\\2", TABLE.2$lineage)
TABLE.2$lineage<-gsub("(.*)order__(.*)","\\1\\2", TABLE.2$lineage)
TABLE.2$lineage<-gsub("(.*)class__(.*)","\\1\\2", TABLE.2$lineage)
TABLE.2$lineage<-gsub("phylum__(.*)","\\1", TABLE.2$lineage)

# fix spacing
TABLE.2$lineage<-gsub(";","; ", TABLE.2$lineage)

# fix Plecoptera_Insecta
TABLE.2$lineage<-gsub("(Plecoptera)_Insecta","\\1", TABLE.2$lineage)

# fix Genus_species
TABLE.2$lineage<-gsub(";(.*)_(.*)", ";\\1 \\2", TABLE.2$lineage)

# sort by ascending lineage
TABLE.3<-TABLE.2[order(TABLE.2$lineage),]

# create factors
TABLE.3$site<-factor(TABLE.3$site, levels=c("ILC", "NZC"),
                     labels=c("Island Lake", "Nimitz"))
TABLE.3$lineage <- as.factor(TABLE.3$lineage)

# create a couple presence-absence heat maps
h<-ggplot(TABLE.3, aes(variable,lineage)) +
  geom_tile(aes(fill=site, alpha=value)) +
  facet_wrap(~site) +
  xlab("Method") +
  ylab("ESVs") +
  scale_y_discrete(limits = rev(levels(TABLE.3$lineage))) +
  theme_bw() +
  scale_alpha_continuous(range = c(0, 1)) +
  theme(
        axis.text.x = element_text(size=6),
        axis.title.y = element_text(size=6),
        axis.title.x = element_text(size=6),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 6),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=5)) +
  guides(fill = guide_legend(nrow=1))

pdf("FigS7_indic_heatmap_ESVs.pdf", width = 8, height=10)
print(h)
dev.off()
