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

# Subsample 1C1E down to 4, and 1C3E down to 4 to create balanced dataset
B1<-Arth_df[Arth_df$cores==1 & Arth_df$extractions==1,]
# Randomly subsample 8 unique gridcoord without replacement so 1C3E can be subsampled to create balanced design
B1.gridcoord<-sample(unique(B1$gridcoord), 4, replace=FALSE)
# Retrieve gridcoord
B1<-B1[B1$gridcoord=="11" |
         B1$gridcoord=="23" |
         B1$gridcoord=="31" |
         B1$gridcoord=="35",]

B2<-Arth_df[Arth_df$cores==1 & Arth_df$extractions==3,]
# Grab same gridcoord as above for fair comparison
B2<-B2[B2$gridcoord=="11" |
         B2$gridcoord=="23" |
         B2$gridcoord=="31" |
         B2$gridcoord=="35",]
# Grab XC3E
B3<-Arth_df[Arth_df$cores!=1 & Arth_df$extractions==3,]
# combine B1, B2, and B3
B4<-rbind(B1,B2,B3)

# Pivot to make matrix for vegan
matrix<-dcast(B4, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)

# Move marker_OTU to row names
rownames(matrix)<-matrix$marker_OTU
matrix<-matrix[,-1]

# Transpose to get sites in rows, ESVs in columns
matrix2<-t(matrix)

# Remove columns that sum to zero
notnull<-matrix2[,colSums(matrix2) !=0]

#remove rows that sum tozero
notnull2<-notnull[rowSums(notnull) !=0,]

#calculate 15th percentile for rrarefy function
percentile<-quantile(rowSums(notnull2), prob=0.15)

set.seed(1234)

#Rarefy down to 15th percentile library size to normalize read depth across samples
df<-rrarefy(notnull2,sample=percentile)

# Convert to presence absence matrix
df[df>0] <-1

###################################################################
##### Do indicator species analysis
###################################################################

# Turn data frame into tables for each experiment
# Used rarefied p-a matrices for indicator species analysis

# Break it down by experiment
X1C1E<-data.frame(df[grepl("_1C1E_", rownames(df)),])
X1C3E<-data.frame(df[grepl("_1C3E_", rownames(df)),])

# Break out the pooled soil cores for a more detailed correlation analysis
X2C3E<-data.frame(df[grepl("_2C3E_",rownames(df)),])
X4C3E<-data.frame(df[grepl("_4C3E_",rownames(df)),])
X6C3E<-data.frame(df[grepl("_6C3E_",rownames(df)),])
X8C3E<-data.frame(df[grepl("_8C3E_",rownames(df)),])
X915C3E<-data.frame(df[(grepl("_9C3E_",rownames(df)) |
                          grepl("_12C3E_",rownames(df)) |
                          grepl("_13C3E_",rownames(df)) |
                          grepl("_14C3E_",rownames(df)) |
                          grepl("_15C3E_",rownames(df)) ),])

# Sort out the groups for each experiment:  
# look at site indicators
groups_1C1E<-c(rep(1,12), 
               rep(2,12))
groups_1C3E<-c(rep(1,12), 
               rep(2,11))
groups_915C3E<-c(rep(1,12),
                 rep(2,12))
groups_2C3E<-c(rep(1,12),
               rep(2,12))
groups_4C3E<-c(rep(1,12),
               rep(2,12))
groups_6C3E<-c(rep(1,12),
               rep(2,12))
groups_8C3E<-c(rep(1,12),
               rep(2,12))

# Do indicator analysis
X1C1E_indval=multipatt(X1C1E, groups_1C1E, control=how(nperm=999))
X1C3E_indval=multipatt(X1C3E, groups_1C3E, control=how(nperm=999))
X915C3E_indval=multipatt(X915C3E, groups_915C3E, control=how(nperm=999))
X2C3E_indval=multipatt(X2C3E, groups_2C3E, control=how(nperm=999))
X4C3E_indval=multipatt(X4C3E, groups_4C3E, control=how(nperm=999))
X6C3E_indval=multipatt(X6C3E, groups_6C3E, control=how(nperm=999))
X8C3E_indval=multipatt(X8C3E, groups_8C3E, control=how(nperm=999))

# Grab the data frame containing the p-values from each analysis
X1C1E_indval_df<-data.frame(X1C1E_indval$sign)
X1C3E_indval_df<-data.frame(X1C3E_indval$sign)
X915C3E_indval_df<-data.frame(X915C3E_indval$sign)
X2C3E_indval_df<-data.frame(X2C3E_indval$sign)
X4C3E_indval_df<-data.frame(X4C3E_indval$sign)
X6C3E_indval_df<-data.frame(X6C3E_indval$sign)
X8C3E_indval_df<-data.frame(X8C3E_indval$sign)

# Grab ILC indicator ESVs with pvalues <= 0.05
ILC_1C1E_indic<-X1C1E_indval_df[which(X1C1E_indval_df$p.value<=0.05 & 
                                        X1C1E_indval_df$s.1==1 &
                                        X1C1E_indval_df$s.2==0),]
ILC_1C3E_indic<-X1C3E_indval_df[which(X1C3E_indval_df$p.value<=0.05 & 
                                        X1C3E_indval_df$s.1==1 &
                                        X1C3E_indval_df$s.2==0),]
ILC_915C3E_indic<-X915C3E_indval_df[which(X915C3E_indval_df$p.value<=0.05 & 
                                            X915C3E_indval_df$s.1==1 &
                                            X915C3E_indval_df$s.2==0),]
ILC_2C3E_indic<-X2C3E_indval_df[which(X2C3E_indval_df$p.value<=0.05 & 
                                        X2C3E_indval_df$s.1==1 &
                                        X2C3E_indval_df$s.2==0),]
ILC_4C3E_indic<-X4C3E_indval_df[which(X4C3E_indval_df$p.value<=0.05 & 
                                        X4C3E_indval_df$s.1==1 &
                                        X4C3E_indval_df$s.2==0),]
ILC_6C3E_indic<-X6C3E_indval_df[which(X6C3E_indval_df$p.value<=0.05 & 
                                        X6C3E_indval_df$s.1==1 &
                                        X6C3E_indval_df$s.2==0),]
ILC_8C3E_indic<-X8C3E_indval_df[which(X8C3E_indval_df$p.value<=0.05 & 
                                        X8C3E_indval_df$s.1==1 &
                                        X8C3E_indval_df$s.2==0),]

# Add site to the data frame
ILC_1C1E_indic$site <- "ILC"
ILC_1C3E_indic$site <- "ILC"
ILC_915C3E_indic$site <- "ILC"
ILC_2C3E_indic$site <- "ILC"
ILC_4C3E_indic$site <- "ILC"
ILC_6C3E_indic$site <- "ILC"
ILC_8C3E_indic$site <- "ILC"

# Add experiment to the data frame
ILC_1C1E_indic$expt <- "1C1E"
ILC_1C3E_indic$expt <- "1C3E"
ILC_915C3E_indic$expt <- "915C3E"
ILC_2C3E_indic$expt <- "2C3E"
ILC_4C3E_indic$expt <- "4C3E"
ILC_6C3E_indic$expt <- "6C3E"
ILC_8C3E_indic$expt <- "8C3E"

# Move rownames to column then reset row names
ILC_1C1E_indic$marker_OTU <- rownames(ILC_1C1E_indic)
ILC_1C3E_indic$marker_OTU <- rownames(ILC_1C3E_indic)
ILC_915C3E_indic$marker_OTU <- rownames(ILC_915C3E_indic)
ILC_2C3E_indic$marker_OTU <- rownames(ILC_2C3E_indic)
ILC_4C3E_indic$marker_OTU <- rownames(ILC_4C3E_indic)
ILC_6C3E_indic$marker_OTU <- rownames(ILC_6C3E_indic)
ILC_8C3E_indic$marker_OTU <- rownames(ILC_8C3E_indic)

rownames(ILC_1C1E_indic) <- NULL
rownames(ILC_1C3E_indic) <- NULL
rownames(ILC_915C3E_indic) <- NULL
rownames(ILC_2C3E_indic) <- NULL
rownames(ILC_4C3E_indic) <- NULL
rownames(ILC_6C3E_indic) <- NULL
rownames(ILC_8C3E_indic) <- NULL

# Grab NZC indicator ESVs with pvalues <= 0.05
NZC_1C1E_indic<-X1C1E_indval_df[which(X1C1E_indval_df$p.value<=0.05 & 
                                        X1C1E_indval_df$s.1==0 &
                                        X1C1E_indval_df$s.2==1),]
NZC_1C3E_indic<-X1C3E_indval_df[which(X1C3E_indval_df$p.value<=0.05 & 
                                        X1C3E_indval_df$s.1==0 &
                                        X1C3E_indval_df$s.2==1),]
NZC_915C3E_indic<-X915C3E_indval_df[which(X915C3E_indval_df$p.value<=0.05 & 
                                            X915C3E_indval_df$s.1==0 &
                                            X915C3E_indval_df$s.2==1),]
NZC_2C3E_indic<-X2C3E_indval_df[which(X2C3E_indval_df$p.value<=0.05 & 
                                        X2C3E_indval_df$s.1==0 &
                                        X2C3E_indval_df$s.2==1),]
NZC_4C3E_indic<-X4C3E_indval_df[which(X4C3E_indval_df$p.value<=0.05 & 
                                        X4C3E_indval_df$s.1==0 &
                                        X4C3E_indval_df$s.2==1),]
NZC_6C3E_indic<-X6C3E_indval_df[which(X6C3E_indval_df$p.value<=0.05 & 
                                        X6C3E_indval_df$s.1==0 &
                                        X6C3E_indval_df$s.2==1),]
NZC_8C3E_indic<-X8C3E_indval_df[which(X8C3E_indval_df$p.value<=0.05 & 
                                        X8C3E_indval_df$s.1==0 &
                                        X8C3E_indval_df$s.2==1),]

# Add site to the data frame
NZC_1C1E_indic$site <- "NZC"
NZC_1C3E_indic$site <- "NZC"
NZC_915C3E_indic$site <- "NZC"
NZC_2C3E_indic$site <- "NZC"
NZC_4C3E_indic$site <- "NZC"
NZC_6C3E_indic$site <- "NZC"
NZC_8C3E_indic$site <- "NZC"

# Add experiment to the data frame
NZC_1C1E_indic$expt <- "1C1E"
NZC_1C3E_indic$expt <- "1C3E"
NZC_915C3E_indic$expt <- "915C3E"
NZC_2C3E_indic$expt <- "2C3E"
NZC_4C3E_indic$expt <- "4C3E"
NZC_6C3E_indic$expt <- "6C3E"
NZC_8C3E_indic$expt <- "8C3E"

# Move rownames to column then reset row names
NZC_1C1E_indic$marker_OTU <- rownames(NZC_1C1E_indic)
NZC_1C3E_indic$marker_OTU <- rownames(NZC_1C3E_indic)
NZC_915C3E_indic$marker_OTU <- rownames(NZC_915C3E_indic)
NZC_2C3E_indic$marker_OTU <- rownames(NZC_2C3E_indic)
NZC_4C3E_indic$marker_OTU <- rownames(NZC_4C3E_indic)
NZC_6C3E_indic$marker_OTU <- rownames(NZC_6C3E_indic)
NZC_8C3E_indic$marker_OTU <- rownames(NZC_8C3E_indic)

rownames(NZC_1C1E_indic) <- NULL
rownames(NZC_1C3E_indic) <- NULL
rownames(NZC_915C3E_indic) <- NULL
rownames(NZC_2C3E_indic) <- NULL
rownames(NZC_4C3E_indic) <- NULL
rownames(NZC_6C3E_indic) <- NULL
rownames(NZC_8C3E_indic) <- NULL

# Combine the data frames
indic <- rbind(
  ILC_1C1E_indic, 
  ILC_1C3E_indic, 
  ILC_915C3E_indic,            
  ILC_2C3E_indic, 
  ILC_4C3E_indic, 
  ILC_6C3E_indic,
  ILC_8C3E_indic, 
  NZC_1C1E_indic, 
  NZC_1C3E_indic, 
  NZC_915C3E_indic,            
  NZC_2C3E_indic, 
  NZC_4C3E_indic, 
  NZC_6C3E_indic,
  NZC_8C3E_indic)

#######################################################
### Create heat trees for indicator taxa
#######################################################

# select indicator taxa based on marker_OTU from original balanced matrix
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

#start lineage with Arthropoda phylum
Zs$lineage <- apply( Zs[,44:49] , 1 , paste , collapse = ";" )

# make sure reads are integers
Zs$reads<-as.numeric(factor(Zs$reads))

# make abundance matrix with ESVs in rows and sites in columns
pivot2<-dcast(Zs, marker_OTU + lineage ~ GRDIname, value.var="reads", fun.aggregate = sum)

# Narrow down this abundance matrix to just the indicator ESVs
indic2<-merge(x=indic, y=pivot2, by="marker_OTU", all.x=TRUE)

# Remove unneeded columns
indic3<-indic2[,-c(2:8)]

# remove duplicate rows
indic4<-unique(indic3)

# make all numeric (double) columns integer instead
indic4[,3:length(names(indic4))] <- lapply(indic4[,3:length(names(indic4))], as.integer)

# set any NAs to zero
indic4[is.na(indic4)] <- 0

#reset rownames
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
#ensure esvs are all integers
esvs2[]<-lapply(esvs2, as.integer)
# Remove columns that sum to zero
esvs2<-esvs2[,colSums(esvs2) !=0]
#remove rows that sum to zero
esvs2<-esvs2[rowSums(esvs2) !=0,]

# merge
esvs<-merge(esvs1,esvs2,by="row.names")
esvs<-esvs[,-1]

# create samples file with columns sample, expt, site
## grab just GRDIname
Y<-data.frame(colnames(esvs))
Y2<-data.frame(Y[-c(1:2),])

# rename column
names(Y2)<-"sample_ID"

## split copy into new dataframe
Y3 <- data.frame(do.call('rbind', strsplit(as.character(Y2$sample_ID),'_',fixed=TRUE)),stringsAsFactors = FALSE)

## rename columns
names(Y3)<-c("site","date","coreExpt","rep")

## copy coreExpt column
Y3$expt<-Y3$coreExpt

# merge
Y4<-merge(Y2, Y3, by="row.names")
Y4<-Y4[,-c(1,4:6)]

## add column for layer
Y4$layer<-sapply(as.character(Y4[,1]), function (x) substrRight(x, 1))

## create site_layer column
Y4$site_layer<-paste(Y4$site, Y4$layer, sep="_")

# rename for readability
samples<-Y4

# make sure columns in esvs are in same order as rows in samples
# split up esvs to grab target vector
esvs2<-esvs[,3:length(names(esvs))]
target<-names(esvs2)
 
# match up row order in samples
samples<-samples[match(target, samples$sample_ID),]

# reset rownames
rownames(samples) <- NULL

# ensure samples are not factors
samples[]<-lapply(samples, as.character)

##### run metacoder
# check format of input files first
str(esvs)
str(samples)

# Add X in front of each expt, because starting with number messes everything up
samples$expt <- paste("X", samples$expt, sep="")

# parse taxonomic info and find reads per ESV per sample, SLOW be PATIENT!
obj <- parse_tax_data(esvs, 
                      class_cols = "lineage", 
                      class_sep = ";", 
                      class_regex = "^(.*)__(.*)$", 
                      class_key = c(tax_rank = "info", 
                                    tax_name = "taxon_name"))

print(obj)

# Get reads per taxon per sample
obj$data$tax_abund <- calc_taxon_abund(obj, "tax_data",
                                      cols = samples$sample_ID)
print(obj)

# number of samples that have reads for each site
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", 
                                   groups = samples$site)
print(obj)

## Create custom palettes

coral_palette <- function ()
{
  return(c("white","gray87","gray48","#F8766D"))
}

cyan_palette <- function () 
{
  return(c("white","gray87","gray48","#00BFC4"))
}

# plot taxonomic data for each sample in heat trees
set.seed(1) # This makes the plot appear the same each time it is run 
ILC<-heat_tree(obj, 
          title = "Island Lake", 
          title_size = 0.05,
          node_label = ifelse(obj$data$tax_occ$ILC45 <=42,"",taxon_names), # present in 50%+ samples n=84
          node_size = n_obs, #default is constant size
          node_label_size_range = c(0.02,0.05), #default c(0.02,0.02)
          node_color = ILC45, 
          node_color_range = coral_palette(), #default is diverging_palette
          node_size_axis_label = "ESVs",
          node_color_axis_label = "Samples",
          layout = "davidson-harel", # The primary layout algorithm
          initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations

set.seed(1) # This makes the plot appear the same each time it is run 
NZC<-heat_tree(obj, 
               title = "Nimitz",
               title_size = 0.05,
               node_label = ifelse(obj$data$tax_occ$NZC85 <=41.5,"",taxon_names), # present in 50%+ samples (n=83)
               node_size = n_obs, #default is constant size
               node_label_size_range = c(0.02,0.05), #default c(0.02,0.02)
               node_color = NZC85, 
               node_color_range = cyan_palette(), #default is diverging_palette
               node_size_axis_label = "ESVs",
               node_color_axis_label = "Samples",
               layout = "davidson-harel", # The primary layout algorithm
               initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations

# print both plots on one page
g<-grid.arrange(ILC, NZC, ncol=1)
ggsave("Fig3_two_heattrees.pdf", g)

# number of samples that have reads for each site
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", 
                                   groups = samples$site_layer)
print(obj)

## Create custom palettes

green_palette <- function ()
{
  return(c("white","gray87","gray48","#4DAF4A"))
}

blue_palette <- function () 
{
  return(c("white","gray87","gray48","#377EB8"))
}

orange_palette <- function()
{
  return(c("white","gray87","gray48","#FF7F00"))
}

par(mar=c(5.1,4.1,4.1,2.1))

# plot taxonomic data for each sample in heat trees
set.seed(1) # This makes the plot appear the same each time it is run 
I_B<-heat_tree(obj, 
               title = "Island Lake - Bryophyte", 
               title_size = 0.05,
               node_label = ifelse(obj$data$tax_occ$ILC45_B <=14,"",taxon_names), # present in 50%+ samples n=28
               node_size = n_obs, #default is constant size
               node_label_size_range = c(0.02,0.05), #default c(0.02,0.02)
               node_color = ILC45_B, 
               node_color_range = green_palette(), #default is diverging_palette
               node_size_axis_label = "ESVs",
               node_color_axis_label = "Samples",
               layout = "davidson-harel", # The primary layout algorithm
               initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations

set.seed(1) # This makes the plot appear the same each time it is run 
I_O<-heat_tree(obj, 
               title = "Island Lake - Organic",
               title_size = 0.05,
               node_label = ifelse(obj$data$tax_occ$ILC45_O <=14,"",taxon_names), # present in 50%+ samples (n=28)
               node_size = n_obs, #default is constant size
               node_label_size_range = c(0.02,0.05), #default c(0.02,0.02)
               node_color = ILC45_O, 
               node_color_range = blue_palette(), #default is diverging_palette
               node_size_axis_label = "ESVs",
               node_color_axis_label = "Samples",
               layout = "davidson-harel", # The primary layout algorithm
               initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations

set.seed(1) # This makes the plot appear the same each time it is run 
I_M<-heat_tree(obj, 
               title = "Island Lake - Mineral",
               title_size = 0.05,
               node_label = ifelse(obj$data$tax_occ$ILC45_M <=14,"",taxon_names), # present in 50%+ samples (n=28)
               node_size = n_obs, #default is constant size
               node_label_size_range = c(0.02,0.05), #default c(0.02,0.02)
               node_color = ILC45_M, 
               node_color_range = orange_palette(), #default is diverging_palette
               node_size_axis_label = "ESVs",
               node_color_axis_label = "Samples",
               layout = "davidson-harel", # The primary layout algorithm
               initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations

# plot taxonomic data for each sample in heat trees
set.seed(1) # This makes the plot appear the same each time it is run 
N_B<-heat_tree(obj, 
               title = "Nimitz - Bryophyte", 
               title_size = 0.05,
               node_label = ifelse(obj$data$tax_occ$NZC85_B <=14,"",taxon_names), # present in 50%+ samples n=28
               node_size = n_obs, #default is constant size
               node_label_size_range = c(0.02,0.05), #default c(0.02,0.02)
               node_color = NZC85_B, 
               node_color_range = green_palette(), #default is diverging_palette
               node_size_axis_label = "ESVs",
               node_color_axis_label = "Samples",
               layout = "davidson-harel", # The primary layout algorithm
               initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations

set.seed(1) # This makes the plot appear the same each time it is run 
N_O<-heat_tree(obj, 
               title = "Nimitz - Organic",
               title_size = 0.05,
               node_label = ifelse(obj$data$tax_occ$NZC85_O <=13.5,"",taxon_names), # present in 50%+ samples (n=28)
               node_size = n_obs, #default is constant size
               node_label_size_range = c(0.02,0.05), #default c(0.02,0.02)
               node_color = NZC85_O, 
               node_color_range = blue_palette(), #default is diverging_palette
               node_size_axis_label = "ESVs",
               node_color_axis_label = "Samples",
               layout = "davidson-harel", # The primary layout algorithm
               initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations

set.seed(1) # This makes the plot appear the same each time it is run 
N_M<-heat_tree(obj, 
               title = "Nimitz - Mineral",
               title_size = 0.05,
               node_label = ifelse(obj$data$tax_occ$NZC85_M <=14,"",taxon_names), # present in 50%+ samples (n=28)
               node_size = n_obs, #default is constant size
               node_label_size_range = c(0.02,0.05), #default c(0.02,0.02)
               node_color = NZC85_M, 
               node_color_range = orange_palette(), #default is diverging_palette
               node_size_axis_label = "ESVs",
               node_color_axis_label = "Samples",
               layout = "davidson-harel", # The primary layout algorithm
               initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations

# print both plots on one page
h<-grid.arrange(I_B, N_B, I_O, N_O, I_M, N_M, ncol=2)
ggsave("FigS10_layers_heattrees.pdf", h, height=10, width=7.5, units = "in")

# number of samples that have reads for each extraction experiment
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", 
                                   groups = samples$expt)
print(obj)

# plot taxonomic data for each sample in heat trees

grey_palette <- function()
{
  return(c("white","gray87","gray48","red4"))
}

set.seed(1) # This makes the plot appear the same each time it is run 
E1<-heat_tree(obj, 
               title = "1 DNA extraction", 
               title_size = 0.05,
               node_label = ifelse(obj$data$tax_occ[[2]] <=12,"",taxon_names), # present in 50%+ samples n=24
               node_size = n_obs, #default is constant size
               node_label_size_range = c(0.02,0.05), #default c(0.02,0.02)
               node_color = obj$data$tax_occ[[2]], 
               node_color_range = grey_palette(), #default is diverging_palette
               node_size_axis_label = "ESVs",
               node_color_axis_label = "Samples",
               layout = "davidson-harel", # The primary layout algorithm
               initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations

set.seed(1) # This makes the plot appear the same each time it is run 
E3<-heat_tree(obj, 
               title = "3 DNA extractions",
               title_size = 0.05,
               node_label = ifelse(obj$data$tax_occ[[3]] <=12,"",taxon_names), # present in 50%+ samples (n=23)
               node_size = n_obs, #default is constant size
               node_label_size_range = c(0.02,0.05), #default c(0.02,0.02)
               node_color = obj$data$tax_occ[[3]], 
               node_color_range = grey_palette(), #default is diverging_palette
               node_size_axis_label = "ESVs",
               node_color_axis_label = "Samples",
               layout = "davidson-harel", # The primary layout algorithm
               initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations

# print both plots on one page
i<-grid.arrange(E1, E3, ncol=1)
ggsave("FigS9_two_extraction_heattrees.pdf", i)

# number of samples that have reads for each site
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", 
                                   groups = samples$expt)
print(obj)

# plot taxonomic data for each sample in heat trees

set.seed(1) # This makes the plot appear the same each time it is run 
C1<-heat_tree(obj, 
               title = "1 core", 
               title_size = 0.05,
               node_label = ifelse(obj$data$tax_occ[[3]] <=11.5,"",taxon_names), # present in 50%+ samples n=24
               node_size = n_obs, #default is constant size
               node_label_size_range = c(0.02,0.05), #default c(0.02,0.02)
               node_color = obj$data$tax_occ[[3]], 
               node_color_range = grey_palette(), #default is diverging_palette
               node_size_axis_label = "ESVs",
               node_color_axis_label = "Samples",
               layout = "davidson-harel", # The primary layout algorithm
               initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations

set.seed(1) # This makes the plot appear the same each time it is run 
C2<-heat_tree(obj, 
               title = "2 pooled cores",
               title_size = 0.05,
               node_label = ifelse(obj$data$tax_occ[[4]] <=12,"",taxon_names), # present in 50%+ samples (n=24)
               node_size = n_obs, #default is constant size
               node_label_size_range = c(0.02,0.05), #default c(0.02,0.02)
               node_color = obj$data$tax_occ[[4]], 
               node_color_range = grey_palette(), #default is diverging_palette
               node_size_axis_label = "ESVs",
               node_color_axis_label = "Samples",
               layout = "davidson-harel", # The primary layout algorithm
               initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations

set.seed(1) # This makes the plot appear the same each time it is run 
C4<-heat_tree(obj, 
               title = "4 pooled cores",
               title_size = 0.05,
               node_label = ifelse(obj$data$tax_occ[[5]] <=12,"",taxon_names), # present in 50%+ samples (n=24)
               node_size = n_obs, #default is constant size
               node_label_size_range = c(0.02,0.05), #default c(0.02,0.02)
               node_color = obj$data$tax_occ[[5]], 
               node_color_range = grey_palette(), #default is diverging_palette
               node_size_axis_label = "ESVs",
               node_color_axis_label = "Samples",
               layout = "davidson-harel", # The primary layout algorithm
               initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations

# plot taxonomic data for each sample in heat trees
set.seed(1) # This makes the plot appear the same each time it is run 
C6<-heat_tree(obj, 
               title = "6 pooled cores", 
               title_size = 0.05,
               node_label = ifelse(obj$data$tax_occ[[6]] <=12,"",taxon_names), # present in 50%+ samples n=24
               node_size = n_obs, #default is constant size
               node_label_size_range = c(0.02,0.05), #default c(0.02,0.02)
               node_color = obj$data$tax_occ[[6]], 
               node_color_range = grey_palette(), #default is diverging_palette
               node_size_axis_label = "ESVs",
               node_color_axis_label = "Samples",
               layout = "davidson-harel", # The primary layout algorithm
               initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations

set.seed(1) # This makes the plot appear the same each time it is run 
C8<-heat_tree(obj, 
               title = "8 pooled cores",
               title_size = 0.05,
               node_label = ifelse(obj$data$tax_occ[[7]] <=12,"",taxon_names), # present in 10%+ samples (n=24)
               node_size = n_obs, #default is constant size
               node_label_size_range = c(0.02,0.05), #default c(0.02,0.02)
               node_color = obj$data$tax_occ[[7]], 
               node_color_range = grey_palette(), #default is diverging_palette
               node_size_axis_label = "ESVs",
               node_color_axis_label = "Samples",
               layout = "davidson-harel", # The primary layout algorithm
               initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations

set.seed(1) # This makes the plot appear the same each time it is run 
C915<-heat_tree(obj, 
               title = "9-15 cores pooled",
               title_size = 0.05,
               node_label = ifelse(obj$data$tax_occ[[8]] <=12,"",taxon_names), # present in 10%+ samples (n=24)
               node_size = n_obs, #default is constant size
               node_label_size_range = c(0.02,0.05), #default c(0.02,0.02)
               node_color = obj$data$tax_occ[[8]], 
               node_color_range = grey_palette(), #default is diverging_palette
               node_size_axis_label = "ESVs",
               node_color_axis_label = "Samples",
               layout = "davidson-harel", # The primary layout algorithm
               initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations

# print both plots on one page
j<-grid.arrange(C1, C2, C4, C6, C8, C915, ncol=2)
ggsave("FigS8_pooled_heattrees.pdf", j, height=10, width=7.5, units = "in")
