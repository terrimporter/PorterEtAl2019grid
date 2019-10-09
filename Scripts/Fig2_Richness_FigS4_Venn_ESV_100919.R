# Teresita M. Porter, Oct. 9, 2019

library(vegan)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(grid)
library(gridExtra)
library(reshape2)
library(dplyr)
library("car")
# to calculate venn counts
source("http://www.bioconductor.org/biocLite.R")
biocLite("limma")
library(limma)
library(ggforce) # to draw circles with ggplot

########################################################
##### NEW FUNCTIUON TO GET TARGET GRIDCOORD DF's
########################################################

# New function to grab target gridcoord rows
# Retrieve target gridcoord rows for each of 4 samples (reps in columns)

get_target_gridcoord<- function (matrix, df) {
  
  out=NULL
  df.list<-c()
  for (i in 1:(ncol(matrix))) {
    gridcoord.list<-matrix[,i]
    
    for (j in 1:length(gridcoord.list)) {
      targetgridcoord<-gridcoord.list[[j]] 
      df.j<-df[df$gridcoord==targetgridcoord,]
   #   print(i)
      df.j$sample<-i
      out=rbind(out, df.j)
    }
    df.list[[i]]<-out
    out<-data.frame()
  }
  return(df.list)
}

########################################################

# edit cores column to reflect bioinformatic pooling
edit_cores<-function(list, x){
  list2<-list()
  df<-data.frame()
  for (i in 1:length(list)) { # loop through df in list
    df<-list[[i]]
    df$cores<-x
    list2[[i]]<-df # add to new list
    df<-data.frame()
  }
  return(list2)
}

########################################################

#read infile prepared by python script
master<-read.table(file="LV2016_2.csv", head=TRUE, sep=",")

# Select phylum Arthropoda only (use this to get overall richness across all expts)
Arth_df<-master[master$phylum=="Arthropoda",]

# OTU matrix for 1C1E experiment (n=8 x 3 layers)
#A<-Arth_df[Arth_df$extractions=="1",]

###############################################################
# Create matrices for Fig 1A, bioinformatic pooling 1-15 cores
# OTU matrix for 1C3E experiment ( n=36 x 3 layers)
B<-Arth_df[(Arth_df$cores=="1") & (Arth_df$extractions=="3"),]

# randomly sample 1 gridcoord from ILC and NZC sites, repeat this sampling 4 times, create matrix for rarefaction
B.ILC.gridcoord.1<-t(as.matrix(sapply(rep(1,4), function (x) sample(unique(B$gridcoord[B$site=="ILC45"]), x, replace=FALSE))))
B.ILC.gridcoord.2<-sapply(rep(2,4), function (x) sample(unique(B$gridcoord[B$site=="ILC45"]), x, replace=FALSE))
B.ILC.gridcoord.4<-sapply(rep(4,4), function (x) sample(unique(B$gridcoord[B$site=="ILC45"]), x, replace=FALSE))
B.ILC.gridcoord.6<-sapply(rep(6,4), function (x) sample(unique(B$gridcoord[B$site=="ILC45"]), x, replace=FALSE))
B.ILC.gridcoord.8<-sapply(rep(8,4), function (x) sample(unique(B$gridcoord[B$site=="ILC45"]), x, replace=FALSE))
B.ILC.gridcoord.15<-sapply(rep(15,4), function (x) sample(unique(B$gridcoord[B$site=="ILC45"]), x, replace=FALSE))
B.NZC.gridcoord.1<-t(as.matrix(sapply(rep(1,4), function (x) sample(unique(B$gridcoord[B$site=="NZC85"]), x, replace=FALSE))))
B.NZC.gridcoord.2<-sapply(rep(2,4), function (x) sample(unique(B$gridcoord[B$site=="NZC85"]), x, replace=FALSE))
B.NZC.gridcoord.4<-sapply(rep(4,4), function (x) sample(unique(B$gridcoord[B$site=="NZC85"]), x, replace=FALSE))
B.NZC.gridcoord.6<-sapply(rep(6,4), function (x) sample(unique(B$gridcoord[B$site=="NZC85"]), x, replace=FALSE))
B.NZC.gridcoord.8<-sapply(rep(8,4), function (x) sample(unique(B$gridcoord[B$site=="NZC85"]), x, replace=FALSE))
B.NZC.gridcoord.15<-sapply(rep(15,4), function (x) sample(unique(B$gridcoord[B$site=="NZC85"]), x, replace=FALSE))

# Retrieve target gridcoord rows for each of 4 samples
B.ILC.samples.1<-get_target_gridcoord(B.ILC.gridcoord.1,B)
B.ILC.samples.2<-get_target_gridcoord(B.ILC.gridcoord.2,B)
B.ILC.samples.4<-get_target_gridcoord(B.ILC.gridcoord.4,B)
B.ILC.samples.6<-get_target_gridcoord(B.ILC.gridcoord.6,B)
B.ILC.samples.8<-get_target_gridcoord(B.ILC.gridcoord.8,B)
B.ILC.samples.15<-get_target_gridcoord(B.ILC.gridcoord.15,B)
B.NZC.samples.1<-get_target_gridcoord(B.NZC.gridcoord.1,B)
B.NZC.samples.2<-get_target_gridcoord(B.NZC.gridcoord.2,B)
B.NZC.samples.4<-get_target_gridcoord(B.NZC.gridcoord.4,B)
B.NZC.samples.6<-get_target_gridcoord(B.NZC.gridcoord.6,B)
B.NZC.samples.8<-get_target_gridcoord(B.NZC.gridcoord.8,B)
B.NZC.samples.15<-get_target_gridcoord(B.NZC.gridcoord.15,B)

# edit cores field
B.ILC.samples.1b<-edit_cores(B.ILC.samples.1, 1)
B.ILC.samples.2b<-edit_cores(B.ILC.samples.2, 2)
B.ILC.samples.4b<-edit_cores(B.ILC.samples.4, 4)
B.ILC.samples.6b<-edit_cores(B.ILC.samples.6, 6)
B.ILC.samples.8b<-edit_cores(B.ILC.samples.8, 8)
B.ILC.samples.15b<-edit_cores(B.ILC.samples.15, 15)
B.NZC.samples.1b<-edit_cores(B.NZC.samples.1, 1)
B.NZC.samples.2b<-edit_cores(B.NZC.samples.2, 2)
B.NZC.samples.4b<-edit_cores(B.NZC.samples.4, 4)
B.NZC.samples.6b<-edit_cores(B.NZC.samples.6, 6)
B.NZC.samples.8b<-edit_cores(B.NZC.samples.8, 8)
B.NZC.samples.15b<-edit_cores(B.NZC.samples.15, 15)

# combine lists
B.combo<-c(B.ILC.samples.1b, B.ILC.samples.2b, B.ILC.samples.4b, B.ILC.samples.6b, B.ILC.samples.8b, B.ILC.samples.15b,
           B.NZC.samples.1b, B.NZC.samples.2b, B.NZC.samples.4b, B.NZC.samples.6b, B.NZC.samples.8b, B.NZC.samples.15b)

# combine each data frame in list into single data frame
B.combo.df<-do.call("rbind", B.combo)

###############################################################
# Create matrices for Fig 1C, compare 1 or 3 DNA extractions
# OTU matrix for 1C1E experiment (n=8 x 3 layers)
B2<-Arth_df[(Arth_df$cores=="1") & (Arth_df$extractions=="1"),]

# OTU matrix for 1C3E experiment subsampled down to 8 reps to compare with 1C3E for balanced design
B3<-Arth_df[(Arth_df$cores=="1") & (Arth_df$extractions=="3"),]

# match the 8 gridcoord already used in 1C1E
B3<-B3[B3$gridcoord=="11" |
         B3$gridcoord=="15" |
         B3$gridcoord=="23" |
         B3$gridcoord=="31" |
         B3$gridcoord=="35" |
         B3$gridcoord=="43" |
         B3$gridcoord=="51" |
         B3$gridcoord=="55",]

###############################################################
# Create matrices for Fig 1B, manual pooling 1-15 cores
# grab 1C3E df (n = 36 x 3 layers)
B

# Subsample down to 4 to make a fair comparison with XC3E
B.gridcoord<-sample(unique(B$gridcoord), 4, replace=FALSE)

# Retrieve gridcoord from B.gridcoord above (it's random so could change)
B<-B[B$gridcoord=="11" |
    B$gridcoord=="16" |
    B$gridcoord=="33" |
    B$gridcoord=="56",]

# OTU matrix for 2C3E experiment
C<-Arth_df[(Arth_df$cores=="2") & (Arth_df$extractions=="3"),]

# OTU matrix for 4C3E experiment
D<-Arth_df[(Arth_df$cores=="4") & (Arth_df$extractions=="3"),]

# OTU matrix for 6C3E experiment
E<-Arth_df[(Arth_df$cores=="6") & (Arth_df$extractions=="3"),]

# OTU matrix for 8C3E experiment
F<-Arth_df[(Arth_df$cores=="8") & (Arth_df$extractions=="3"),]

# OTU matrix for 15C3E experiment
G<-Arth_df[((Arth_df$cores=="9") |
            (Arth_df$cores=="12") | (Arth_df$cores=="13") | 
            (Arth_df$cores=="14") | (Arth_df$cores=="15"))
           & (Arth_df$extractions=="3"),]

#split out ILC and NZC datasets
ILC_Arth<-Arth_df[Arth_df$site=="ILC45",]
NZC_Arth<-Arth_df[Arth_df$site=="NZC85",]
ILC_B<-B[B$site=="ILC45",]
NZC_B<-B[B$site=="NZC85",]
ILC_B.combo.df<-B.combo.df[B.combo.df$site=="ILC45",]
NZC_B.combo.df<-B.combo.df[B.combo.df$site=="NZC85",]
ILC_B2<-B2[B2$site=="ILC45",]
NZC_B2<-B2[B2$site=="NZC85",]
ILC_B3<-B3[B3$site=="ILC45",]
NZC_B3<-B3[B3$site=="NZC85",]
ILC_C<-C[C$site=="ILC45",]
NZC_C<-C[C$site=="NZC85",]
ILC_D<-D[D$site=="ILC45",]
NZC_D<-D[D$site=="NZC85",]
ILC_E<-E[E$site=="ILC45",]
NZC_E<-E[E$site=="NZC85",]
ILC_F<-F[F$site=="ILC45",]
NZC_F<-F[F$site=="NZC85",]
ILC_G<-G[G$site=="ILC45",]
NZC_G<-G[G$site=="NZC85",]

# Pivot to make matrices for vegan, samples in rows, OTUs in columns, reads in cells
matrixArth_ILC<-dcast(ILC_Arth, GRDIname ~ marker_OTU, value.var="reads", fun.aggregate = sum)
matrixB_ILC<-dcast(ILC_B, GRDIname ~ marker_OTU, value.var="reads", fun.aggregate = sum)
matrixB.combo.df_ILC<-dcast(ILC_B.combo.df, paste(site,layer,cores,sample,sep="_") ~ marker_OTU, value.var="reads", fun.aggregate = sum)
names(matrixB.combo.df_ILC)[1]<-"GRDIname"
matrixB2_ILC<-dcast(ILC_B2, GRDIname ~ marker_OTU, value.var="reads", fun.aggregate = sum)
matrixB3_ILC<-dcast(ILC_B3, GRDIname ~ marker_OTU, value.var="reads", fun.aggregate = sum)
matrixC_ILC<-dcast(ILC_C, GRDIname ~ marker_OTU, value.var="reads", fun.aggregate = sum)
matrixD_ILC<-dcast(ILC_D, GRDIname ~ marker_OTU, value.var="reads", fun.aggregate = sum)
matrixE_ILC<-dcast(ILC_E, GRDIname ~ marker_OTU, value.var="reads", fun.aggregate = sum)
matrixF_ILC<-dcast(ILC_F, GRDIname ~ marker_OTU, value.var="reads", fun.aggregate = sum)
matrixG_ILC<-dcast(ILC_G, GRDIname ~ marker_OTU, value.var="reads", fun.aggregate = sum)

matrixArth_NZC<-dcast(NZC_Arth, GRDIname ~ marker_OTU, value.var="reads", fun.aggregate = sum)
matrixB_NZC<-dcast(NZC_B, GRDIname ~ marker_OTU, value.var="reads", fun.aggregate = sum)
matrixB.combo.df_NZC<-dcast(NZC_B.combo.df, paste(site,layer,cores,sample,sep="_") ~ marker_OTU, value.var="reads", fun.aggregate = sum)
names(matrixB.combo.df_NZC)[1]<-"GRDIname"
matrixB2_NZC<-dcast(NZC_B2, GRDIname ~ marker_OTU, value.var="reads", fun.aggregate = sum)
matrixB3_NZC<-dcast(NZC_B3, GRDIname ~ marker_OTU, value.var="reads", fun.aggregate = sum)
matrixC_NZC<-dcast(NZC_C, GRDIname ~ marker_OTU, value.var="reads", fun.aggregate = sum)
matrixD_NZC<-dcast(NZC_D, GRDIname ~ marker_OTU, value.var="reads", fun.aggregate = sum)
matrixE_NZC<-dcast(NZC_E, GRDIname ~ marker_OTU, value.var="reads", fun.aggregate = sum)
matrixF_NZC<-dcast(NZC_F, GRDIname ~ marker_OTU, value.var="reads", fun.aggregate = sum)
matrixG_NZC<-dcast(NZC_G, GRDIname ~ marker_OTU, value.var="reads", fun.aggregate = sum)

# Move GRDIname to row names
rownames(matrixArth_ILC)<-matrixArth_ILC$GRDIname
matrixArth_ILC<-matrixArth_ILC[,-1]
rownames(matrixB_ILC)<-matrixB_ILC$GRDIname
matrixB_ILC<-matrixB_ILC[,-1]
rownames(matrixB.combo.df_ILC)<-matrixB.combo.df_ILC$GRDIname
matrixB.combo.df_ILC<-matrixB.combo.df_ILC[,-1]
rownames(matrixB2_ILC)<-matrixB2_ILC$GRDIname
matrixB2_ILC<-matrixB2_ILC[,-1]
rownames(matrixB3_ILC)<-matrixB3_ILC$GRDIname
matrixB3_ILC<-matrixB3_ILC[,-1]
rownames(matrixC_ILC)<-matrixC_ILC$GRDIname
matrixC_ILC<-matrixC_ILC[,-1]
rownames(matrixD_ILC)<-matrixD_ILC$GRDIname
matrixD_ILC<-matrixD_ILC[,-1]
rownames(matrixE_ILC)<-matrixE_ILC$GRDIname
matrixE_ILC<-matrixE_ILC[,-1]
rownames(matrixF_ILC)<-matrixF_ILC$GRDIname
matrixF_ILC<-matrixF_ILC[,-1]
rownames(matrixG_ILC)<-matrixG_ILC$GRDIname
matrixG_ILC<-matrixG_ILC[,-1]
rownames(matrixArth_NZC)<-matrixArth_NZC$GRDIname
matrixArth_NZC<-matrixArth_NZC[,-1]
rownames(matrixB_NZC)<-matrixB_NZC$GRDIname
matrixB_NZC<-matrixB_NZC[,-1]
rownames(matrixB.combo.df_NZC)<-matrixB.combo.df_NZC$GRDIname
matrixB.combo.df_NZC<-matrixB.combo.df_NZC[,-1]
rownames(matrixB2_NZC)<-matrixB2_NZC$GRDIname
matrixB2_NZC<-matrixB2_NZC[,-1]
rownames(matrixB3_NZC)<-matrixB3_NZC$GRDIname
matrixB3_NZC<-matrixB3_NZC[,-1]
rownames(matrixC_NZC)<-matrixC_NZC$GRDIname
matrixC_NZC<-matrixC_NZC[,-1]
rownames(matrixD_NZC)<-matrixD_NZC$GRDIname
matrixD_NZC<-matrixD_NZC[,-1]
rownames(matrixE_NZC)<-matrixE_NZC$GRDIname
matrixE_NZC<-matrixE_NZC[,-1]
rownames(matrixF_NZC)<-matrixF_NZC$GRDIname
matrixF_NZC<-matrixF_NZC[,-1]
rownames(matrixG_NZC)<-matrixG_NZC$GRDIname
matrixG_NZC<-matrixG_NZC[,-1]

#remove columns with only zeros
ILC_Arth_notnull<-matrixArth_ILC[,colSums(matrixArth_ILC) !=0]
ILC_B_notnull<-matrixB_ILC[,colSums(matrixB_ILC) !=0]
ILC_B.combo.df_notnull<-matrixB.combo.df_ILC[,colSums(matrixB.combo.df_ILC) !=0]
ILC_B2_notnull<-matrixB2_ILC[,colSums(matrixB2_ILC) !=0]
ILC_B3_notnull<-matrixB3_ILC[,colSums(matrixB3_ILC) !=0]
ILC_C_notnull<-matrixC_ILC[,colSums(matrixC_ILC) !=0]
ILC_D_notnull<-matrixD_ILC[,colSums(matrixD_ILC) !=0]
ILC_E_notnull<-matrixE_ILC[,colSums(matrixE_ILC) !=0]
ILC_F_notnull<-matrixF_ILC[,colSums(matrixF_ILC) !=0]
ILC_G_notnull<-matrixG_ILC[,colSums(matrixG_ILC) !=0]

NZC_Arth_notnull<-matrixArth_NZC[,colSums(matrixArth_NZC) !=0]
NZC_B_notnull<-matrixB_NZC[,colSums(matrixB_NZC) !=0]
NZC_B.combo.df_notnull<-matrixB.combo.df_NZC[,colSums(matrixB.combo.df_NZC) !=0]
NZC_B2_notnull<-matrixB2_NZC[,colSums(matrixB2_NZC) !=0]
NZC_B3_notnull<-matrixB3_NZC[,colSums(matrixB3_NZC) !=0]
NZC_C_notnull<-matrixC_NZC[,colSums(matrixC_NZC) !=0]
NZC_D_notnull<-matrixD_NZC[,colSums(matrixD_NZC) !=0]
NZC_E_notnull<-matrixE_NZC[,colSums(matrixE_NZC) !=0]
NZC_F_notnull<-matrixF_NZC[,colSums(matrixF_NZC) !=0]
NZC_G_notnull<-matrixG_NZC[,colSums(matrixG_NZC) !=0]

#remove rows with only zeros
ILC_Arth_notnull2<-ILC_Arth_notnull[rowSums(ILC_Arth_notnull) !=0,]
ILC_B_notnull2<-ILC_B_notnull[rowSums(ILC_B_notnull) !=0,]
ILC_B.combo.df_notnull2<-ILC_B.combo.df_notnull[rowSums(ILC_B.combo.df_notnull) !=0,]
ILC_B2_notnull2<-ILC_B2_notnull[rowSums(ILC_B2_notnull) !=0,]
ILC_B3_notnull2<-ILC_B3_notnull[rowSums(ILC_B3_notnull) !=0,]
ILC_C_notnull2<-ILC_C_notnull[rowSums(ILC_C_notnull) !=0,]
ILC_D_notnull2<-ILC_D_notnull[rowSums(ILC_D_notnull) !=0,]
ILC_E_notnull2<-ILC_E_notnull[rowSums(ILC_E_notnull) !=0,]
ILC_F_notnull2<-ILC_F_notnull[rowSums(ILC_F_notnull) !=0,]
ILC_G_notnull2<-ILC_G_notnull[rowSums(ILC_G_notnull) !=0,]

NZC_Arth_notnull2<-NZC_Arth_notnull[rowSums(NZC_Arth_notnull) !=0,]
NZC_B_notnull2<-NZC_B_notnull[rowSums(NZC_B_notnull) !=0,]
NZC_B.combo.df_notnull2<-NZC_B.combo.df_notnull[rowSums(NZC_B.combo.df_notnull) !=0,]
NZC_B2_notnull2<-NZC_B2_notnull[rowSums(NZC_B2_notnull) !=0,]
NZC_B3_notnull2<-NZC_B3_notnull[rowSums(NZC_B3_notnull) !=0,]
NZC_C_notnull2<-NZC_C_notnull[rowSums(NZC_C_notnull) !=0,]
NZC_D_notnull2<-NZC_D_notnull[rowSums(NZC_D_notnull) !=0,]
NZC_E_notnull2<-NZC_E_notnull[rowSums(NZC_E_notnull) !=0,]
NZC_F_notnull2<-NZC_F_notnull[rowSums(NZC_F_notnull) !=0,]
NZC_G_notnull2<-NZC_G_notnull[rowSums(NZC_G_notnull) !=0,]

#calculate 15th percentile for rrarefy function
ILC_Arth_15percentile<-quantile(rowSums(ILC_Arth_notnull2), prob=0.15)
ILC_B_15percentile<-quantile(rowSums(ILC_B_notnull2), prob=0.15)
ILC_B.combo.df_15percentile<-quantile(rowSums(ILC_B.combo.df_notnull2), prob=0.15)
ILC_B2_15percentile<-quantile(rowSums(ILC_B2_notnull2), prob=0.15)
ILC_B3_15percentile<-quantile(rowSums(ILC_B3_notnull2), prob=0.15)
ILC_C_15percentile<-quantile(rowSums(ILC_C_notnull2), prob=0.15)
ILC_D_15percentile<-quantile(rowSums(ILC_D_notnull2), prob=0.15)
ILC_E_15percentile<-quantile(rowSums(ILC_E_notnull2), prob=0.15)
ILC_F_15percentile<-quantile(rowSums(ILC_F_notnull2), prob=0.15)
ILC_G_15percentile<-quantile(rowSums(ILC_G_notnull2), prob=0.15)

NZC_Arth_15percentile<-quantile(rowSums(NZC_Arth_notnull2), prob=0.15)
NZC_B_15percentile<-quantile(rowSums(NZC_B_notnull2), prob=0.15)
NZC_B.combo.df_15percentile<-quantile(rowSums(NZC_B.combo.df_notnull2), prob=0.15)
NZC_B2_15percentile<-quantile(rowSums(NZC_B2_notnull2), prob=0.15)
NZC_B3_15percentile<-quantile(rowSums(NZC_B3_notnull2), prob=0.15)
NZC_C_15percentile<-quantile(rowSums(NZC_C_notnull2), prob=0.15)
NZC_D_15percentile<-quantile(rowSums(NZC_D_notnull2), prob=0.15)
NZC_E_15percentile<-quantile(rowSums(NZC_E_notnull2), prob=0.15)
NZC_F_15percentile<-quantile(rowSums(NZC_F_notnull2), prob=0.15)
NZC_G_15percentile<-quantile(rowSums(NZC_G_notnull2), prob=0.15)

set.seed(1234)

###################################################################
##### Rarefy the dataset down to the 15th percentile
###################################################################

### Only do this for OTUs becaue too much missing data for genus and family ranks ###
#Rarefy original ILC matrices down to 15th percentile library size to normalize read depth across samples
# sites in rows, OTUs in columns, reads in cells
ILC_Arth_df<-rrarefy(ILC_Arth_notnull2,sample=ILC_Arth_15percentile)
ILC_B_df<-rrarefy(ILC_B_notnull2,sample=ILC_B_15percentile)
ILC_B.combo.df_df<-rrarefy(ILC_B.combo.df_notnull2,sample=ILC_B.combo.df_15percentile)
ILC_B2_df<-rrarefy(ILC_B2_notnull2,sample=ILC_B2_15percentile)
ILC_B3_df<-rrarefy(ILC_B3_notnull2,sample=ILC_B3_15percentile)
ILC_C_df<-rrarefy(ILC_C_notnull2,sample=ILC_C_15percentile)
ILC_D_df<-rrarefy(ILC_D_notnull2,sample=ILC_D_15percentile)
ILC_E_df<-rrarefy(ILC_E_notnull2,sample=ILC_E_15percentile)
ILC_F_df<-rrarefy(ILC_F_notnull2,sample=ILC_F_15percentile)
ILC_G_df<-rrarefy(ILC_G_notnull2,sample=ILC_G_15percentile)

NZC_Arth_df<-rrarefy(NZC_Arth_notnull2,sample=NZC_Arth_15percentile)
NZC_B_df<-rrarefy(NZC_B_notnull2,sample=NZC_B_15percentile)
NZC_B.combo.df_df<-rrarefy(NZC_B.combo.df_notnull2,sample=NZC_B.combo.df_15percentile)
NZC_B2_df<-rrarefy(NZC_B2_notnull2,sample=NZC_B2_15percentile)
NZC_B3_df<-rrarefy(NZC_B3_notnull2,sample=NZC_B3_15percentile)
NZC_C_df<-rrarefy(NZC_C_notnull2,sample=NZC_C_15percentile)
NZC_D_df<-rrarefy(NZC_D_notnull2,sample=NZC_D_15percentile)
NZC_E_df<-rrarefy(NZC_E_notnull2,sample=NZC_E_15percentile)
NZC_F_df<-rrarefy(NZC_F_notnull2,sample=NZC_F_15percentile)
NZC_G_df<-rrarefy(NZC_G_notnull2,sample=NZC_G_15percentile)

##################################################################################
# Plot richness for 1C3E bioinformatically pooled cores (not accumulation curves)
##################################################################################

# do specnum for each site and layer separately
ILC_B.combo.df_specnumber_B<-data.frame(specnumber(ILC_B.combo.df_df[grepl("B_", rownames(ILC_B.combo.df_df)),]))
ILC_B.combo.df_specnumber_O<-data.frame(specnumber(ILC_B.combo.df_df[grepl("O_", rownames(ILC_B.combo.df_df)),]))
ILC_B.combo.df_specnumber_M<-data.frame(specnumber(ILC_B.combo.df_df[grepl("M_", rownames(ILC_B.combo.df_df)),]))
NZC_B.combo.df_specnumber_B<-data.frame(specnumber(NZC_B.combo.df_df[grepl("B_", rownames(NZC_B.combo.df_df)),]))
NZC_B.combo.df_specnumber_O<-data.frame(specnumber(NZC_B.combo.df_df[grepl("O_", rownames(NZC_B.combo.df_df)),]))
NZC_B.combo.df_specnumber_M<-data.frame(specnumber(NZC_B.combo.df_df[grepl("M_", rownames(NZC_B.combo.df_df)),]))

# Reformat lists for ggplot
ILC_B.combo.df_specnumber_B[,2]<-"ILC"
ILC_B.combo.df_specnumber_B[,3]<-"Bryophyte"
foo <- data.frame(do.call('rbind', strsplit(as.character(rownames(ILC_B.combo.df_specnumber_B)),'_',fixed=TRUE)))
ILC_B.combo.df_specnumber_B[,4]<-foo$X3
ILC_B.combo.df_specnumber_B[,5]<-foo$X4
names(ILC_B.combo.df_specnumber_B)<-c("richness","site","layer","cores","sample")

ILC_B.combo.df_specnumber_O[,2]<-"ILC"
ILC_B.combo.df_specnumber_O[,3]<-"Organic"
foo <- data.frame(do.call('rbind', strsplit(as.character(rownames(ILC_B.combo.df_specnumber_O)),'_',fixed=TRUE)))
ILC_B.combo.df_specnumber_O[,4]<-foo$X3
ILC_B.combo.df_specnumber_O[,5]<-foo$X4
names(ILC_B.combo.df_specnumber_O)<-c("richness","site","layer","cores","sample")

ILC_B.combo.df_specnumber_M[,2]<-"ILC"
ILC_B.combo.df_specnumber_M[,3]<-"Mineral"
foo <- data.frame(do.call('rbind', strsplit(as.character(rownames(ILC_B.combo.df_specnumber_M)),'_',fixed=TRUE)))
ILC_B.combo.df_specnumber_M[,4]<-foo$X3
ILC_B.combo.df_specnumber_M[,5]<-foo$X4
names(ILC_B.combo.df_specnumber_M)<-c("richness","site","layer","cores","sample")

NZC_B.combo.df_specnumber_B[,2]<-"NZC"
NZC_B.combo.df_specnumber_B[,3]<-"Bryophyte"
foo <- data.frame(do.call('rbind', strsplit(as.character(rownames(NZC_B.combo.df_specnumber_B)),'_',fixed=TRUE)))
NZC_B.combo.df_specnumber_B[,4]<-foo$X3
NZC_B.combo.df_specnumber_B[,5]<-foo$X4
names(NZC_B.combo.df_specnumber_B)<-c("richness","site","layer","cores","sample")

NZC_B.combo.df_specnumber_O[,2]<-"NZC"
NZC_B.combo.df_specnumber_O[,3]<-"Organic"
foo <- data.frame(do.call('rbind', strsplit(as.character(rownames(NZC_B.combo.df_specnumber_O)),'_',fixed=TRUE)))
NZC_B.combo.df_specnumber_O[,4]<-foo$X3
NZC_B.combo.df_specnumber_O[,5]<-foo$X4
names(NZC_B.combo.df_specnumber_O)<-c("richness","site","layer","cores","sample")

NZC_B.combo.df_specnumber_M[,2]<-"NZC"
NZC_B.combo.df_specnumber_M[,3]<-"Mineral"
foo <- data.frame(do.call('rbind', strsplit(as.character(rownames(NZC_B.combo.df_specnumber_M)),'_',fixed=TRUE)))
NZC_B.combo.df_specnumber_M[,4]<-foo$X3
NZC_B.combo.df_specnumber_M[,5]<-foo$X4
names(NZC_B.combo.df_specnumber_M)<-c("richness","site","layer","cores","sample")

summary_bio<-rbind(ILC_B.combo.df_specnumber_B, ILC_B.combo.df_specnumber_O, ILC_B.combo.df_specnumber_M, 
                NZC_B.combo.df_specnumber_B, NZC_B.combo.df_specnumber_O, NZC_B.combo.df_specnumber_M)

# Create factors
summary_bio$layer<-factor(summary_bio$layer, levels=c("Bryophyte", "Organic", "Mineral"))
summary_bio$site<-factor(summary_bio$site, levels=c("ILC", "NZC"),
                         labels=c("Island Lake","Nimitz"))
summary_bio$cores<-factor(summary_bio$cores, levels=c("1","2","4","6","8","15"))

# create boxplot plot for 1C3E, bioinformatically pooled cores
p1<-ggplot(summary_bio, aes(x=cores, y=richness, fill=layer)) +
  ggtitle("A)") +
  facet_wrap(~site+layer) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  labs(x="Bioinformatically pooled individual samples",y="ESV Richness") +
  scale_fill_manual(values=c("#4DAF4A","#377EB8","#FF7F00")) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(), 
        legend.position = "none")

###################################################################
##### Calculate species accumulation for 1C3E & XC3E
# pooled cores vs OTU richness
###################################################################

# do specnum for each site and layer separately
ILC_B_specnumber_B<-data.frame(specnumber(ILC_B_df[grepl("B$", rownames(ILC_B_df)),]))
ILC_B_specnumber_O<-data.frame(specnumber(ILC_B_df[grepl("O$", rownames(ILC_B_df)),]))
ILC_B_specnumber_M<-data.frame(specnumber(ILC_B_df[grepl("M$", rownames(ILC_B_df)),]))
NZC_B_specnumber_B<-data.frame(specnumber(NZC_B_df[grepl("B$", rownames(NZC_B_df)),]))
NZC_B_specnumber_O<-data.frame(specnumber(NZC_B_df[grepl("O$", rownames(NZC_B_df)),]))
NZC_B_specnumber_M<-data.frame(specnumber(NZC_B_df[grepl("M$", rownames(NZC_B_df)),]))

ILC_C_specnumber_B<-data.frame(specnumber(ILC_C_df[grepl("B$", rownames(ILC_C_df)),]))
ILC_C_specnumber_O<-data.frame(specnumber(ILC_C_df[grepl("O$", rownames(ILC_C_df)),]))
ILC_C_specnumber_M<-data.frame(specnumber(ILC_C_df[grepl("M$", rownames(ILC_C_df)),]))
NZC_C_specnumber_B<-data.frame(specnumber(NZC_C_df[grepl("B$", rownames(NZC_C_df)),]))
NZC_C_specnumber_O<-data.frame(specnumber(NZC_C_df[grepl("O$", rownames(NZC_C_df)),]))
NZC_C_specnumber_M<-data.frame(specnumber(NZC_C_df[grepl("M$", rownames(NZC_C_df)),]))

ILC_D_specnumber_B<-data.frame(specnumber(ILC_D_df[grepl("B$", rownames(ILC_D_df)),]))
ILC_D_specnumber_O<-data.frame(specnumber(ILC_D_df[grepl("O$", rownames(ILC_D_df)),]))
ILC_D_specnumber_M<-data.frame(specnumber(ILC_D_df[grepl("M$", rownames(ILC_D_df)),]))
NZC_D_specnumber_B<-data.frame(specnumber(NZC_D_df[grepl("B$", rownames(NZC_D_df)),]))
NZC_D_specnumber_O<-data.frame(specnumber(NZC_D_df[grepl("O$", rownames(NZC_D_df)),]))
NZC_D_specnumber_M<-data.frame(specnumber(NZC_D_df[grepl("M$", rownames(NZC_D_df)),]))

ILC_E_specnumber_B<-data.frame(specnumber(ILC_E_df[grepl("B$", rownames(ILC_E_df)),]))
ILC_E_specnumber_O<-data.frame(specnumber(ILC_E_df[grepl("O$", rownames(ILC_E_df)),]))
ILC_E_specnumber_M<-data.frame(specnumber(ILC_E_df[grepl("M$", rownames(ILC_E_df)),]))
NZC_E_specnumber_B<-data.frame(specnumber(NZC_E_df[grepl("B$", rownames(NZC_E_df)),]))
NZC_E_specnumber_O<-data.frame(specnumber(NZC_E_df[grepl("O$", rownames(NZC_E_df)),]))
NZC_E_specnumber_M<-data.frame(specnumber(NZC_E_df[grepl("M$", rownames(NZC_E_df)),]))

ILC_F_specnumber_B<-data.frame(specnumber(ILC_F_df[grepl("B$", rownames(ILC_F_df)),]))
ILC_F_specnumber_O<-data.frame(specnumber(ILC_F_df[grepl("O$", rownames(ILC_F_df)),]))
ILC_F_specnumber_M<-data.frame(specnumber(ILC_F_df[grepl("M$", rownames(ILC_F_df)),]))
NZC_F_specnumber_B<-data.frame(specnumber(NZC_F_df[grepl("B$", rownames(NZC_F_df)),]))
NZC_F_specnumber_O<-data.frame(specnumber(NZC_F_df[grepl("O$", rownames(NZC_F_df)),]))
NZC_F_specnumber_M<-data.frame(specnumber(NZC_F_df[grepl("M$", rownames(NZC_F_df)),]))

ILC_G_specnumber_B<-data.frame(specnumber(ILC_G_df[grepl("B$", rownames(ILC_G_df)),]))
ILC_G_specnumber_O<-data.frame(specnumber(ILC_G_df[grepl("O$", rownames(ILC_G_df)),]))
ILC_G_specnumber_M<-data.frame(specnumber(ILC_G_df[grepl("M$", rownames(ILC_G_df)),]))
NZC_G_specnumber_B<-data.frame(specnumber(NZC_G_df[grepl("B$", rownames(NZC_G_df)),]))
NZC_G_specnumber_O<-data.frame(specnumber(NZC_G_df[grepl("O$", rownames(NZC_G_df)),]))
NZC_G_specnumber_M<-data.frame(specnumber(NZC_G_df[grepl("M$", rownames(NZC_G_df)),]))

# Reformat lists for ggplot
ILC_B_specnumber_B[,2]<-"ILC"
ILC_B_specnumber_B[,3]<-"Bryophyte"
names(ILC_B_specnumber_B)<-c("richness","site","layer")

ILC_B_specnumber_O[,2]<-"ILC"
ILC_B_specnumber_O[,3]<-"Organic"
names(ILC_B_specnumber_O)<-c("richness","site","layer")

ILC_B_specnumber_M[,2]<-"ILC"
ILC_B_specnumber_M[,3]<-"Mineral"
names(ILC_B_specnumber_M)<-c("richness","site","layer")

NZC_B_specnumber_B[,2]<-"NZC"
NZC_B_specnumber_B[,3]<-"Bryophyte"
names(NZC_B_specnumber_B)<-c("richness","site","layer")

NZC_B_specnumber_O[,2]<-"NZC"
NZC_B_specnumber_O[,3]<-"Organic"
names(NZC_B_specnumber_O)<-c("richness","site","layer")

NZC_B_specnumber_M[,2]<-"NZC"
NZC_B_specnumber_M[,3]<-"Mineral"
names(NZC_B_specnumber_M)<-c("richness","site","layer")

summary1<-rbind(ILC_B_specnumber_B, ILC_B_specnumber_O, ILC_B_specnumber_M, 
                NZC_B_specnumber_B, NZC_B_specnumber_O, NZC_B_specnumber_M)
summary1$cores<-"1"
summary1$extractions<-"3"

ILC_C_specnumber_B[,2]<-"ILC"
ILC_C_specnumber_B[,3]<-"Bryophyte"
names(ILC_C_specnumber_B)<-c("richness","site","layer")

ILC_C_specnumber_O[,2]<-"ILC"
ILC_C_specnumber_O[,3]<-"Organic"
names(ILC_C_specnumber_O)<-c("richness","site","layer")

ILC_C_specnumber_M[,2]<-"ILC"
ILC_C_specnumber_M[,3]<-"Mineral"
names(ILC_C_specnumber_M)<-c("richness","site","layer")

NZC_C_specnumber_B[,2]<-"NZC"
NZC_C_specnumber_B[,3]<-"Bryophyte"
names(NZC_C_specnumber_B)<-c("richness","site","layer")

NZC_C_specnumber_O[,2]<-"NZC"
NZC_C_specnumber_O[,3]<-"Organic"
names(NZC_C_specnumber_O)<-c("richness","site","layer")

NZC_C_specnumber_M[,2]<-"NZC"
NZC_C_specnumber_M[,3]<-"Mineral"
names(NZC_C_specnumber_M)<-c("richness","site","layer")

summary2<-rbind(ILC_C_specnumber_B, ILC_C_specnumber_O, ILC_C_specnumber_M, 
                NZC_C_specnumber_B, NZC_C_specnumber_O, NZC_C_specnumber_M)
summary2$cores<-"2"
summary2$extractions<-"3"

ILC_D_specnumber_B[,2]<-"ILC"
ILC_D_specnumber_B[,3]<-"Bryophyte"
names(ILC_D_specnumber_B)<-c("richness","site","layer")

ILC_D_specnumber_O[,2]<-"ILC"
ILC_D_specnumber_O[,3]<-"Organic"
names(ILC_D_specnumber_O)<-c("richness","site","layer")

ILC_D_specnumber_M[,2]<-"ILC"
ILC_D_specnumber_M[,3]<-"Mineral"
names(ILC_D_specnumber_M)<-c("richness","site","layer")

NZC_D_specnumber_B[,2]<-"NZC"
NZC_D_specnumber_B[,3]<-"Bryophyte"
names(NZC_D_specnumber_B)<-c("richness","site","layer")

NZC_D_specnumber_O[,2]<-"NZC"
NZC_D_specnumber_O[,3]<-"Organic"
names(NZC_D_specnumber_O)<-c("richness","site","layer")

NZC_D_specnumber_M[,2]<-"NZC"
NZC_D_specnumber_M[,3]<-"Mineral"
names(NZC_D_specnumber_M)<-c("richness","site","layer")

summary3<-rbind(ILC_D_specnumber_B, ILC_D_specnumber_O, ILC_D_specnumber_M, 
                NZC_D_specnumber_B, NZC_D_specnumber_O, NZC_D_specnumber_M)
summary3$cores<-"4"
summary3$extractions<-"3"

ILC_E_specnumber_B[,2]<-"ILC"
ILC_E_specnumber_B[,3]<-"Bryophyte"
names(ILC_E_specnumber_B)<-c("richness","site","layer")

ILC_E_specnumber_O[,2]<-"ILC"
ILC_E_specnumber_O[,3]<-"Organic"
names(ILC_E_specnumber_O)<-c("richness","site","layer")

ILC_E_specnumber_M[,2]<-"ILC"
ILC_E_specnumber_M[,3]<-"Mineral"
names(ILC_E_specnumber_M)<-c("richness","site","layer")

NZC_E_specnumber_B[,2]<-"NZC"
NZC_E_specnumber_B[,3]<-"Bryophyte"
names(NZC_E_specnumber_B)<-c("richness","site","layer")

NZC_E_specnumber_O[,2]<-"NZC"
NZC_E_specnumber_O[,3]<-"Organic"
names(NZC_E_specnumber_O)<-c("richness","site","layer")

NZC_E_specnumber_M[,2]<-"NZC"
NZC_E_specnumber_M[,3]<-"Mineral"
names(NZC_E_specnumber_M)<-c("richness","site","layer")

summary4<-rbind(ILC_E_specnumber_B, ILC_E_specnumber_O, ILC_E_specnumber_M, 
                NZC_E_specnumber_B, NZC_E_specnumber_O, NZC_E_specnumber_M)
summary4$cores<-"6"
summary4$extractions<-"3"

ILC_F_specnumber_B[,2]<-"ILC"
ILC_F_specnumber_B[,3]<-"Bryophyte"
names(ILC_F_specnumber_B)<-c("richness","site","layer")

ILC_F_specnumber_O[,2]<-"ILC"
ILC_F_specnumber_O[,3]<-"Organic"
names(ILC_F_specnumber_O)<-c("richness","site","layer")

ILC_F_specnumber_M[,2]<-"ILC"
ILC_F_specnumber_M[,3]<-"Mineral"
names(ILC_F_specnumber_M)<-c("richness","site","layer")

NZC_F_specnumber_B[,2]<-"NZC"
NZC_F_specnumber_B[,3]<-"Bryophyte"
names(NZC_F_specnumber_B)<-c("richness","site","layer")

NZC_F_specnumber_O[,2]<-"NZC"
NZC_F_specnumber_O[,3]<-"Organic"
names(NZC_F_specnumber_O)<-c("richness","site","layer")

NZC_F_specnumber_M[,2]<-"NZC"
NZC_F_specnumber_M[,3]<-"Mineral"
names(NZC_F_specnumber_M)<-c("richness","site","layer")

summary5<-rbind(ILC_F_specnumber_B, ILC_F_specnumber_O, ILC_F_specnumber_M, 
                NZC_F_specnumber_B, NZC_F_specnumber_O, NZC_F_specnumber_M)
summary5$cores<-"8"
summary5$extractions<-"3"

ILC_G_specnumber_B[,2]<-"ILC"
ILC_G_specnumber_B[,3]<-"Bryophyte"
names(ILC_G_specnumber_B)<-c("richness","site","layer")

ILC_G_specnumber_O[,2]<-"ILC"
ILC_G_specnumber_O[,3]<-"Organic"
names(ILC_G_specnumber_O)<-c("richness","site","layer")

ILC_G_specnumber_M[,2]<-"ILC"
ILC_G_specnumber_M[,3]<-"Mineral"
names(ILC_G_specnumber_M)<-c("richness","site","layer")

NZC_G_specnumber_B[,2]<-"NZC"
NZC_G_specnumber_B[,3]<-"Bryophyte"
names(NZC_G_specnumber_B)<-c("richness","site","layer")

NZC_G_specnumber_O[,2]<-"NZC"
NZC_G_specnumber_O[,3]<-"Organic"
names(NZC_G_specnumber_O)<-c("richness","site","layer")

NZC_G_specnumber_M[,2]<-"NZC"
NZC_G_specnumber_M[,3]<-"Mineral"
names(NZC_G_specnumber_M)<-c("richness","site","layer")

summary6<-rbind(ILC_G_specnumber_B, ILC_G_specnumber_O, ILC_G_specnumber_M, 
                NZC_G_specnumber_B, NZC_G_specnumber_O, NZC_G_specnumber_M)
summary6$cores<-"9_15"
summary6$extractions<-"3"

summary_num<-rbind(summary1, summary2, summary3, summary4, summary5, summary6)

# Create factors
summary_num$layer<-factor(summary_num$layer, levels=c("Bryophyte", "Organic", "Mineral"))
summary_num$site<-factor(summary_num$site, levels=c("ILC", "NZC"),
                          labels=c("Island Lake","Nimitz"))

# create boxplot plot for 1C3E vs XC3E
p2<-ggplot(summary_num, aes(x=cores, y=richness, fill=layer)) +
  ggtitle("B)") +
  facet_wrap(~site+layer) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  labs(x="Manually pooled samples",y="ESV Richness") +
  scale_fill_manual(values=c("#4DAF4A","#377EB8","#FF7F00")) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(), 
        legend.position = "none")

###################################################################
##### Calculate species accumulation for 1C1E and 1C3E
# Extractions vs Richness
###################################################################

# do specnum for each site, layer, and extraction separately
ILC.1<-ILC_B2_df[grepl("^ILC45_", rownames(ILC_B2_df)),]
NZC.1<-NZC_B2_df[grepl("^NZC85_", rownames(NZC_B2_df)),]

ILC_1<-ILC.1[grepl("_1C1E_", rownames(ILC.1)),]
NZC_1<-NZC.1[grepl("_1C1E_", rownames(NZC.1)),]

ILC.3<-ILC_B3_df[grepl("^ILC45_", rownames(ILC_B3_df)),]
NZC.3<-NZC_B3_df[grepl("^NZC85_", rownames(NZC_B3_df)),]

ILC_3<-ILC.3[grepl("_1C3E_", rownames(ILC.3)),]
NZC_3<-NZC.3[grepl("_1C3E_", rownames(NZC.3)),]

ILC_1_B_specnumber<-data.frame(specnumber(ILC_1[grepl("B$", rownames(ILC_1)),]))
ILC_1_O_specnumber<-data.frame(specnumber(ILC_1[grepl("O$", rownames(ILC_1)),]))
ILC_1_M_specnumber<-data.frame(specnumber(ILC_1[grepl("M$", rownames(ILC_1)),]))

ILC_3_B_specnumber<-data.frame(specnumber(ILC_3[grepl("B$", rownames(ILC_3)),]))
ILC_3_O_specnumber<-data.frame(specnumber(ILC_3[grepl("O$", rownames(ILC_3)),]))
ILC_3_M_specnumber<-data.frame(specnumber(ILC_3[grepl("M$", rownames(ILC_3)),]))

NZC_1_B_specnumber<-data.frame(specnumber(NZC_1[grepl("B$", rownames(NZC_1)),]))
NZC_1_O_specnumber<-data.frame(specnumber(NZC_1[grepl("O$", rownames(NZC_1)),]))
NZC_1_M_specnumber<-data.frame(specnumber(NZC_1[grepl("M$", rownames(NZC_1)),]))

NZC_3_B_specnumber<-data.frame(specnumber(NZC_3[grepl("B$", rownames(NZC_3)),]))
NZC_3_O_specnumber<-data.frame(specnumber(NZC_3[grepl("O$", rownames(NZC_3)),]))
NZC_3_M_specnumber<-data.frame(specnumber(NZC_3[grepl("M$", rownames(NZC_3)),]))

# Reformat lists for ggplot
foo <- data.frame(do.call('rbind', strsplit(as.character(rownames(ILC_1_B_specnumber)),'_',fixed=TRUE)))
ILC_1_B_specnumber[,2]<-foo$X1
ILC_1_B_specnumber[,3]<-foo$X3
ILC_1_B_specnumber[,4]<-substring(foo$X4, 3, 3) 
ILC_1_B_specnumber[,5]<-substring(foo$X4, 0, 2)
names(ILC_1_B_specnumber)<-c("richness","site", "extractions","layers", "gridcoord")

foo <- data.frame(do.call('rbind', strsplit(as.character(rownames(ILC_1_O_specnumber)),'_',fixed=TRUE)))
ILC_1_O_specnumber[,2]<-foo$X1
ILC_1_O_specnumber[,3]<-foo$X3
ILC_1_O_specnumber[,4]<-substring(foo$X4, 3, 3) 
ILC_1_O_specnumber[,5]<-substring(foo$X4, 0, 2)
names(ILC_1_O_specnumber)<-c("richness","site","extractions", "layers","gridcoord")

foo <- data.frame(do.call('rbind', strsplit(as.character(rownames(ILC_1_M_specnumber)),'_',fixed=TRUE)))
ILC_1_M_specnumber[,2]<-foo$X1
ILC_1_M_specnumber[,3]<-foo$X3
ILC_1_M_specnumber[,4]<-substring(foo$X4, 3, 3)
ILC_1_M_specnumber[,5]<-substring(foo$X4, 0, 2)
names(ILC_1_M_specnumber)<-c("richness","site","extractions", "layers","gridcoord")

foo <- data.frame(do.call('rbind', strsplit(as.character(rownames(ILC_3_B_specnumber)),'_',fixed=TRUE)))
ILC_3_B_specnumber[,2]<-foo$X1
ILC_3_B_specnumber[,3]<-foo$X3
ILC_3_B_specnumber[,4]<-substring(foo$X4, 3, 3)
ILC_3_B_specnumber[,5]<-substring(foo$X4, 0, 2)
names(ILC_3_B_specnumber)<-c("richness","site","extractions", "layers","gridcoord")

foo <- data.frame(do.call('rbind', strsplit(as.character(rownames(ILC_3_O_specnumber)),'_',fixed=TRUE)))
ILC_3_O_specnumber[,2]<-foo$X1
ILC_3_O_specnumber[,3]<-foo$X3
ILC_3_O_specnumber[,4]<-substring(foo$X4, 3, 3)
ILC_3_O_specnumber[,5]<-substring(foo$X4, 0, 2)
names(ILC_3_O_specnumber)<-c("richness","site","extractions", "layers","gridcoord")

foo <- data.frame(do.call('rbind', strsplit(as.character(rownames(ILC_3_M_specnumber)),'_',fixed=TRUE)))
ILC_3_M_specnumber[,2]<-foo$X1
ILC_3_M_specnumber[,3]<-foo$X3
ILC_3_M_specnumber[,4]<-substring(foo$X4, 3, 3)
ILC_3_M_specnumber[,5]<-substring(foo$X4, 0, 2)
names(ILC_3_M_specnumber)<-c("richness","site","extractions", "layers","gridcoord")

foo <- data.frame(do.call('rbind', strsplit(as.character(rownames(NZC_1_B_specnumber)),'_',fixed=TRUE)))
NZC_1_B_specnumber[,2]<-foo$X1
NZC_1_B_specnumber[,3]<-foo$X3
NZC_1_B_specnumber[,4]<-substring(foo$X4, 3, 3)
NZC_1_B_specnumber[,5]<-substring(foo$X4, 0, 2)
names(NZC_1_B_specnumber)<-c("richness","site","extractions", "layers","gridcoord")

foo <- data.frame(do.call('rbind', strsplit(as.character(rownames(NZC_1_O_specnumber)),'_',fixed=TRUE)))
NZC_1_O_specnumber[,2]<-foo$X1
NZC_1_O_specnumber[,3]<-foo$X3
NZC_1_O_specnumber[,4]<-substring(foo$X4, 3, 3)
NZC_1_O_specnumber[,5]<-substring(foo$X4, 0, 2)
names(NZC_1_O_specnumber)<-c("richness","site","extractions", "layers","gridcoord")

foo <- data.frame(do.call('rbind', strsplit(as.character(rownames(NZC_1_M_specnumber)),'_',fixed=TRUE)))
NZC_1_M_specnumber[,2]<-foo$X1
NZC_1_M_specnumber[,3]<-foo$X3
NZC_1_M_specnumber[,4]<-substring(foo$X4, 3, 3)
NZC_1_M_specnumber[,5]<-substring(foo$X4, 0, 2)
names(NZC_1_M_specnumber)<-c("richness","site","extractions", "layers","gridcoord")

foo <- data.frame(do.call('rbind', strsplit(as.character(rownames(NZC_3_B_specnumber)),'_',fixed=TRUE)))
NZC_3_B_specnumber[,2]<-foo$X1
NZC_3_B_specnumber[,3]<-foo$X3
NZC_3_B_specnumber[,4]<-substring(foo$X4, 3, 3)
NZC_3_B_specnumber[,5]<-substring(foo$X4, 0, 2)
names(NZC_3_B_specnumber)<-c("richness","site","extractions", "layers","gridcoord")

foo <- data.frame(do.call('rbind', strsplit(as.character(rownames(NZC_3_O_specnumber)),'_',fixed=TRUE)))
NZC_3_O_specnumber[,2]<-foo$X1
NZC_3_O_specnumber[,3]<-foo$X3
NZC_3_O_specnumber[,4]<-substring(foo$X4, 3, 3)
NZC_3_O_specnumber[,5]<-substring(foo$X4, 0, 2)
names(NZC_3_O_specnumber)<-c("richness","site","extractions", "layers","gridcoord")

foo <- data.frame(do.call('rbind', strsplit(as.character(rownames(NZC_3_M_specnumber)),'_',fixed=TRUE)))
NZC_3_M_specnumber[,2]<-foo$X1
NZC_3_M_specnumber[,3]<-foo$X3
NZC_3_M_specnumber[,4]<-substring(foo$X4, 3, 3)
NZC_3_M_specnumber[,5]<-substring(foo$X4, 0, 2)
names(NZC_3_M_specnumber)<-c("richness","site","extractions", "layers","gridcoord")

summary_ext<-rbind(ILC_1_B_specnumber, ILC_1_O_specnumber, ILC_1_M_specnumber, 
                   ILC_3_B_specnumber, ILC_3_O_specnumber, ILC_3_M_specnumber,
                   NZC_1_B_specnumber, NZC_1_O_specnumber, NZC_1_M_specnumber, 
                   NZC_3_B_specnumber, NZC_3_O_specnumber, NZC_3_M_specnumber)

# For each experiment, calc richness diff
extraction_df<-dcast(summary_ext, site + layers + gridcoord ~ extractions, value.var="richness", fun.aggregate = sum)
extraction_df$diff<-(extraction_df[,5] - extraction_df[,4])

# Create factors
extraction_df$layers<-factor(extraction_df$layers, levels=c("B", "O", "M"), labels=c("Bryophyte","Organic","Mineral"))
extraction_df$site<-factor(extraction_df$site, levels=c("ILC45", "NZC85"),
                           labels=c("Island Lake","Nimitz"))

# create boxplot plot for 1C3E vs XC3E
p3<-ggplot(extraction_df, aes(x=layers, y=diff, fill=layers)) +
  ggtitle("C)") +
  facet_wrap(~site) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  labs(x="Layer",y="ESV Richness\n3 extractions - 1 extraction") +
  scale_fill_manual(values=c("#4DAF4A","#377EB8","#FF7F00")) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(), 
        legend.position = "none")

lay<-rbind(c(1),
           c(1),
           c(2),
           c(2),
           c(3))
  
g<-grid.arrange(p1, p2, p3, layout_matrix=lay)

ggsave("F2_richness_ESV.pdf", g, width=8, height=10, units="in")

#####################################
# Check for normality with Shapiro-Wilk’s test, sig result means not normal
shapiro.test(summary_ext$richness)
# W = 0.94999, p-value = 0.001259

#visual inspection, sometimes small sample sizes can pass normality tests
qqPlot(summary_ext$richness)

# Use Kruskal-Wallis test to check for any significant differences among extractions
new<-summary_ext
new$extractions<-as.factor(new$extractions)
kruskal.test(richness ~ extractions, data = new)
# Kruskal-Wallis chi-squared = 0.41817, df = 1, p-value = 0.5179

# Use multiple pairwise-comparison bewteen groups to check for specific diffs across cores just in case 
# p.adjust method Benjamini & Hochberg (1995)
pairwise.wilcox.test(new$richness, new$extractions,
                     p.adjust.method = "BH")
# n/s

#################################################################
### Create venn diagrams from rarefied data

# calc TOTAL richness from rarefied ILC_Arth_df and NZC_Arth_df
ILC.norm.tot.richness = specnumber(colSums(ILC_Arth_df))
# richness of 2,108 ESVs across all rarefied ILC data
NZC.norm.tot.richness = specnumber(colSums(NZC_Arth_df))
# richness of 2,052 ESVs across all rarefied NZC data

# ILC normalized library size, pooled across all samples
# Bryophyte presence-absence counts
b <- data.frame(rowSums(t(ILC_Arth_df[grepl("B$", rownames(ILC_Arth_df)),])))
b[b>0] <-1
# Organic presence-absence counts
o <- data.frame(rowSums(t(ILC_Arth_df[grepl("O$", rownames(ILC_Arth_df)),])))
o[o>0] <-1
# Mineral presence-absence counts
m <- data.frame(rowSums(t(ILC_Arth_df[grepl("M$", rownames(ILC_Arth_df)),])))
m[m>0] <-1
# create conunts for Venn
ILC.venndata <- cbind(b,o,m)
names(ILC.venndata) <- c("Bryophyte","Organic","Mineral")
ILC.venncounts <- vennCounts(ILC.venndata)
# ignore ESVs with zero counts after rarefaction
ILC.venncounts[1,4] <- 0

# draw venn diagram with ggplot instead
# df with coordinates for circles
df.venn <- data.frame(x = c(-0.866, 0.866, 0),
                      y = c(1, 1, -0.5),
                      labels = c('Bryophyte', 'Organic', 'Mineral'))
df.venn$labels <- factor(df.venn$labels, levels=c('Bryophyte', 'Organic', 'Mineral'))

# transform venncounts into df, add coordinates for labels
class(ILC.venncounts) <- 'matrix'
df.ILC.venncounts <- as.data.frame(ILC.venncounts)[-1,] %>%
  mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),
         y = c(-0.8, 1.2, 0, 1.2, 0, 1.4, 0.4))

# draw venn with ggplot

ILCrichness <- paste("Island Lake ESVs =",
                     ILC.norm.tot.richness)

v1 <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  ggtitle(paste0(ILCrichness)) +
  geom_circle(aes(fill=labels), size = 1) +
  geom_circle(size = 1, show.legend = FALSE, color = "black", alpha=0.3) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#4DAF4A","#377EB8", "#FF7F00")) +
  scale_color_manual(values=c("#4DAF4A","#377EB8", "#FF7F00"), guide = FALSE) +
  labs(fill = NULL) +
  annotate("text", x = df.ILC.venncounts$x, y = df.ILC.venncounts$y, label = df.ILC.venncounts$Counts, size = 5)

# NZC normalized library size, pooled across all samples
# Bryophyte presence-absence counts
b <- data.frame(rowSums(t(NZC_Arth_df[grepl("B$", rownames(NZC_Arth_df)),])))
b[b>0] <-1
# Organic presence-absence counts
o <- data.frame(rowSums(t(NZC_Arth_df[grepl("O$", rownames(NZC_Arth_df)),])))
o[o>0] <-1
# Mineral presence-absence counts
m <- data.frame(rowSums(t(NZC_Arth_df[grepl("M$", rownames(NZC_Arth_df)),])))
m[m>0] <-1
# create conunts for Venn
NZC.venndata <- cbind(b,o,m)
names(NZC.venndata) <- c("Bryophyte","Organic","Mineral")
NZC.venncounts <- vennCounts(NZC.venndata)
# ignore OTUs with zero counts after rarefaction
NZC.venncounts[1,4] <- 0

# draw venn diagram with ggplot instead, get coord from above
df.venn

# transform venncounts into df, add coordinates for labels
class(NZC.venncounts) <- 'matrix'
df.NZC.venncounts <- as.data.frame(NZC.venncounts)[-1,] %>%
  mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),
         y = c(-0.8, 1.2, 0, 1.2, 0, 1.4, 0.4))

# draw venn with ggplot

NZCrichness <- paste("Nimitz ESVs =",
                     NZC.norm.tot.richness)

v2 <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  ggtitle(paste0(NZCrichness)) +
  geom_circle(aes(fill=labels), size = 1) +
  geom_circle(size = 1, show.legend = FALSE, color = "black", alpha=0.3) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#4DAF4A","#377EB8", "#FF7F00")) +
  scale_color_manual(values=c("#4DAF4A","#377EB8", "#FF7F00"), guide = FALSE) +
  labs(fill = NULL) +
  annotate("text", x = df.NZC.venncounts$x, y = df.NZC.venncounts$y, label = df.NZC.venncounts$Counts, size = 5)

lay <- rbind(c(1,2),
             c(1,2))

g <- grid.arrange(v1, v2, layout_matrix = lay)

ggsave("FigS4_Venn.pdf", g)

#################################################################
### Calc greatest richness detected including up to 15 individually collected cores versus up to 15 pooled cores

## Bioinformatically pooled
# convert to presence-absence ILC
ILC_B.combo.df_df[ILC_B.combo.df_df > 0] <- 1
# check 15 indvidually pooled cores
median(rowSums(ILC_B.combo.df_df[grepl("_15_", rownames(ILC_B.combo.df_df)),]))
# 595.5 
sd(rowSums(ILC_B.combo.df_df[grepl("_15_", rownames(ILC_B.combo.df_df)),]))
# 153.1

# convert to presence-absence NZC
NZC_B.combo.df_df[NZC_B.combo.df_df > 0] <- 1
# check 15 indvidually pooled cores
median(rowSums(NZC_B.combo.df_df[grepl("_15_", rownames(NZC_B.combo.df_df)),]))
# 510
sd(rowSums(NZC_B.combo.df_df[grepl("_15_", rownames(NZC_B.combo.df_df)),]))
# 128.1

## Manually pooled
# convert to presence-absence ILC
ILC_G_df[ILC_G_df > 0] <- 1
# check 15 indvidually pooled cores
median(rowSums(ILC_G_df[grepl("_15C3E_", rownames(ILC_G_df)),]))
# 167
sd(rowSums(ILC_G_df[grepl("_15C3E_", rownames(ILC_G_df)),]))
# 74.5

# convert to presence-absence NZC
NZC_G_df[NZC_G_df > 0] <- 1
# check 15 indvidually pooled cores
median(rowSums(NZC_G_df[grepl("_15C3E_", rownames(NZC_G_df)),]))
# 126
sd(rowSums(NZC_G_df[grepl("_15C3E_", rownames(NZC_G_df)),]))
# 84.6

#########
# Figure out how many single cores would need to be sampled to match results from largest class of pooled cores
# summary_num, pool data across sites and layers, get median richness
median(summary_num$richness[summary_num$cores=="9_15" & 
                          (summary_num$layer == "Bryophyte" |
                          summary_num$layer == "Organic")])
# 205.5 ESVs
sd(summary_num$richness[summary_num$cores=="9_15" & 
                              (summary_num$layer == "Bryophyte" |
                                 summary_num$layer == "Organic")])
# 34.0 ESVs

# summary_bio, pool data across sites and layer, get median richness
median(summary_bio$richness[summary_bio$cores=="1" & 
                         (summary_bio$layer == "Bryophyte" |
                          summary_bio$layer == "Organic")])
# 177 ESVs
sd(summary_bio$richness[summary_bio$cores=="1" & 
                        (summary_bio$layer == "Bryophyte" |
                         summary_bio$layer == "Organic")])
# 43.6

# Use multiple pairwise-comparison bewteen groups to check for specific diffs across cores just in case 
# p.adjust method Benjamini & Hochberg (1995)

# grab data for 15 field pooled cores
composite915 <- data.frame(summary_num$richness[summary_num$cores=="9_15" & 
                       (summary_num$layer == "Bryophyte" |
                          summary_num$layer == "Organic")])
composite915$type <- "composite915"
names(composite915)[1] <- "richness"

# grab richness for 15 individual cores
individual15 <- data.frame(summary_bio$richness[summary_bio$cores=="15" & 
                       (summary_bio$layer == "Bryophyte" |
                          summary_bio$layer == "Organic")])
individual15$type <- "individual15"
names(individual15)[1] <- "richness"

# put into one df for statistical comparison
comp <- rbind(composite915, individual15)

pairwise.wilcox.test(comp$richness, comp$type,
                     p.adjust.method = "BH")
# 1.5e-06 

median(composite915$richness)
# 205.5
median(individual15$richness)
# 598.5

