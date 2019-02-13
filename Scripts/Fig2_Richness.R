# Teresita M. Porter, Feb. 13, 2019

library(vegan)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(grid)
library(gridExtra)
library(reshape2)
library(dplyr)
library("car")

###################################################################

#read infile prepared by python script
master<-read.table(file="LV2016_2.csv", head=TRUE, sep=",")

# Select phylum Arthropoda only
Arth_df<-master[master$phylum=="Arthropoda",]

# ESV matrix for 1C1E experiment
A<-Arth_df[Arth_df$extractions=="1",]

# ESV matrix for 1C3E experiment
B<-Arth_df[(Arth_df$cores=="1") & (Arth_df$extractions=="3"),]

# ESV matrix for 1C3E experiment subsampled down to 4 reps to compare with XC3E for balanced design
B2<-Arth_df[(Arth_df$cores=="1") & (Arth_df$extractions=="3"),]

# Randomly subsample 4 unique gridcoord without replacement 
B2.gridcoord<-sample(unique(B2$gridcoord), 4, replace=FALSE)

# Retrieve gridcoord
B2<-B2[B2$gridcoord=="11" |
         B2$gridcoord=="16" |
         B2$gridcoord=="33" |
         B2$gridcoord=="56",]

# ESV matrix for 1C3E experiment subsampled down to 8 reps to compare with XC3E for balanced design
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

# ESV matrix for 2C3E experiment
C<-Arth_df[(Arth_df$cores=="2") & (Arth_df$extractions=="3"),]

# ESV matrix for 4C3E experiment
D<-Arth_df[(Arth_df$cores=="4") & (Arth_df$extractions=="3"),]

# ESV matrix for 6C3E experiment
E<-Arth_df[(Arth_df$cores=="6") & (Arth_df$extractions=="3"),]

# ESV matrix for 8C3E experiment
F<-Arth_df[(Arth_df$cores=="8") & (Arth_df$extractions=="3"),]

# ESV matrix for 15C3E experiment
G<-Arth_df[((Arth_df$cores=="9") |
            (Arth_df$cores=="12") | (Arth_df$cores=="13") | 
            (Arth_df$cores=="14") | (Arth_df$cores=="15"))
           & (Arth_df$extractions=="3"),]

#split out ILC and NZC datasets
ILC_Arth<-Arth_df[Arth_df$site=="ILC45",]
NZC_Arth<-Arth_df[Arth_df$site=="NZC85",]
ILC_A<-A[A$site=="ILC45",]
NZC_A<-A[A$site=="NZC85",]
ILC_B<-B[B$site=="ILC45",]
NZC_B<-B[B$site=="NZC85",]
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

# Pivot to make matrices for vegan, ESVs in rows, sites in columns, reads in cells
matrixArth_ILC<-dcast(ILC_Arth, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)
matrixA_ILC<-dcast(ILC_A, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)
matrixB_ILC<-dcast(ILC_B, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)
matrixB2_ILC<-dcast(ILC_B2, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)
matrixB3_ILC<-dcast(ILC_B3, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)
matrixC_ILC<-dcast(ILC_C, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)
matrixD_ILC<-dcast(ILC_D, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)
matrixE_ILC<-dcast(ILC_E, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)
matrixF_ILC<-dcast(ILC_F, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)
matrixG_ILC<-dcast(ILC_G, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)
matrixArth_NZC<-dcast(NZC_Arth, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)
matrixA_NZC<-dcast(NZC_A, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)
matrixB_NZC<-dcast(NZC_B, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)
matrixB2_NZC<-dcast(NZC_B2, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)
matrixB3_NZC<-dcast(NZC_B3, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)
matrixC_NZC<-dcast(NZC_C, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)
matrixD_NZC<-dcast(NZC_D, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)
matrixE_NZC<-dcast(NZC_E, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)
matrixF_NZC<-dcast(NZC_F, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)
matrixG_NZC<-dcast(NZC_G, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)

# Move marker_OTU to row names
rownames(matrixArth_ILC)<-matrixArth_ILC$marker_OTU
matrixArth_ILC<-matrixArth_ILC[,-1]
rownames(matrixA_ILC)<-matrixA_ILC$marker_OTU
matrixA_ILC<-matrixA_ILC[,-1]
rownames(matrixB_ILC)<-matrixB_ILC$marker_OTU
matrixB_ILC<-matrixB_ILC[,-1]
rownames(matrixB2_ILC)<-matrixB2_ILC$marker_OTU
matrixB2_ILC<-matrixB2_ILC[,-1]
rownames(matrixB3_ILC)<-matrixB3_ILC$marker_OTU
matrixB3_ILC<-matrixB3_ILC[,-1]
rownames(matrixC_ILC)<-matrixC_ILC$marker_OTU
matrixC_ILC<-matrixC_ILC[,-1]
rownames(matrixD_ILC)<-matrixD_ILC$marker_OTU
matrixD_ILC<-matrixD_ILC[,-1]
rownames(matrixE_ILC)<-matrixE_ILC$marker_OTU
matrixE_ILC<-matrixE_ILC[,-1]
rownames(matrixF_ILC)<-matrixF_ILC$marker_OTU
matrixF_ILC<-matrixF_ILC[,-1]
rownames(matrixG_ILC)<-matrixG_ILC$marker_OTU
matrixG_ILC<-matrixG_ILC[,-1]
rownames(matrixArth_ILC)<-matrixArth_ILC$marker_OTU
matrixArth_NZC<-matrixArth_NZC[,-1]
rownames(matrixA_ILC)<-matrixA_ILC$marker_OTU
matrixA_NZC<-matrixA_NZC[,-1]
rownames(matrixB_NZC)<-matrixB_NZC$marker_OTU
matrixB_NZC<-matrixB_NZC[,-1]
rownames(matrixB2_NZC)<-matrixB2_NZC$marker_OTU
matrixB2_NZC<-matrixB2_NZC[,-1]
rownames(matrixB3_NZC)<-matrixB3_NZC$marker_OTU
matrixB3_NZC<-matrixB3_NZC[,-1]
rownames(matrixC_NZC)<-matrixC_NZC$marker_OTU
matrixC_NZC<-matrixC_NZC[,-1]
rownames(matrixD_NZC)<-matrixD_NZC$marker_OTU
matrixD_NZC<-matrixD_NZC[,-1]
rownames(matrixE_NZC)<-matrixE_NZC$marker_OTU
matrixE_NZC<-matrixE_NZC[,-1]
rownames(matrixF_NZC)<-matrixF_NZC$marker_OTU
matrixF_NZC<-matrixF_NZC[,-1]
rownames(matrixG_NZC)<-matrixG_NZC$marker_OTU
matrixG_NZC<-matrixG_NZC[,-1]

# Transpose to get sites in rows, ESVs in columns, reads in cells
matrixArth_ILC2<-t(matrixArth_ILC)
matrixA_ILC2<-t(matrixA_ILC)
matrixB_ILC2<-t(matrixB_ILC)
matrixB2_ILC2<-t(matrixB2_ILC)
matrixB3_ILC2<-t(matrixB3_ILC)
matrixC_ILC2<-t(matrixC_ILC)
matrixD_ILC2<-t(matrixD_ILC)
matrixE_ILC2<-t(matrixE_ILC)
matrixF_ILC2<-t(matrixF_ILC)
matrixG_ILC2<-t(matrixG_ILC)
matrixArth_NZC2<-t(matrixArth_NZC)
matrixA_NZC2<-t(matrixA_NZC)
matrixB_NZC2<-t(matrixB_NZC)
matrixB2_NZC2<-t(matrixB2_NZC)
matrixB3_NZC2<-t(matrixB3_NZC)
matrixC_NZC2<-t(matrixC_NZC)
matrixD_NZC2<-t(matrixD_NZC)
matrixE_NZC2<-t(matrixE_NZC)
matrixF_NZC2<-t(matrixF_NZC)
matrixG_NZC2<-t(matrixG_NZC)

#remove columns with only zeros
ILC_Arth_notnull<-matrixArth_ILC2[,colSums(matrixArth_ILC2) !=0]
ILC_A_notnull<-matrixA_ILC2[,colSums(matrixA_ILC2) !=0]
ILC_B_notnull<-matrixB_ILC2[,colSums(matrixB_ILC2) !=0]
ILC_B2_notnull<-matrixB2_ILC2[,colSums(matrixB2_ILC2) !=0]
ILC_B3_notnull<-matrixB3_ILC2[,colSums(matrixB3_ILC2) !=0]
ILC_C_notnull<-matrixC_ILC2[,colSums(matrixC_ILC2) !=0]
ILC_D_notnull<-matrixD_ILC2[,colSums(matrixD_ILC2) !=0]
ILC_E_notnull<-matrixE_ILC2[,colSums(matrixE_ILC2) !=0]
ILC_F_notnull<-matrixF_ILC2[,colSums(matrixF_ILC2) !=0]
ILC_G_notnull<-matrixG_ILC2[,colSums(matrixG_ILC2) !=0]
NZC_Arth_notnull<-matrixArth_NZC2[,colSums(matrixArth_NZC2) !=0]
NZC_A_notnull<-matrixA_NZC2[,colSums(matrixA_NZC2) !=0]
NZC_B_notnull<-matrixB_NZC2[,colSums(matrixB_NZC2) !=0]
NZC_B2_notnull<-matrixB2_NZC2[,colSums(matrixB2_NZC2) !=0]
NZC_B3_notnull<-matrixB3_NZC2[,colSums(matrixB3_NZC2) !=0]
NZC_C_notnull<-matrixC_NZC2[,colSums(matrixC_NZC2) !=0]
NZC_D_notnull<-matrixD_NZC2[,colSums(matrixD_NZC2) !=0]
NZC_E_notnull<-matrixE_NZC2[,colSums(matrixE_NZC2) !=0]
NZC_F_notnull<-matrixF_NZC2[,colSums(matrixF_NZC2) !=0]
NZC_G_notnull<-matrixG_NZC2[,colSums(matrixG_NZC2) !=0]

#remove rows with only zeros
ILC_Arth_notnull2<-ILC_Arth_notnull[rowSums(ILC_Arth_notnull) !=0,]
ILC_A_notnull2<-ILC_A_notnull[rowSums(ILC_A_notnull) !=0,]
ILC_B_notnull2<-ILC_B_notnull[rowSums(ILC_B_notnull) !=0,]
ILC_B2_notnull2<-ILC_B2_notnull[rowSums(ILC_B2_notnull) !=0,]
ILC_B3_notnull2<-ILC_B3_notnull[rowSums(ILC_B3_notnull) !=0,]
ILC_C_notnull2<-ILC_C_notnull[rowSums(ILC_C_notnull) !=0,]
ILC_D_notnull2<-ILC_D_notnull[rowSums(ILC_D_notnull) !=0,]
ILC_E_notnull2<-ILC_E_notnull[rowSums(ILC_E_notnull) !=0,]
ILC_F_notnull2<-ILC_F_notnull[rowSums(ILC_F_notnull) !=0,]
ILC_G_notnull2<-ILC_G_notnull[rowSums(ILC_G_notnull) !=0,]
NZC_Arth_notnull2<-NZC_Arth_notnull[rowSums(NZC_Arth_notnull) !=0,]
NZC_A_notnull2<-NZC_A_notnull[rowSums(NZC_A_notnull) !=0,]
NZC_B_notnull2<-NZC_B_notnull[rowSums(NZC_B_notnull) !=0,]
NZC_B2_notnull2<-NZC_B2_notnull[rowSums(NZC_B2_notnull) !=0,]
NZC_B3_notnull2<-NZC_B3_notnull[rowSums(NZC_B3_notnull) !=0,]
NZC_C_notnull2<-NZC_C_notnull[rowSums(NZC_C_notnull) !=0,]
NZC_D_notnull2<-NZC_D_notnull[rowSums(NZC_D_notnull) !=0,]
NZC_E_notnull2<-NZC_E_notnull[rowSums(NZC_E_notnull) !=0,]
NZC_F_notnull2<-NZC_F_notnull[rowSums(NZC_F_notnull) !=0,]
NZC_G_notnull2<-NZC_G_notnull[rowSums(NZC_G_notnull) !=0,]

#calculate 15th percentile for rrarefy function
ILC_Arth_15percentile<-quantile(rowSums(ILC_Arth_notnull2), prob=0.15)
ILC_A_15percentile<-quantile(rowSums(ILC_A_notnull2), prob=0.15)
ILC_B_15percentile<-quantile(rowSums(ILC_B_notnull2), prob=0.15)
ILC_B2_15percentile<-quantile(rowSums(ILC_B2_notnull2), prob=0.15)
ILC_B3_15percentile<-quantile(rowSums(ILC_B3_notnull2), prob=0.15)
ILC_C_15percentile<-quantile(rowSums(ILC_C_notnull2), prob=0.15)
ILC_D_15percentile<-quantile(rowSums(ILC_D_notnull2), prob=0.15)
ILC_E_15percentile<-quantile(rowSums(ILC_E_notnull2), prob=0.15)
ILC_F_15percentile<-quantile(rowSums(ILC_F_notnull2), prob=0.15)
ILC_G_15percentile<-quantile(rowSums(ILC_G_notnull2), prob=0.15)
NZC_Arth_15percentile<-quantile(rowSums(NZC_Arth_notnull2), prob=0.15)
NZC_A_15percentile<-quantile(rowSums(NZC_A_notnull2), prob=0.15)
NZC_B_15percentile<-quantile(rowSums(NZC_B_notnull2), prob=0.15)
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

#Rarefy original ILC matrices down to 15th percentile library size to normalize read depth across samples
ILC_Arth_df<-rrarefy(ILC_Arth_notnull2,sample=ILC_Arth_15percentile)
ILC_A_df<-rrarefy(ILC_A_notnull2,sample=ILC_A_15percentile)
ILC_B_df<-rrarefy(ILC_B_notnull2,sample=ILC_B_15percentile)
ILC_B2_df<-rrarefy(ILC_B2_notnull2,sample=ILC_B2_15percentile)
ILC_B3_df<-rrarefy(ILC_B3_notnull2,sample=ILC_B3_15percentile)
ILC_C_df<-rrarefy(ILC_C_notnull2,sample=ILC_C_15percentile)
ILC_D_df<-rrarefy(ILC_D_notnull2,sample=ILC_D_15percentile)
ILC_E_df<-rrarefy(ILC_E_notnull2,sample=ILC_E_15percentile)
ILC_F_df<-rrarefy(ILC_F_notnull2,sample=ILC_F_15percentile)
ILC_G_df<-rrarefy(ILC_G_notnull2,sample=ILC_G_15percentile)
NZC_Arth_df<-rrarefy(NZC_Arth_notnull2,sample=NZC_Arth_15percentile)
NZC_A_df<-rrarefy(NZC_A_notnull2,sample=NZC_A_15percentile)
NZC_B_df<-rrarefy(NZC_B_notnull2,sample=NZC_B_15percentile)
NZC_B2_df<-rrarefy(NZC_B2_notnull2,sample=NZC_B2_15percentile)
NZC_B3_df<-rrarefy(NZC_B3_notnull2,sample=NZC_B3_15percentile)
NZC_C_df<-rrarefy(NZC_C_notnull2,sample=NZC_C_15percentile)
NZC_D_df<-rrarefy(NZC_D_notnull2,sample=NZC_D_15percentile)
NZC_E_df<-rrarefy(NZC_E_notnull2,sample=NZC_E_15percentile)
NZC_F_df<-rrarefy(NZC_F_notnull2,sample=NZC_F_15percentile)
NZC_G_df<-rrarefy(NZC_G_notnull2,sample=NZC_G_15percentile)

###################################################################
##### Calculate species accumulation for 1C3E
# Cores vs ESV richness
###################################################################

# do specaccum for each site and layer separately
ILC_B_specaccum_B<-specaccum(ILC_B_df[grepl("B$", rownames(ILC_B_df)),], "random")
ILC_B_specaccum_O<-specaccum(ILC_B_df[grepl("O$", rownames(ILC_B_df)),], "random")
ILC_B_specaccum_M<-specaccum(ILC_B_df[grepl("M$", rownames(ILC_B_df)),], "random")
NZC_B_specaccum_B<-specaccum(NZC_B_df[grepl("B$", rownames(NZC_B_df)),], "random")
NZC_B_specaccum_O<-specaccum(NZC_B_df[grepl("O$", rownames(NZC_B_df)),], "random")
NZC_B_specaccum_M<-specaccum(NZC_B_df[grepl("M$", rownames(NZC_B_df)),], "random")

# Reformat lists for ggplot
df<-data.frame(matrix(ncol=5,nrow=36))
df[,1] <- data.frame(matrix(unlist(ILC_B_specaccum_B$sites)))
df[,2] <- data.frame(matrix(unlist(ILC_B_specaccum_B$richness)))
df[,3] <- data.frame(matrix(unlist(ILC_B_specaccum_B$sd)))
df[,4]<-"ILC"
df[,5]<-"Bryophyte"
names(df)<-c("samples","richness","sd","site","layer")

df2<-data.frame(matrix(ncol=5,nrow=36))
df2[,1] <- data.frame(matrix(unlist(ILC_B_specaccum_O$sites)))
df2[,2] <- data.frame(matrix(unlist(ILC_B_specaccum_O$richness)))
df2[,3] <- data.frame(matrix(unlist(ILC_B_specaccum_O$sd)))
df2[,4]<-"ILC"
df2[,5]<-"Organic"
names(df2)<-c("samples","richness","sd","site","layer")

df3<-data.frame(matrix(ncol=5,nrow=36))
df3[,1] <- data.frame(matrix(unlist(ILC_B_specaccum_M$sites)))
df3[,2] <- data.frame(matrix(unlist(ILC_B_specaccum_M$richness)))
df3[,3] <- data.frame(matrix(unlist(ILC_B_specaccum_M$sd)))
df3[,4]<-"ILC"
df3[,5]<-"Mineral"
names(df3)<-c("samples","richness","sd","site","layer")

df4<-data.frame(matrix(ncol=5,nrow=35))
df4[,1] <- data.frame(matrix(unlist(NZC_B_specaccum_B$sites)))
df4[,2] <- data.frame(matrix(unlist(NZC_B_specaccum_B$richness)))
df4[,3] <- data.frame(matrix(unlist(NZC_B_specaccum_B$sd)))
df4[,4]<-"NZC"
df4[,5]<-"Bryophyte"
names(df4)<-c("samples","richness","sd","site","layer")

df5<-data.frame(matrix(ncol=5,nrow=32))
df5[,1] <- data.frame(matrix(unlist(NZC_B_specaccum_O$sites)))
df5[,2] <- data.frame(matrix(unlist(NZC_B_specaccum_O$richness)))
df5[,3] <- data.frame(matrix(unlist(NZC_B_specaccum_O$sd)))
df5[,4]<-"NZC"
df5[,5]<-"Organic"
names(df5)<-c("samples","richness","sd","site","layer")

df6<-data.frame(matrix(ncol=5,nrow=36))
df6[,1] <- data.frame(matrix(unlist(NZC_B_specaccum_M$sites)))
df6[,2] <- data.frame(matrix(unlist(NZC_B_specaccum_M$richness)))
df6[,3] <- data.frame(matrix(unlist(NZC_B_specaccum_M$sd)))
df6[,4]<-"NZC"
df6[,5]<-"Mineral"
names(df6)<-c("samples","richness","sd","site","layer")

summary_acc<-rbind(df, df2, df3, df4, df5, df6)
summary_acc$cores<-"1"
summary_acc$extractions<-"3"

# Create factors
summary_acc$layer<-factor(summary_acc$layer, levels=c("Bryophyte","Organic","Mineral"))
summary_acc$site<-factor(summary_acc$site, levels=c("ILC","NZC"),
                         labels=c("Island Lake","Nimitz"))

# create specaccum plot for 1C3E
p1<-ggplot(summary_acc, aes(x=samples, y=richness, color=layer)) +
  ggtitle("A)") +
  geom_line() +
  geom_ribbon(aes(ymin=richness-sd, ymax=richness+sd, fill=layer, color=NULL), alpha=0.3) +
  facet_wrap(~site) +
  labs(x="Samples", y="ESV Richness") +
  scale_color_manual(values = c("#4DAF4A","#377EB8","#FF7F00")) +
  scale_fill_manual(values = c("#4DAF4A","#377EB8","#FF7F00")) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "right",
        legend.title = element_blank())

###################################################################
##### Calculate species accumulation for 1C3E & XC3E
# pooled cores vs ESV richness
###################################################################

# do specnum for each site and layer separately
ILC_B2_specnumber_B<-data.frame(specnumber(ILC_B2_df[grepl("B$", rownames(ILC_B2_df)),]))
ILC_B2_specnumber_O<-data.frame(specnumber(ILC_B2_df[grepl("O$", rownames(ILC_B2_df)),]))
ILC_B2_specnumber_M<-data.frame(specnumber(ILC_B2_df[grepl("M$", rownames(ILC_B2_df)),]))
NZC_B2_specnumber_B<-data.frame(specnumber(NZC_B2_df[grepl("B$", rownames(NZC_B2_df)),]))
NZC_B2_specnumber_O<-data.frame(specnumber(NZC_B2_df[grepl("O$", rownames(NZC_B2_df)),]))
NZC_B2_specnumber_M<-data.frame(specnumber(NZC_B2_df[grepl("M$", rownames(NZC_B2_df)),]))

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
ILC_B2_specnumber_B[,2]<-"ILC"
ILC_B2_specnumber_B[,3]<-"Bryophyte"
names(ILC_B2_specnumber_B)<-c("richness","site","layer")

ILC_B2_specnumber_O[,2]<-"ILC"
ILC_B2_specnumber_O[,3]<-"Organic"
names(ILC_B2_specnumber_O)<-c("richness","site","layer")

ILC_B2_specnumber_M[,2]<-"ILC"
ILC_B2_specnumber_M[,3]<-"Mineral"
names(ILC_B2_specnumber_M)<-c("richness","site","layer")

NZC_B2_specnumber_B[,2]<-"NZC"
NZC_B2_specnumber_B[,3]<-"Bryophyte"
names(NZC_B2_specnumber_B)<-c("richness","site","layer")

NZC_B2_specnumber_O[,2]<-"NZC"
NZC_B2_specnumber_O[,3]<-"Organic"
names(NZC_B2_specnumber_O)<-c("richness","site","layer")

NZC_B2_specnumber_M[,2]<-"NZC"
NZC_B2_specnumber_M[,3]<-"Mineral"
names(NZC_B2_specnumber_M)<-c("richness","site","layer")

summary1<-rbind(ILC_B2_specnumber_B, ILC_B2_specnumber_O, ILC_B2_specnumber_M, 
                NZC_B2_specnumber_B, NZC_B2_specnumber_O, NZC_B2_specnumber_M)
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
p3<-ggplot(summary_num, aes(x=cores, y=richness, fill=layer)) +
  ggtitle("B)") +
  facet_wrap(~site+layer) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  labs(x="Pooled samples",y="ESV Richness") +
  scale_fill_manual(values=c("#4DAF4A","#377EB8","#FF7F00")) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(), 
        legend.position = "none")

###################################################################
##### Calculate species accumulation for 1C1E
# Extractions vs Richness
###################################################################

# do specnum for each site and layer separately
ILC_A_specnumber_B<-data.frame(specnumber(ILC_A_df[grepl("B$", rownames(ILC_A_df)),]))
ILC_A_specnumber_O<-data.frame(specnumber(ILC_A_df[grepl("O$", rownames(ILC_A_df)),]))
ILC_A_specnumber_M<-data.frame(specnumber(ILC_A_df[grepl("M$", rownames(ILC_A_df)),]))
NZC_A_specnumber_B<-data.frame(specnumber(NZC_A_df[grepl("B$", rownames(NZC_A_df)),]))
NZC_A_specnumber_O<-data.frame(specnumber(NZC_A_df[grepl("O$", rownames(NZC_A_df)),]))
NZC_A_specnumber_M<-data.frame(specnumber(NZC_A_df[grepl("M$", rownames(NZC_A_df)),]))

ILC_A_specnumber_B[,2]<-"ILC"
ILC_A_specnumber_B[,3]<-"Bryophyte"
names(ILC_A_specnumber_B)<-c("richness","site","layer")

ILC_A_specnumber_O[,2]<-"ILC"
ILC_A_specnumber_O[,3]<-"Organic"
names(ILC_A_specnumber_O)<-c("richness","site","layer")

ILC_A_specnumber_M[,2]<-"ILC"
ILC_A_specnumber_M[,3]<-"Mineral"
names(ILC_A_specnumber_M)<-c("richness","site","layer")

NZC_A_specnumber_B[,2]<-"NZC"
NZC_A_specnumber_B[,3]<-"Bryophyte"
names(NZC_A_specnumber_B)<-c("richness","site","layer")

NZC_A_specnumber_O[,2]<-"NZC"
NZC_A_specnumber_O[,3]<-"Organic"
names(NZC_A_specnumber_O)<-c("richness","site","layer")

NZC_A_specnumber_M[,2]<-"NZC"
NZC_A_specnumber_M[,3]<-"Mineral"
names(NZC_A_specnumber_M)<-c("richness","site","layer")

summary7<-rbind(ILC_A_specnumber_B, ILC_A_specnumber_O, ILC_A_specnumber_M, 
                NZC_A_specnumber_B, NZC_A_specnumber_O, NZC_A_specnumber_M)
summary7$cores<-"1"
summary7$extractions<-"1"

ILC_B3_specnumber_B<-data.frame(specnumber(ILC_B3_df[grepl("B$", rownames(ILC_B3_df)),]))
ILC_B3_specnumber_O<-data.frame(specnumber(ILC_B3_df[grepl("O$", rownames(ILC_B3_df)),]))
ILC_B3_specnumber_M<-data.frame(specnumber(ILC_B3_df[grepl("M$", rownames(ILC_B3_df)),]))
NZC_B3_specnumber_B<-data.frame(specnumber(NZC_B3_df[grepl("B$", rownames(NZC_B3_df)),]))
NZC_B3_specnumber_O<-data.frame(specnumber(NZC_B3_df[grepl("O$", rownames(NZC_B3_df)),]))
NZC_B3_specnumber_M<-data.frame(specnumber(NZC_B3_df[grepl("M$", rownames(NZC_B3_df)),]))

ILC_B3_specnumber_B[,2]<-"ILC"
ILC_B3_specnumber_B[,3]<-"Bryophyte"
names(ILC_B3_specnumber_B)<-c("richness","site","layer")

ILC_B3_specnumber_O[,2]<-"ILC"
ILC_B3_specnumber_O[,3]<-"Organic"
names(ILC_B3_specnumber_O)<-c("richness","site","layer")

ILC_B3_specnumber_M[,2]<-"ILC"
ILC_B3_specnumber_M[,3]<-"Mineral"
names(ILC_B3_specnumber_M)<-c("richness","site","layer")

NZC_B3_specnumber_B[,2]<-"NZC"
NZC_B3_specnumber_B[,3]<-"Bryophyte"
names(NZC_B3_specnumber_B)<-c("richness","site","layer")

NZC_B3_specnumber_O[,2]<-"NZC"
NZC_B3_specnumber_O[,3]<-"Organic"
names(NZC_B3_specnumber_O)<-c("richness","site","layer")

NZC_B3_specnumber_M[,2]<-"NZC"
NZC_B3_specnumber_M[,3]<-"Mineral"
names(NZC_B3_specnumber_M)<-c("richness","site","layer")

summary8<-rbind(ILC_B3_specnumber_B, ILC_B3_specnumber_O, ILC_B3_specnumber_M, 
                NZC_B3_specnumber_B, NZC_B3_specnumber_O, NZC_B3_specnumber_M)
summary8$cores<-"1"
summary8$extractions<-"3"

summary_ext<-rbind(summary8, summary7)

# Create factors
summary_ext$layer<-factor(summary_ext$layer, levels=c("Bryophyte","Organic","Mineral"))
summary_ext$site<-factor(summary_ext$site, levels=c("ILC","NZC"),
                         labels=c("Island Lake","Nimitz"))

# Create boxplots to compare extractions by site &n layer
p4<-ggplot(summary_ext, aes(x=extractions,y=richness, fill=layer)) +
  ggtitle("C)") +
  geom_boxplot() +
  geom_jitter() +
  labs(x="Extractions",y="ESV Richness") +
  facet_wrap(~site+layer) +
  scale_fill_manual(values=c("#4DAF4A","#377EB8","#FF7F00")) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(), 
        legend.position = "none")

lay<-rbind(c(1),
           c(2),
           c(2),
           c(3),
           c(3))
  
g<-grid.arrange(p1, p3, p4, layout_matrix=lay)

ggsave("F3_richness.pdf", g, width=8, height=10, units="in")

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
