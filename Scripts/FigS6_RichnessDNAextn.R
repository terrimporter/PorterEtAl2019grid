# Teresita M. Porter, July 8, 2019

library(vegan)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(grid)
library(gridExtra)
library(reshape2)
library(dplyr)
library("car")
library(data.table) #rownames to col

########################################################
##### NEW FUNCTION TO GET TARGET GRIDCOORD DF's
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

#read infile
master<-read.table(file="LV2016_2.csv", head=TRUE, sep=",")

# Select phylum Arthropoda only
Arth_df<-master[master$phylum=="Arthropoda",]

# ESV matrix for 1C1E experiment
A<-Arth_df[Arth_df$extractions=="1",]

# ESV matrix for 1C3E experiment subsampled down to 8 reps to compare with 1C1E for balanced design
B3<-Arth_df[(Arth_df$cores=="1") & (Arth_df$extractions=="3"),]
# match the 8 gridcoord already used in 1C1E
# Retrieve gridcoord
B3<-B3[B3$gridcoord=="11" |
         B3$gridcoord=="15" |
         B3$gridcoord=="23" |
         B3$gridcoord=="31" |
         B3$gridcoord=="35" |
         B3$gridcoord=="43" |
         B3$gridcoord=="51" |
         B3$gridcoord=="55",]

#split out ILC and NZC datasets
ILC_A<-A[A$site=="ILC45",]
NZC_A<-A[A$site=="NZC85",]
ILC_B3<-B3[B3$site=="ILC45",]
NZC_B3<-B3[B3$site=="NZC85",]

# Pivot to make matrices for vegan, ESVs in rows, sites in columns, reads in cells
matrixA_ILC<-dcast(ILC_A, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)
matrixB3_ILC<-dcast(ILC_B3, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)
matrixA_NZC<-dcast(NZC_A, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)
matrixB3_NZC<-dcast(NZC_B3, marker_OTU ~ GRDIname, value.var="reads", fun.aggregate = sum)

# Move marker_OTU to row names
rownames(matrixA_ILC)<-matrixA_ILC$marker_OTU
matrixA_ILC<-matrixA_ILC[,-1]
rownames(matrixB3_ILC)<-matrixB3_ILC$marker_OTU
matrixB3_ILC<-matrixB3_ILC[,-1]
rownames(matrixA_ILC)<-matrixA_ILC$marker_OTU
matrixA_NZC<-matrixA_NZC[,-1]
rownames(matrixB3_NZC)<-matrixB3_NZC$marker_OTU
matrixB3_NZC<-matrixB3_NZC[,-1]

# Transpose to get sites in rows, ESVs in columns, reads in cells
matrixA_ILC2<-t(matrixA_ILC)
matrixB3_ILC2<-t(matrixB3_ILC)
matrixA_NZC2<-t(matrixA_NZC)
matrixB3_NZC2<-t(matrixB3_NZC)

#remove columns with only zeros
ILC_A_notnull<-matrixA_ILC2[,colSums(matrixA_ILC2) !=0]
ILC_B3_notnull<-matrixB3_ILC2[,colSums(matrixB3_ILC2) !=0]
NZC_A_notnull<-matrixA_NZC2[,colSums(matrixA_NZC2) !=0]
NZC_B3_notnull<-matrixB3_NZC2[,colSums(matrixB3_NZC2) !=0]

#remove rows with only zeros
ILC_A_notnull2<-ILC_A_notnull[rowSums(ILC_A_notnull) !=0,]
ILC_B3_notnull2<-ILC_B3_notnull[rowSums(ILC_B3_notnull) !=0,]
NZC_A_notnull2<-NZC_A_notnull[rowSums(NZC_A_notnull) !=0,]
NZC_B3_notnull2<-NZC_B3_notnull[rowSums(NZC_B3_notnull) !=0,]

#calculate 15th percentile for rrarefy function
ILC_A_15percentile<-quantile(rowSums(ILC_A_notnull2), prob=0.15)
ILC_B3_15percentile<-quantile(rowSums(ILC_B3_notnull2), prob=0.15)
NZC_A_15percentile<-quantile(rowSums(NZC_A_notnull2), prob=0.15)
NZC_B3_15percentile<-quantile(rowSums(NZC_B3_notnull2), prob=0.15)

set.seed(1234)

###################################################################
##### Rarefy the dataset down to the 15th percentile
###################################################################

ILC_A_df<-rrarefy(ILC_A_notnull2,sample=ILC_A_15percentile)
ILC_B3_df<-rrarefy(ILC_B3_notnull2,sample=ILC_B3_15percentile)
NZC_A_df<-rrarefy(NZC_A_notnull2,sample=NZC_A_15percentile)
NZC_B3_df<-rrarefy(NZC_B3_notnull2,sample=NZC_B3_15percentile)

###################################################################
##### Calculate richness for 1C1E
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

# Add column for grid coords so that difference between 1 and 3 can be compared for each gridcoord
foo <- data.frame(do.call('rbind', strsplit(as.character(rownames(summary_ext)),'_',fixed=TRUE)))
summary_ext$gridcoord<-foo$X4

# Create boxplots to compare extractions by site & layer
p1<-ggplot(summary_ext, aes(x=extractions,y=richness, fill=layer)) +
  geom_boxplot() +
  geom_jitter() +
  labs(x="DNA extractions",y="ESV Richness") +
  facet_wrap(~site+layer) +
  scale_fill_manual(values=c("#4DAF4A","#377EB8","#FF7F00")) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(), 
        legend.position = "none")

ggsave("FigS6_DNAextnRichness.pdf", p1, width=8, height=5, units="in")

# Check for normality with Shapiro-Wilk’s test, sig result means not normal
shapiro.test(summary_ext$richness)
# W = 0.9507, p-value = 0.001396

#visual inspection, sometimes small sample sizes can pass normality tests
qqPlot(summary_ext$richness)

# Use multiple pairwise-comparison bewteen groups to check for specific richness differences across cores just in case 
# p.adjust method Benjamini & Hochberg (1995)
pairwise.wilcox.test(summary_ext$richness, summary_ext$site, p.adjust.method = "BH")
# n/s
pairwise.wilcox.test(summary_ext$richness, summary_ext$layer, p.adjust.method = "BH")
# B-O n/s, B-M sig, O-M sig
pairwise.wilcox.test(summary_ext$richness, as.factor(summary_ext$extractions), p.adjust.method = "BH")
# n/s
