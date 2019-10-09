# Teresita M. Porter, Oct. 9, 2019

library(vegan)
library(reshape2)
library(gridExtra)
library(grid)
library(ggplot2)
library(plyr)
library(scales)

###################################################################

# Read infiles
SSU<-read.csv(file="18Sv4_rdp.csv", head=TRUE)
ITS.tmp<-read.csv(file="ITS2_rdp.csv", head=TRUE)
ITS<-ITS.tmp[grepl("^ITS2", ITS.tmp$ITS2_GlobalESV),]
BR5.tmp<-read.csv(file="BR5_rdp.csv", head=TRUE)
BR5.tmp2<-BR5.tmp[!grepl("Cannot", BR5.tmp$BR5_GlobalESV),]
F230R.tmp<-read.csv(file="F230R_rdp.csv", head=TRUE)
F230R.tmp2<-F230R.tmp[!grepl("Cannot", F230R.tmp$F230R_GlobalESV),]

# Select phylum Arthropoda only
BR5<-BR5.tmp2[BR5.tmp2$Phylum=="Arthropoda",]
F230R<-F230R.tmp2[F230R.tmp2$Phylum=="Arthropoda",]

# Get all ESV counts
length(unique(SSU$X18Sv4_GlobalESV))
#18,480
length(unique(ITS$ITS2_GlobalESV))
#16,596
length(unique(BR5$BR5_GlobalESV))
#3,968
length(unique(F230R$F230R_GlobalESV))
#9,521

# Get all ESV read counts
sum(SSU$ESVsize)
#13,106,296
sum(ITS$ESVsize)
#20,427,842
sum(BR5$ESVsize)
#371,622
sum(F230R$ESVsize)
#7,919,626

# Summarize ESVs in all detected phyla
SSU.2 <- dcast(SSU, Phylum ~ . , value.var="X18Sv4_GlobalESV", fun.aggregate=length)
names(SSU.2)<-c("Phyla","X18Sv4_GlobalESV")
ITS.2 <- dcast(ITS, Phylum ~ . , value.var="ITS2_GlobalESV", fun.aggregate=length)
names(ITS.2)<-c("Phyla","ITS2_GlobalESV")
BR5.2 <- dcast(BR5, Phylum ~ . , value.var="BR5_GlobalESV", fun.aggregate=length)
names(BR5.2)<-c("Phyla","BR5_GlobalESV")
F230R.2 <- dcast(F230R, Phylum ~ . , value.var="F230R_GlobalESV", fun.aggregate=length)
names(F230R.2)<-c("Phyla","F230R_GlobalESV")

# Sort descending
SSU.2<-SSU.2[order(-SSU.2$X18Sv4_GlobalESV),]
ITS.2<-ITS.2[order(-ITS.2$ITS2_GlobalESV),]
BR5.2<-BR5.2[order(-BR5.2$BR5_GlobalESV),]
F230R.2<-F230R.2[order(-F230R.2$F230R_GlobalESV),]

# Keep top 10
SSU.2<-SSU.2[1:10,]
ITS.2<-ITS.2[1:10,]
BR5.2<-BR5.2[1:10,]
F230R.2<-F230R.2[1:10,]

# Summarize reads in all detected phyla
SSU.3 <- dcast(SSU, Phylum ~ . , value.var="ESVsize", fun.aggregate=sum)
names(SSU.3)<-c("Phyla","ESVsize")
ITS.3 <- dcast(ITS, Phylum ~ . , value.var="ESVsize", fun.aggregate=sum)
names(ITS.3)<-c("Phyla","ESVsize")
BR5.3 <- dcast(BR5, Phylum ~ . , value.var="ESVsize", fun.aggregate=sum)
names(BR5.3)<-c("Phyla","ESVsize")
F230R.3 <- dcast(F230R, Phylum ~ . , value.var="ESVsize", fun.aggregate=sum)
names(F230R.3)<-c("Phyla","ESVsize")

# Sort descending
SSU.3<-SSU.3[order(-SSU.3$ESVsize),]
ITS.3<-ITS.3[order(-ITS.3$ESVsize),]
BR5.3<-BR5.3[order(-BR5.3$ESVsize),]
F230R.3<-F230R.3[order(-F230R.3$ESVsize),]

# Keep top 10
SSU.3<-SSU.3[1:10,]
ITS.3<-ITS.3[1:10,]
BR5.3<-BR5.3[1:10,]
F230R.3<-F230R.3[1:10,]

# Create phylum summary table
SSU.phylumTable<-merge(SSU.2, SSU.3, by="Phyla")
ITS.phylumTable<-merge(ITS.2, ITS.3, by="Phyla")
BR5.phylumTable<-merge(BR5.2, BR5.3, by="Phyla")
F230R.phylumTable<-merge(F230R.2, F230R.3, by="Phyla")

# Create long form form ggplot
SSU.phylumTable.long<-melt(SSU.phylumTable, id=c("Phyla"))
ITS.phylumTable.long<-melt(ITS.phylumTable, id=c("Phyla"))
BR5.phylumTable.long<-melt(BR5.phylumTable, id=c("Phyla"))
F230R.phylumTable.long<-melt(F230R.phylumTable, id=c("Phyla"))

# Create 100% stacked horizontal bar plot to summarize phyla diversity
p1<-ggplot(data=SSU.phylumTable.long, aes(x=variable, y=value, fill=Phyla)) +
  ggtitle("18Sv4") +
  geom_bar(position="fill", stat="identity") +
  labs(y="Proportion", x="Top 10 phyla") +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),
        axis.title.x=element_blank()) +
  # annotate top 3 ESVs
  annotate("text",x=1, y=0.73,label="33,481") +
  annotate("text",x=1, y=0.55,label="31,513") +
  annotate("text",x=1, y=0.12,label="30,606") +
  # annotate top 3 Reads
  annotate("text",x=2, y=0.12,label="3,561,942") +
  annotate("text",x=2, y=0.73,label="1,962,238") +
  annotate("text",x=2, y=0.61,label="1,848,910") 

ggsave("FigS1_phylumdist.pdf")

# Select phylum Arthropoda only
Arth_df<-A[A$phylum=="Arthropoda",]

# Get Arthropoda ESVs counts
length(unique(Arth_df$marker_OTU))
#3,598
length(unique(Arth_df$marker_OTU[Arth_df$marker=="BE"]))
#775
length(unique(Arth_df$marker_OTU[Arth_df$marker=="F230"]))
#2,823

# Count number of reads in arthropod ESVs
sum(Arth_df$reads)
#2,692,708
sum(Arth_df$reads[Arth_df$marker=="BE"])
#294,070
sum(Arth_df$reads[Arth_df$marker=="F230"])
#2,398,638
