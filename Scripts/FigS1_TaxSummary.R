# Teresita M. Porter, Feb. 13, 2019

library(vegan)
library(reshape2)
library(gridExtra)
library(grid)
library(ggplot2)
library(plyr)
library(scales)

###################################################################

# Read infile prepared by python script
## LV2016_2.csv is the taxonomic assignment table
A<-read.table(file="LV2016_2.csv", head=TRUE, sep=",")

# Get all ESV counts for Table S2
length(unique(A$marker_OTU))
#67,626
length(unique(A$marker_OTU[A$marker=="BE"]))
#50,857
length(unique(A$marker_OTU[A$marker=="F230"]))
#16,769

# Get all ESV read counts for Table S2
sum(A$reads)
#19,562,246
sum(A$reads[A$marker=="BE"])
#8,269,741
sum(A$reads[A$marker=="F230"])
#11,292,505

# Summarize ESVs in all detected phyla
A2 <- dcast(A, phylum ~ . , value.var="marker_OTU", fun.aggregate=length)
names(A2)<-c("Phyla","ESVs")

# Sort descending
A2<-A2[order(-A2$ESVs),]

# Keep top 10
A2<-A2[1:10,]

# Summarize reads in all detected phyla
A3 <- dcast(A, phylum ~ . , value.var="reads", fun.aggregate=sum)
names(A3)<-c("Phyla","Reads")

# Sort descending
A3<-A3[order(-A3$Reads),]

# Keep top 10
A3<-A3[1:10,]

# Create phylum summary table
phylumTable<-merge(A2,A3,by="Phyla")

# Create long form form ggplot
phylumTable.long<-melt(phylumTable, id=c("Phyla"))

# Create 100% stacked horizontal bar plot to summarize phyla diversity
p<-ggplot(data=phylumTable.long, aes(x=variable, y=value, fill=Phyla)) +
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
  annotate("text",x=1, y=0.15,label="131,269") +
  annotate("text",x=1, y=0.60,label="277,216") +
  annotate("text",x=1, y=0.96,label="40,982") +
  # annotate top 3 Reads
  annotate("text",x=2, y=0.15,label="5,824,520") +
  annotate("text",x=2, y=0.60,label="7,326,559") +
  annotate("text",x=2, y=0.96,label="2,692,708") 

ggsave("FigS1_phylumdist.pdf")

# Select phylum Arthropoda only
Arth_df<-A[A$phylum=="Arthropoda",]

# Get Arthropoda ESVs counts for Table S3
length(unique(Arth_df$marker_OTU))
#3,598
length(unique(Arth_df$marker_OTU[Arth_df$marker=="BE"]))
#775
length(unique(Arth_df$marker_OTU[Arth_df$marker=="F230"]))
#2,823

# Count number of reads in arthropod ESVs for Table S3
sum(Arth_df$reads)
#2,692,708
sum(Arth_df$reads[Arth_df$marker=="BE"])
#294,070
sum(Arth_df$reads[Arth_df$marker=="F230"])
#2,398,638
