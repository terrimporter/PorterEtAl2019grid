# Teresita M. Porter, July 8, 2019

library(vegan)
library(reshape2)
library(gridExtra)
library(grid)
library(ggplot2)
library(plyr)
library(scales)

###################################################################

# Read infile
A<-read.table(file="LV2016_2.csv", head=TRUE, sep=",")

# Select phylum Arthropoda only
B<-A[A$phylum=="Arthropoda",]

# Calc total taxa and total taxa confidently id'd
species<-length(unique(B$species))
species_good<-length(unique(B$species[B$sBP>=0.70]))
genus<-length(unique(B$genus))
genus_good<-length(unique(B$genus[B$gBP>=0.30]))
family<-length(unique(B$family))
family_good<-length(unique(B$family[B$fBP>=0.20]))

# create df for ggplot
df<-data.frame("rank"=c("species","species","genus","genus","family","family"),
               "status"=rep(c("all","good"),3),
               "value"=c(species, species_good, genus, genus_good, family, family_good))

# create factors
df$rank = factor(df$rank, levels=c("species","genus","family"), 
                 labels=c("Species","Genus","Family"))
df$status = factor(df$status, levels=c("all","good"),
                   labels=c("All taxa","Confidently identified taxa"))

# create bar plot with two series
p<-ggplot(df, aes(fill=status, y=value, x=rank)) +
  geom_bar(position="dodge",stat="identity") +
  scale_x_discrete(limits = rev(levels(rank))) +
  labs(x="Rank", y="Unique Arthropod Taxa") +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title=element_blank(),
        legend.text=element_text(size=14),
        legend.position = "bottom",
        legend.spacing.x = unit(0.5, 'cm'),
        axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        axis.title.x=element_blank())
ggsave("FigS3_confidentids.pdf")
