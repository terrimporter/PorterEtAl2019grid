# Teresita M. Porter, Feb. 13, 2019

library(vegan)
library(reshape2)
library(gridExtra)
library(grid)
library(ggplot2)
library(plyr)

###################################################################

# Read infile prepared by python script
## LV2016_2.csv is the taxonomic assignment table
A<-read.table(file="LV2016_2.csv", head=TRUE, sep=",")

# Select phylum Arthropoda only
Arth_df<-A[A$phylum=="Arthropoda",]

## Summarize order data
order_ESVs<-aggregate(Arth_df$marker_OTU, by=list(Arth_df$order, Arth_df$layer, Arth_df$site), FUN=length)

## Rename columns
names(order_ESVs)<-c("order","layer","site","ESVs")

## Get uniques
unique_orders<-unique(order_ESVs$order)

# Create factors
order_ESVs$order<-factor(order_ESVs$order, levels=rev(unique_orders))
order_ESVs$layer<-factor(order_ESVs$layer, levels=c("B","O","M"), labels=c("Bryophyte","Organic","Mineral"))
order_ESVs$site<-factor(order_ESVs$site, levels=c("ILC45", "NZC85"), labels=c("Island Lake", "Nimitz"))

# Horizontal stacked bar plot with facet wrap
o<-ggplot(data=order_ESVs, aes(x=order, y=ESVs, fill=layer)) +
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~site) +
  coord_flip() +
  labs(y="ESVs", x="Arthropod orders") +
  scale_fill_manual(name="Layer",values=c("#4DAF4A","#377EB8","#FF7F00")) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=10),
        legend.title = element_blank(),
        legend.text = element_text(size=12, margin = margin(r = 10, unit = "pt")),
        legend.position = "bottom", 
        strip.text = element_text(size = 14))
ggsave("SF2_order.pdf", width=8.5, height=11, units="in")
