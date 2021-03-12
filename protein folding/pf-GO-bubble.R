rm(list = ls())
setwd("D:/ArchR/protein folding")
a<-read.csv("GO.csv")
head(a)
#a<-a[,-6:-7]
#########bubble plot############
library(ggplot2)
library(RColorBrewer)
library(ggrepel)

ggplot(a,aes(x=celltype,y=Description))+
  geom_point(aes(size=Count,color=-log10(p.adjust)),alpha=0.6)+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())+
  scale_color_gradient(low = "blue",high = "red")



ggplot(q,aes(x=cell,y=Description))+
  geom_point(aes(size=Count,color=-log10(p.adjust)),alpha=0.6)+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())+
  scale_color_gradient(low = "blue",high = "red")

