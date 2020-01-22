#### Set up working directory ####
data.pa<-"./input/wet.input/"
meta.pa<-"./input/meta.data/"
options(warn=-1)
#### End ####

#### Libraries ####
library(ggplot2)
library(dplyr)
library(ggbeeswarm)
library(tidyverse)
#### End ####

#### color set up using ggsci ####
guan.color<-alpha(c("#7493cb","#f89a1d","#73cedf","#99991EFF","#BFA8CC","#60BC34","#FAD3D3","#E62187","#15B2EA","#269C79","#FFB246"),1)
#barplot(rep(1,11),rep(1,11),col=guan.color,space=0,horiz = T,axes = F,xlim=c(0,5))
#text(x=rep(0.5,11),y=seq(0.5,10.5,by=1),1:11,cex=1)
#### End ####

#### Fig.1b VLP number ####
meta<-read.delim(paste0(meta.pa,"meta.data.txt"),colClasses = c(rep("character",12),rep("numeric",4),"character",rep("numeric",2),rep("character",4)),stringsAsFactors = T)
other<-read.delim(paste0(data.pa,"other.data.txt"),colClasses =c(rep("character",2),rep("numeric",9)))
tmp<-other%>%left_join(meta%>%filter(Experiment=="Shotgun"))

g<-ggplot(tmp,aes(x=Study_group,y=(VLP_count_per_gram+100000)/100000))+
  geom_violin(aes(fill=Study_group))+
  geom_beeswarm(data=. %>% filter(Study_group=="Month 0",VLP_count_per_gram==0),cex=1,alpha=0.8,size=3,shape=21,colour="white",fill="black",groupOnX = T)+
  geom_beeswarm(data=. %>% filter(Study_group=="Month 0",VLP_count_per_gram>0),cex=1,alpha=0.8,size=3,shape=21,colour="white",fill="black",groupOnX = T)+
  geom_beeswarm(data=. %>% filter(Study_group=="Month 1"),cex=4,alpha=0.8,size=3,shape=21,colour="white",fill="black",groupOnX = T)+
  geom_beeswarm(data=. %>% filter(Study_group=="Month 4"),cex=4,alpha=0.8,size=3,shape=21,colour="white",fill="black",groupOnX = T)+
  #geom_hline(yintercept = 6e1,linetype="dashed")+
  scale_y_continuous(trans='log10', breaks=c(1e2,1e4,1e6),labels = c(expression(10^7),expression(10^9),expression(10^11)), expand = c(0.05,0))+
  expand_limits(y = c(1, 1e7))+
  scale_x_discrete(breaks=c("Month 0","Month 1","Month 4"),
                   labels=c("0", "1", "4"))+
  scale_fill_manual(values = guan.color[1:3])+
  #annotation_logticks(sides = "l")+
  ylab("VLP counts\nper gram feces")+
  xlab("Month")+
  annotate("text", x=c(1.5,2.5,2,3.4), y=c(1.5e6,1.5e6,3e7,35),size=5, label= c(expression(paste(italic(P)," < 0.0001")),"NS",expression(paste(italic(P)," < 0.0001")),"LoQ"))+
  annotate("segment", x = c(1.05,1.05,2.05), xend = c(1.95,2.95,2.95), y = c(5e5,1e7,5e5), yend = c(5e5,1e7,5e5))+
  theme_classic()+
  theme(axis.title = element_text(size=35,colour="black"),
        axis.text = element_text(size=30,colour = "black"),
        legend.position = "none"
  )
pdf("output/Fig.1b.pdf", width = 6, height = 5)
g
dev.off()

wilcox.test(tmp[tmp$Study_group=="Month 0","VLP_count_per_gram"],tmp[tmp$Study_group=="Month 1","VLP_count_per_gram"],paired = T)
wilcox.test(tmp[tmp$Study_group=="Month 0","VLP_count_per_gram"],tmp[tmp$Study_group=="Month 4","VLP_count_per_gram"],paired = T)
wilcox.test(tmp[tmp$Study_group=="Month 1","VLP_count_per_gram"],tmp[tmp$Study_group=="Month 4","VLP_count_per_gram"],paired = T)
#### End ####

#### Fig.1c total 16s copy number ####
meta<-read.delim(paste0(meta.pa,"meta.data.txt"),colClasses = c(rep("character",12),rep("numeric",4),"character",rep("numeric",2),rep("character",4)),stringsAsFactors = T)
other<-read.delim(paste0(data.pa,"other.data.txt"),colClasses =c(rep("character",2),rep("numeric",9)))
tmp<-other%>%left_join(meta%>%filter(Experiment=="Shotgun"))

g<-ggplot(tmp,aes(Study_group,Total_microbal_DNA_16s_copy_per_gram+1))+
  geom_violin(aes(fill=Study_group))+
  geom_hline(yintercept = 2e3,linetype="dashed")+
  geom_beeswarm(data=. %>% filter(Total_microbal_DNA_16s_copy_per_gram==0),cex=1,size=3,shape=21,colour="white",fill="black",alpha=0.8,groupOnX =T )+
  geom_beeswarm(data=. %>% filter(Total_microbal_DNA_16s_copy_per_gram>0),cex=4,size=3,shape=21,colour="white",fill="black",alpha=0.8,groupOnX = T)+
  scale_y_continuous(trans='log10', breaks=c(1e3,1e6,1e9,1e12),labels = c(expression(10^3),expression(10^6),expression(10^9),expression(10^12)), expand = c(0.05,0))+
  expand_limits(y = c(1, 1e12))+
  scale_x_discrete(breaks=c("Month 0","Month 1","Month 4"),
                   labels=c("0", "1", "4"))+
  ylab("16S copy number\nper gram feces")+
  xlab("Month")+
  scale_fill_manual(values = guan.color[1:3])+
  annotate("text", x=c(1.5,2.5,2,3.4), y=c(1.5e11,1.5e11,3e12,800),size=5, label= c(expression(paste(italic(P)," < 0.0001")),"NS",expression(paste(italic(P)," = 0.0003")),"LoQ"))+
  annotate("segment", x = c(1.05,1.05,2.05), xend = c(1.95,2.95,2.95), y = c(5e10,1e12,5e10), yend = c(5e10,1e12,5e10))+
  theme_classic()+
  theme(axis.title = element_text(size=35,colour="black"),
        axis.text = element_text(size=30,colour = "black"),
        legend.position = "none"
  )
pdf("output/Fig.1c.pdf", width = 6, height = 5)
g
dev.off()

wilcox.test(tmp[tmp$Study_group=="Month 0","Total_microbal_DNA_16s_copy_per_gram"],tmp[tmp$Study_group=="Month 1","Total_microbal_DNA_16s_copy_per_gram"],paired = T)
wilcox.test(tmp[tmp$Study_group=="Month 0","Total_microbal_DNA_16s_copy_per_gram"],tmp[tmp$Study_group=="Month 4","Total_microbal_DNA_16s_copy_per_gram"],paired = T)
wilcox.test(tmp[tmp$Study_group=="Month 1","Total_microbal_DNA_16s_copy_per_gram"],tmp[tmp$Study_group=="Month 4","Total_microbal_DNA_16s_copy_per_gram"],paired = T)


#### End ####

#### Extended Data Fig.10 16s qPCR removal ####
meta<-read.delim(paste0(meta.pa,"meta.data.txt"),colClasses = c(rep("character",12),rep("numeric",4),"character",rep("numeric",2),rep("character",4)),stringsAsFactors = T)
other<-read.delim(paste0(data.pa,"other.data.txt"),colClasses =c(rep("character",2),rep("numeric",9)))
tmp<-other%>%left_join(meta%>%filter(Experiment=="Shotgun"))%>%
  dplyr::select(Study_group,VLP_16s_copy_per_gram,Total_microbal_DNA_16s_copy_per_gram)%>%
  dplyr::rename(Before_VLP_purification=Total_microbal_DNA_16s_copy_per_gram)%>%
  dplyr::rename(After_VLP_purification=VLP_16s_copy_per_gram)%>%
  gather(key=variable,value=value,-Study_group)

tmp$variable<-factor(tmp$variable,levels=c("Before_VLP_purification","After_VLP_purification"))
  
g<-ggplot(tmp,aes(x=Study_group,y=value+1,color=variable,shape=variable))+
  geom_point(size =3,position=position_jitterdodge())+
  scale_y_continuous(breaks=c(1,101,10001,1000001,100000001,10000000001),labels = c("0",expression(10^2),expression(10^4),expression(10^6),expression(10^8),expression(10^10)), expand = c(0.05,0))+
  
  scale_x_discrete(breaks=c("Month 0","Month 1","Month 4"),
                   labels=c("0", "1", "4"))+
  ylab("16S copy number\nper gram feces")+
  xlab("Month")+
  coord_trans(y = "log10")+
  stat_summary(fun.y = mean, geom = "crossbar",mapping=aes(ymin=..y.., ymax=..y..),
               width = 1,position = position_dodge())+
  annotate("text", x=c(1,2,3), y=c(3e11,3e11,3e11),size=5, label= c(expression(paste(italic(P)," = 0.68")),expression(paste(italic(P)," = 0.0009")),expression(paste(italic(P)," = 0.0003"))))+
  annotate("segment", x = c(0.8,1.8,2.8), xend = c(1.2,2.2,3.2), y = c(1e11,1e11,1e11), yend = c(1e11,1e11,1e11))+
  theme_classic()+
  theme(axis.title = element_text(size=35,colour="black"),
        axis.text = element_text(size=30,colour = "black"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=20,colour = "black")
  )
pdf("output/Extended.Data.Fig10.pdf",width=10,height=5)
g
dev.off()

wilcox.test(tmp[tmp$variable=="After_VLP_purification"&tmp$Study_group=="Month 0",]$value,tmp[tmp$variable=="Before_VLP_purification"&tmp$Study_group=="Month 0",]$value,paired = T)
wilcox.test(tmp[tmp$variable=="After_VLP_purification"&tmp$Study_group=="Month 1",]$value,tmp[tmp$variable=="Before_VLP_purification"&tmp$Study_group=="Month 1",]$value,paired = T)
wilcox.test(tmp[tmp$variable=="After_VLP_purification"&tmp$Study_group=="Month 4",]$value,tmp[tmp$variable=="Before_VLP_purification"&tmp$Study_group=="Month 4",]$value,paired = T)


#### End ####







# rm(list=ls())
# options(warn=0)




