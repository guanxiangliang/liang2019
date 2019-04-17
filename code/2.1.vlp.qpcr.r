#### Set up working directory ####
data.pa<-"./input/wet.input/"
meta.pa<-"./input/meta.data/"
options(warn=-1)
#### End ####

#### Libraries ####
library(ggplot2)
library(ggbeeswarm)
library(ggsci)
library(dplyr)
library(reshape)
#### End ####

#### color set up using ggsci ####
# all colror
color=pal_ucscgb(alpha = 1)(20)
# Age group color(three colors)# add 7 for older 
age.col=color[c(5,2,17)]
# Feeding type color
feed.col=c(color[c(11,12,13)],"grey80")
# Cohort color
cohort.col=color[c(20,4)]
#color tron
color.tron=pal_tron(alpha = 1)(7)
#color nature used with the data
#### End ####

#### Fig.1B VLP number ####
meta <- read.delim(paste(meta.pa,"meta.data.txt",sep=""))
data <- read.delim(paste(data.pa,"vlp.qpcr.txt",sep=""))
tmp<-meta[,c("study_group","library_id")]
tmp<-merge(data,tmp,all.x=T,by="library_id")
tmp[tmp$count_per_gram<6.67E+06,"count_per_gram"]<-0 # VLP numbers lower than 1 count per field were removed
g1<-ggplot(tmp,aes(x=study_group,y=count_per_gram+1))+
  geom_violin(aes(fill=study_group))+
  geom_beeswarm(data=. %>% filter(study_group==c("Month 0"),count_per_gram==0),cex=1,alpha=0.8,size=3,shape=21,colour="white",fill="black",groupOnX = T)+
  geom_beeswarm(data=. %>% filter(study_group==c("Month 0"),count_per_gram>0),cex=1,alpha=0.8,size=3,shape=21,colour="white",fill="black",groupOnX = T)+
  geom_beeswarm(data=. %>% filter(study_group==c("Month 1")),cex=4,alpha=0.8,size=3,shape=21,colour="white",fill="black",groupOnX = T)+
  geom_beeswarm(data=. %>% filter(study_group==c("Month 4")),cex=4,alpha=0.8,size=3,shape=21,colour="white",fill="black",groupOnX = T)+
  scale_y_continuous(trans='log10', breaks=c(1,1e4,1e8,1e12),labels = c("0",expression(10^4),expression(10^8),expression(10^12)), expand = c(0.05,0))+
  expand_limits(y = c(1, 1e12))+
  scale_x_discrete(breaks=c("Month 0","Month 1","Month 4"),
                   labels=c("0", "1", "4"))+
  scale_fill_manual(values = age.col)+
  ylab("VLP counts\nper gram feces")+
  xlab("Month")+
  annotate("text", x=c(1.5,2.5,2), y=c(1.5e11,1.5e11,3e12),size=5, label= c(expression(paste(italic(P)," < 0.0001")),"NS",expression(paste(italic(P)," < 0.0001"))))+
  annotate("segment", x = c(1.05,1.05,2.05), xend = c(1.95,2.95,2.95), y = c(5e10,1e12,5e10), yend = c(5e10,1e12,5e10))+
  theme_classic()+
  theme(axis.title = element_text(size=35,colour="black"),
        axis.text = element_text(size=30,colour = "black"),
        legend.position = "none"
  )
ggsave(filename="output/Fig.1B.pdf", width = 6, height = 5)
#### End ####

#### Fig.1C total 16s copy number ####
meta <- read.delim(paste(meta.pa,"meta.data.txt",sep=""))
data <- read.delim(paste(data.pa,"vlp.qpcr.txt",sep=""))
tmp<-meta[,c("study_group","library_id")]
tmp<-merge(data,tmp,all.x=T,by="library_id")
g1<-ggplot(tmp,aes(x=study_group,y=total_copy_per_gram+1))+
  geom_violin(aes(fill=study_group))+
  geom_beeswarm(data=. %>% filter(total_copy_per_gram==0),cex=1,size=3,shape=21,colour="white",fill="black",alpha=0.8,groupOnX =T )+
  geom_beeswarm(data=. %>% filter(total_copy_per_gram>0),cex=4,size=3,shape=21,colour="white",fill="black",alpha=0.8,groupOnX = T)+
  scale_y_continuous(trans='log10', breaks=c(1,1e4,1e8,1e12),labels = c("0",expression(10^4),expression(10^8),expression(10^12)), expand = c(0.05,0))+
  expand_limits(y = c(1, 1e12))+
  scale_x_discrete(breaks=c("Month 0","Month 1","Month 4"),
                   labels=c("0", "1", "4"))+
  ylab("16S copy number\nper gram feces")+
  xlab("Month")+
  scale_fill_manual(values = age.col)+
  annotate("text", x=c(1.5,2.5,2), y=c(1.5e11,1.5e11,3e12),size=5, label= c(expression(paste(italic(P)," < 0.0001")),"NS",expression(paste(italic(P)," = 0.0003"))))+
  annotate("segment", x = c(1.05,1.05,2.05), xend = c(1.95,2.95,2.95), y = c(5e10,1e12,5e10), yend = c(5e10,1e12,5e10))+
  theme_classic()+
  theme(axis.title = element_text(size=35,colour="black"),
        axis.text = element_text(size=30,colour = "black"),
        legend.position = "none"
  )
ggsave("output/Fig.1C.pdf",width = 6, height = 5)
#### End ####

#### Fig.S5 16s qPCR removal ####
meta <- read.delim(paste(meta.pa,"meta.data.txt",sep=""))
data <- read.delim(paste(data.pa,"vlp.qpcr.txt",sep=""))
tmp<-meta[,c("study_group","library_id")]
tmp<-merge(data,tmp,all.x=T)
mtmp<-melt(tmp[,c("library_id","vlp_copy_per_gram","total_copy_per_gram","study_group")])
mtmp$variable<-gsub(mtmp$variable,pattern = "total_copy_per_gram",replacement = "Before VLP purification")
mtmp$variable<-gsub(mtmp$variable,pattern = "vlp_copy_per_gram",replacement = "After VLP purification")
mtmp$variable<-factor(mtmp$variable,levels=c("Before VLP purification","After VLP purification"))
g1<-ggplot(mtmp,aes(x=study_group,y=value+1,color=variable,shape=variable))+
  geom_point(size =3,position=position_jitterdodge())+
  scale_y_continuous(breaks=c(1,101,10001,1000001,100000001,10000000001),labels = c("0",expression(10^2),expression(10^4),expression(10^6),expression(10^8),expression(10^10)), expand = c(0.05,0))+
  
  scale_x_discrete(breaks=c("Month 0","Month 1","Month 4"),
                   labels=c("0", "1", "4"))+
  ylab("16S copy number\nper gram feces")+
  xlab("Month")+
  coord_trans(y = "log10")+
  stat_summary(fun.y = mean, geom = "crossbar",mapping=aes(ymin=..y.., ymax=..y..),
               width = 1,position = position_dodge())+
  annotate("text", x=c(1,2,3), y=c(3e11,3e11,3e11),size=5, label= c(expression(paste(italic(P)," = 0.34")),expression(paste(italic(P)," < 0.0001")),expression(paste(italic(P)," = 0.002"))))+
  annotate("segment", x = c(0.8,1.8,2.8), xend = c(1.2,2.2,3.2), y = c(1e11,1e11,1e11), yend = c(1e11,1e11,1e11))+
  theme_classic()+
  theme(axis.title = element_text(size=35,colour="black"),
        axis.text = element_text(size=30,colour = "black"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=20,colour = "black")
  )
ggsave("output/Fig.S5.pdf",width=10,height=5)
wilcox.test(mtmp[mtmp$variable=="After VLP purification"&mtmp$study_group=="Month 0",]$value,mtmp[mtmp$variable=="Before VLP purification"&mtmp$study_group=="Month 0",]$value)
wilcox.test(mtmp[mtmp$variable=="After VLP purification"&mtmp$study_group=="Month 1",]$value,mtmp[mtmp$variable=="Before VLP purification"&mtmp$study_group=="Month 1",]$value)
wilcox.test(mtmp[mtmp$variable=="After VLP purification"&mtmp$study_group=="Month 4",]$value,mtmp[mtmp$variable=="Before VLP purification"&mtmp$study_group=="Month 4",]$value)


#### End ####

#### Fig.3D adeno qPCR ####
meta <- read.delim(paste(meta.pa,"meta.data.txt",sep=""))
data <- read.delim(paste(data.pa,"adeno.qpcr.txt",sep=""))
tmp<-meta[,c("Infant_feeding_type","library_id","Cohort")]
tmp<-merge(data,tmp,all.x=T,by="library_id")
tmp<-tmp %>% filter(Cohort=="Validation")
tmp$Infant_feeding_type<-gsub(tmp$Infant_feeding_type,pattern ="Mixed",replacement = "Breastmilk")
tmp$inf<-"Uninfected"
tmp[tmp$copy_number_rxn>=5,"inf"]<-"Infected"
tmp<-as.data.frame(tmp %>% count(inf,Infant_feeding_type,Cohort))
tmp1<-tmp%>%filter(Infant_feeding_type=="Formula")
tmp2<-tmp%>%filter(Infant_feeding_type=="Breastmilk")
tmp1<- tmp1%>% mutate(freq=n/sum(n))
tmp2<- tmp2%>% mutate(freq=n/sum(n))
tmp<-rbind(tmp1,tmp2)
ci.f<-binom.test(tmp[tmp$Infant_feeding_type=="Formula","n"][1],sum(tmp[tmp$Infant_feeding_type=="Formula","n"]))$conf[2] # calvculate 95% confidence interval proportion
ci.b<-binom.test(tmp[tmp$Infant_feeding_type=="Breastmilk","n"][1],sum(tmp[tmp$Infant_feeding_type=="Breastmilk","n"]))$conf[2] # calvculate 95% confidence interval proportion
cil.f<-binom.test(tmp[tmp$Infant_feeding_type=="Formula","n"][1],sum(tmp[tmp$Infant_feeding_type=="Formula","n"]))$conf[1] # calvculate 95% confidence interval proportion
cil.b<-binom.test(tmp[tmp$Infant_feeding_type=="Breastmilk","n"][1],sum(tmp[tmp$Infant_feeding_type=="Breastmilk","n"]))$conf[1] # calvculate 95% confidence interval proportion
tmp<-tmp %>% filter(inf=="Infected")
tmp$ci<-NA
tmp[tmp$Infant_feeding_type=="Formula",]$ci<-ci.f
tmp[tmp$Infant_feeding_type=="Breastmilk",]$ci<-ci.b
tmp$cil<-NA
tmp[tmp$Infant_feeding_type=="Formula",]$cil<-cil.f
tmp[tmp$Infant_feeding_type=="Breastmilk",]$cil<-cil.b
tmp$Infant_feeding_type<-factor(tmp$Infant_feeding_type,levels=c("Formula","Breastmilk"))

g1<-ggplot(tmp,aes(y = freq*100, x = Infant_feeding_type,fill=Infant_feeding_type)) +
  
  geom_bar( color="black",stat="identity",width = 0.5)+
  geom_errorbar(aes(ymin=cil*100, ymax=ci*100), width=0.2,size=0.5,
                position=position_dodge(.9)) +
  scale_fill_manual(values=feed.col[1:2])+
  expand_limits(y = c(0, 37))+
  labs(title="",x ="", y = "Percentage of subjects (%)\nwith Adenoviridae")+
  scale_x_discrete(labels=c( "Formula","Breastmilk\n+ Mixed"))+
  facet_grid(.~Cohort,scales = "free", space = "free")+
  #annotate("text", x = c(1.5), y = c(37), size=5, label= c(expression(paste("n = 86, ",italic(P)," = 0.049"))))+
  #annotate("segment", x = c(1.05), xend = c(1.95), y = c(36), yend = c(36))+
  theme_classic()+
  theme(axis.title = element_text(size=25,colour="black"),
        axis.text.x = element_text(size=20,colour = "black"),
        axis.text.y = element_text(size=20,colour = "black",margin = margin(t = -20, unit = "pt")),
        strip.text = element_text(size=25),
        strip.background = element_rect(fill=cohort.col[2]),
        legend.position = "none"
  )

ggsave("output/Fig.3D.pdf", width = 5, height = 6)

#### End ####

# rm(list=ls())
# options(warn=0)




