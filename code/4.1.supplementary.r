#### Set up working directory ####
data.pa<-"./input/wet.input/"
meta.pa<-"./input/meta.data/"
database.pa<-"./input/database/"
options(warn=0)
#### End ####

#### Libraries ####
library(ggplot2)
library(dplyr)
library(vegan)
# library(plyr) #rbind.fill function
library(reshape)
library(ggsci)
library(tidyr) # separate function
#### End ####

#### color set up using ggsci ####
guan.color<-alpha(c("#7493cb","#f89a1d","#73cedf","#99991EFF","#BFA8CC","#60BC34","#FAD3D3","#E62187","#15B2EA","#269C79","#FFB246"),1)
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

#### extended fig.1a-e Total microbial DNA shotgun data ####
meta <- read.delim(paste(meta.pa,"meta.data.txt",sep=""))
# This is from kraken and QC combination processed in Mac
tmp<-read.delim("input/wet.input/shotgun.metagenome/shotgun.hiseq.reads.information.txt")
tmp[is.na(tmp)]<-0
tmp<-cbind(as.data.frame(tmp[-6,1]),prop.table(data.matrix(tmp[-6,-1]), 2))
shot.reads<-tmp

data<-melt(shot.reads)
names(data)<-c("Source","library_id","Percentage")
# merge with prop data come from # total shot gun human reads vs 16s
data<-merge(data,meta,by="library_id")
data<-arrange(data,Source,Percentage) # sort by Source then by Percentage
data$library_id<-factor(data$library_id,levels=unique(data[data$Source=="Bacteria",]$library_id),ordered = T) # change sample id order by bacteria percentage
data$Source<-factor(data$Source,levels=c("Bacteria","Archaea","Viruses","Human","Unassigned"))

#data<-merge(data,prop.com[prop.com$library=="DNA",],by.x="sample_id",by.y="sample_id")
# (1) contig proportion
g<-ggplot(data[!data$Study_group=="Negative Control",],aes(y = Percentage*100, x = interaction(library_id,Study_group), fill = Source))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=color[1:6])+
  labs(x ="",y="")+
  facet_wrap(~Study_group,scales ="free_x",nrow = 1 )+
  theme_classic()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=15,colour="black"),
    axis.ticks.x = element_blank(),
    legend.title=element_blank(),
    legend.text = element_text(size=15,colour="black"),
    strip.background = element_blank(),
    strip.text = element_text(size=20,colour="black") 
  )

pdf("output/extended.fig1a.pdf", width = 6, height = 3.5)
g
dev.off()


kraken <- read.delim("input/wet.input/shotgun.metagenome/sunbeam_output/classify/all_samples.tsv", skip = 1)
names(kraken) <- gsub(x = names(kraken),
                      pattern = "\\.taxa",
                      replacement = "\\")
names(kraken)[1]<-"tax_id" ## kraken is table with 1 is tax id, 2:157 are libraries, 158 is names
kraken.shot<-separate(data = kraken, col = Consensus.Lineage, 
                      into = c("kingdom","phyla","class",
                               "order","family","genus","species"),
                      sep = "\\;")
kraken.shot$kingdom <- gsub(x = kraken.shot$kingdom,
                            pattern = "\\k__",
                            replacement = "\\")
kraken.shot$phyla <- gsub(x = kraken.shot$phyla,
                          pattern = "\\ p__",
                          replacement = "\\")
kraken.shot$class <- gsub(x = kraken.shot$class,
                          pattern = "\\ c__",
                          replacement = "\\")
kraken.shot$order <- gsub(x = kraken.shot$order,
                          pattern = "\\ o__",
                          replacement = "\\")
kraken.shot$family <- gsub(x = kraken.shot$family,
                           pattern = "\\ f__",
                           replacement = "\\")
kraken.shot$genus <- gsub(x = kraken.shot$genus,
                          pattern = "\\ g__",
                          replacement = "\\")
kraken.shot$species <- gsub(x = kraken.shot$species,
                            pattern = "\\ s__",
                            replacement = "\\")
kraken.shot<-kraken.shot[,c(-1,-(76:80))]
kraken.shot<-kraken.shot%>%gather(key=library_id,value=read,-kingdom,-phyla)%>%
  dplyr::mutate(phyla=ifelse(phyla=="Bacteroidetes"|phyla=="Actinobacteria"|phyla=="Firmicutes"|phyla=="Proteobacteria",phyla,"Others"))

sdata<-left_join(kraken.shot,meta)
sdata$library_id<-factor(sdata$library_id,levels=unique(data[data$Source=="Bacteria",]$library_id),ordered = T)

g<-ggplot(sdata%>%filter(Study_group!="Negative Control"),aes(y =read, x = interaction(library_id,Study_group),fill=phyla))+
  geom_bar(stat="identity")+
  labs(x ="",y="")+
  scale_y_continuous(breaks=c(3e6),labels = c("3M"))+
  labs(title = "")+
  theme_classic()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=15,colour="black"),
    axis.ticks.x = element_blank(),
    legend.title=element_blank(),
    legend.text = element_text(size=15,colour="black"),
    strip.background = element_blank(),
    strip.text = element_text(size=20,colour="black"))+
  facet_wrap(~Study_group,scales ="free_x",nrow=1 )
pdf("output/extended.fig1c.pdf", width = 6, height = 3.5)
g
dev.off()


path<-"input/wet.input/discovery.virome/sunbeam_output/qc/log/decontam/"
name<-list.files(path,pattern = "_1")
l=c()
i=1
for (x in name) {
  y<-read.delim(paste(path,x,sep=""))
  y<-y%>%mutate(human.ratio=human38/(nonhost+host))
  y$human38<-strsplit(strsplit(x,"\\/")[[1]][length(strsplit(x,"\\/")[[1]])],"\\_")[[1]][1]
  l[[i]]<-y[,c(1,5)]
  i=i+1
}
data<-do.call(rbind,l)
mdata<-merge(data,meta[meta$library_type=="DNA"&meta$Study_group=="Month 0",],by.x="human38",by.y="library_id")
g<-ggplot(mdata,aes(x=Stool_hours_to_collection,y=human.ratio*100))+
  geom_point(colour="black",size=3,alpha=0.8)+
  geom_smooth(method='lm',formula=y~x,colour="black",linetype = "dashed")+
  annotate("text", x=102, y=95,size=10, label= expression(paste(italic(P)," = 0.04")))+
  annotate("text", x=112, y=85,size=10, label= expression(paste(italic(R)," = -0.45")),hjust = 0.63)+
  scale_x_continuous(trans='log2', expand = c(0.05,0))+
  labs(y="Human DNA\npercentage (%)",x="Hours to collection")+
  theme_classic()+
  theme( 
    axis.title = element_text(size=35,colour="black"),
    axis.text = element_text(size=30,colour = "black"),
    axis.title.y=element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
  )
pdf("output/extended.fig1b.pdf",width=10,height=7) 
g
dev.off()

cor.test(mdata$Stool_hours_to_collection,mdata$human.ratio,method = "spearman")


path <- "input/wet.input/shotgun.metagenome/subset.analysis/sunbeam_output/qc/decontam"
name<-list.files(path=path,pattern = ".metaphlan",full.names = T)

read_metaphlan<-function(x){
  y<-read.delim(x)
  return(y)
}
l<-lapply(name,read_metaphlan)
metaphlan.data<-Reduce(function(x, y) merge(x, y, all=TRUE), l)
metaphlan.data[is.na(metaphlan.data)]<-0



tmp<-metaphlan.data%>%gather(key=library_id,value=read,-species)%>%
  spread(key=species,value=read)

tmp<-as.data.frame(tmp)

rownames(tmp)<-tmp[,1]
tmp<-tmp[,-1]


#decontam function
#tmp1<-shot.meta[,1:4]
#tmp1$sample_id<-factor(tmp1$sample_id,levels=rownames(tmp),ordered = T)
#tmp1<-tmp1[order(tmp1$sample_id),]
#tmp1<-tmp1$final_dna_concentration_ng_ul
#tmp1<-ifelse(tmp1==0,0.00001,tmp1)
# 
tmp1<-meta%>%filter(Experiment=="Shotgun")%>%select(library_id,Study_group)
tmp1$library_id<-factor(tmp1$library_id,levels=rownames(tmp),ordered = T)
tmp1<-tmp1[order(tmp1$library_id),]

tmp1$is.neg<-tmp1$Study_group=="Negative Control"
tmp1<-tmp1$is.neg
library(phyloseq)
library(decontam)
#contam.freq <- isContaminant(as.matrix(tmp), method="frequency", threshold=0.5, conc=tmp1,normalize=FALSE)
contam.pre <- isContaminant(as.matrix(tmp), method="prevalence", threshold=0.5, neg=tmp1,normalize=FALSE)
contam.pre$id<-names(tmp)

tmp<-tmp%>%dplyr::select((contam.pre%>%filter(contaminant=="FALSE"))$id)
data<-tmp

sha<-data.frame(diversity(data,index = "shannon"))
names(sha)<-"shannon"
sha$library_id<-rownames(sha)
tmp<-left_join(sha,meta%>%dplyr::select(Study_group,library_id))
tmp<-tmp%>%mutate(age=ifelse(Study_group=="Month 0",0,ifelse(Study_group=="Month 1",1,4)))
g<-ggplot(tmp%>%filter(Study_group!="Negative Control"),aes(x=Study_group,y=shannon,fill=Study_group))+
  geom_boxplot()+
  scale_fill_manual(values = guan.color[1:3])+
  xlab("")+
  ylab("Bacterial shannon diversity")+
  annotate("text", x=c(1.5,2.5,2), y=c(2,2,2.3),size=5, label= c(expression(paste(italic(P)," = 0.0003")),expression(paste(italic(P)," = 0.39")),expression(paste(italic(P)," = 0.0002"))))+
  annotate("segment", x = c(1.05,2.05,1.05), xend = c(1.95,2.95,2.95), y = c(1.9,1.9,2.2), yend = c(1.9,1.9,2.2))+
  
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(size=20,color="black"),
        axis.text = element_text(size=15,color="black"))
pdf("output/extended.fig1e.pdf",width = 5,height = 5)
g
dev.off()


rich<-data.frame(rowSums(data>0))
names(rich)<-"richness"
rich$library_id<-rownames(rich)
tmp<-left_join(rich,meta%>%dplyr::select(Study_group,library_id))
tmp<-tmp%>%mutate(age=ifelse(Study_group=="Month 0",0,ifelse(Study_group=="Month 1",1,4)))
g<-ggplot(tmp%>%filter(Study_group!="Negative Control"),aes(x=Study_group,y=richness,fill=Study_group))+
  geom_boxplot()+
  scale_fill_manual(values = guan.color[1:3])+
  xlab("")+
  ylab("Bacterial richness")+
  annotate("text", x=c(1.5,2.5,2), y=c(9,9,9.5),size=5, label= c(expression(paste(italic(P)," = 0.0001")),expression(paste(italic(P)," = 0.59")),expression(paste(italic(P)," = 0.0002"))))+
  annotate("segment", x = c(1.05,2.05,1.05), xend = c(1.95,2.95,2.95), y = c(8.8,8.8,9.3), yend = c(8.8,8.8,9.3))+
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(size=20,color="black"),
        axis.text = element_text(size=15,color="black"))
pdf("output/extended.fig1d.pdf",width = 5,height = 5)
g
dev.off()
#### End ####

#### extended fig.4 Herv reads####
#### Data of Herv, Sine, and Line reads ratio for Human 
meta <- read.delim(paste(meta.pa,"meta.data.txt",sep=""))
tmp1<-read.delim(paste(database.pa,"/hg38/sine.id",sep=""),header=F)
tmp2<-read.delim(paste(database.pa,"/hg38/line.id",sep=""),header=F)
tmp3<-read.delim(paste(database.pa,"/hg38/erv.id",sep=""),header=F)
tmp1$re<-"sine"
tmp2$re<-"line"
tmp3$re<-"herv"
repeatid<-rbind(tmp1,tmp2,tmp3)
names(repeatid)<-c("id","re")

names<-list.files(paste(data.pa,"discovery.virome/my.analysis/contigs",sep=""))
l<-c()
j=1
pa1<-paste(data.pa,"discovery.virome/my.analysis/host.reads/repeat.count/",sep="")
pa2<-paste(data.pa,"discovery.virome/my.analysis/host.reads/total.human.reads/",sep="")
for (i in names){
  tryCatch({
    datapath <- file.path(paste(pa1,i,".count",sep=""))
    repeatc <- read.delim(datapath,header=F)
    names(repeatc)<-c("id","count")
    repeatcount<-merge(repeatc,repeatid,by="id")
    repeatc<-aggregate(.~re,repeatcount[,c("count","re")],sum)
    names(repeatc)<-c("re",i)
    datapath <- file.path(paste(pa2,i,".reads",sep=""))
    hervt <- read.delim(datapath,header=F)
    repeatc[4,i]<-(sum(hervt$V3)+sum(hervt$V4))
    repeatc[4,1]<-c("total")
    l[[j]] <- repeatc
    j=j+1
  }, error=function(e){})
}
repeatp<-Reduce(function(x, y) merge(x, y, all=TRUE,by="re"),l)
repeatp<-as.data.frame(t(repeatp))
names(repeatp)<-as.vector(t(repeatp[1,])) # use the first row as header
repeatp<-repeatp[-1,]
repeatp$herv<-as.numeric(as.character(repeatp$herv))
repeatp$line<-as.numeric(as.character(repeatp$line))
repeatp$sine<-as.numeric(as.character(repeatp$sine))
repeatp$total<-as.numeric(as.character(repeatp$total))
repeatp$library_id<-row.names(repeatp)
repeatp<-mutate(repeatp,
                hervp=herv/total,
                linep=line/total,
                sinep=sine/total)


repeatp<-merge(repeatp,meta,by="library_id")


library(ggpubr)

g<-ggbarplot(repeatp%>%filter(Study_group!="Germ-free Mice"), x = "Study_group", y = "hervp", 
             add = c("mean_se", "jitter"),
             color="library_type",
             xlab = FALSE,
             ylab=c("Herv reads percentage"),
             palette = c(guan.color[5:6]),
             position = position_dodge(0.8)
)
p<-ggpar(g,legend.title = "")

pdf("output/extended.fig4.pdf",width =5,height = 3,useDingbats = F)
p
dev.off()
#### End ####




# rm(list=ls())
# options(warn=0)
