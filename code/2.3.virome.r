#### Set up working directory ####
data.pa<-"./input/wet.input/"
meta.pa<-"./input/meta.data/"
database.pa<-"./input/database/"
uniprottax<-uniprottax<-read.delim(paste(database.pa,"uniprot.virome/uniprot.virome.meta",sep=""),header=F)
names(uniprottax)<-c("uniprot.hit.id","taxid","name")
uniprottax$taxid<-as.numeric(as.character(uniprottax$taxid))
host <- read.delim(paste0(meta.pa,"viral.host.txt")) # virus host data
options(warn=-1)
#### End ####

#### Libraries ####
library(ggplot2)
library(tidyr) # separate function
library(taxonomizr)### convert accession id to tax id and to names ###
library(dplyr)
library(ggbeeswarm)
library(RColorBrewer)
library(grid)
library(pheatmap)
library(reshape)
library(vegan) 
library(ape)



# library(gtable) # need for gtable filter
# 
# 

#### End ####

#### color set up using ggsci ####
guan.color<-alpha(c("#7493cb","#f89a1d","#73cedf","#99991EFF","#BFA8CC","#60BC34","#FAD3D3","#E62187","#15B2EA","#269C79","#FFB246"),1)
# # all colror
# color=pal_ucscgb(alpha = 1)(20)
# # Age group color(three colors)# add 7 for older 
# age.col=color[c(5,2,17)]
# # Feeding type color
# feed.col=c(color[c(11,12,13)],"grey80")
# # Cohort color
# cohort.col=color[c(20,4)]
# #color tron
# color.tron=pal_tron(alpha = 1)(7)
# #color nature used with the data
#### End ####

#### Fig.1d discovery virome richness ####
meta<-read.delim(paste0(meta.pa,"meta.data.txt"),colClasses = c(rep("character",12),rep("numeric",4),"character",rep("numeric",2),rep("character",4)),stringsAsFactors = T)
contig.nt<-read.delim("input/wet.input/discovery.virome/my.analysis/cross.assembly/contig.nt")
names(contig.nt)<-c("contig.id", "nt.top.hit", "pident", "qlen","slen","alignment.len" ,"mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
contig.nt<-contig.nt%>%
  distinct(contig.id,.keep_all = T)%>%
  separate(col=nt.top.hit,into=c("1", "2","3","ntname"),sep="\\|")
taxaid<-accessionToTaxa(contig.nt$ntname,"input/database/ac2tax/accessionTaxa.sql")
contig.nt<-cbind(contig.nt,getTaxonomy(taxaid,"input/database/ac2tax/accessionTaxa.sql"))
contig.nt<-contig.nt%>%dplyr::select(contig.id,pident,qlen,slen,alignment.len,superkingdom,family)
contig.nt$superkingdom<-as.character(contig.nt$superkingdom)


# DNA contig
pa="input/wet.input/discovery.virome/my.analysis/cross.assembly"
name<-list.files(pa, full.names = T,pattern = "D")
read_cross<-function(x){
  # filter viral contig
  id<-strsplit(x,"\\/")[[1]][length(strsplit(x,"\\/")[[1]])]
  # circularity
  data<-read.delim(paste(x,"/",id,".contig.3k.cir",sep=""),sep="\t",header=T)
  # total orf number
  tmp<-read.delim(paste(x,"/",id,".orf.number",sep=""),sep="\t",header=F)
  tmp<-separate(data = tmp, col = V1, into = c("contig.id"), sep = "\\_")
  tmp<-data.frame(table(tmp$contig))
  names(tmp)<-c("contig.id","orf.number")
  data<-left_join(data,tmp)
  # uniprot
  uniprot<-read.delim(paste(x,"/",id,".3k.uniprot",sep=""),header=F)
  names(uniprot)<-c("contig.id", "uniprot.hit.id", "pident", "qlen","slen","alignment.len" ,"mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  uniprot<- uniprot  %>% 
    separate(col = contig.id, into = c("contig.id","orfs.id"), sep = "\\_") 
  uniprot<-inner_join(uniprot,uniprottax,by="uniprot.hit.id")
  uniprot<-cbind(uniprot,getTaxonomy(uniprot$taxid,"/media/lorax/users/guanxiang/1_igram/liang2019/input/database/ac2tax/accessionTaxa.sql"))
  # uniprot orf number >50% belong to virus
  tmp<-uniprot%>%group_by(contig.id)%>%distinct(contig.id,orfs.id,.keep_all= TRUE)%>%
    dplyr::summarise(uniprot.total.orf.number=n())
  data<-left_join(data,tmp)
  data[is.na(data$uniprot.total.orf.number),"uniprot.total.orf.number"]<-0
  
  tmp<-uniprot%>%distinct(contig.id,orfs.id,.keep_all = T)%>%dplyr::select("contig.id","species")
  tmp<-table(tmp)
  tmp1<-data.frame(colnames(tmp)[apply(tmp,1,which.max)])
  tmp<-uniprot%>%distinct(contig.id,orfs.id,.keep_all = T)%>%dplyr::select("contig.id","species")
  tmp<-data.frame(table(tmp))
  tmp2<-data.frame(tapply(tmp$Freq, tmp$contig.id, max))
  tmp3<-cbind(tmp1,tmp2)
  names(tmp3)<-c("uniprot.species","uniprot.species.orf.number")
  tmp3$contig.id<-rownames(tmp3)
  rownames(tmp3)<-NULL
  data<-left_join(data,tmp3,by="contig.id")
  data$uniprot.species<-as.character(data$uniprot.species)
  data$uniprot.species.orf.number<-as.numeric(data$uniprot.species.orf.number)
  data[is.na(data$uniprot.species.orf.number),"uniprot.species.orf.number"]<-0
  data<-data%>%mutate(uniprot.species=ifelse(uniprot.species.orf.number==0,"Others",uniprot.species))
  
  
  tmp<-uniprot%>%distinct(contig.id,orfs.id,.keep_all = T)%>%dplyr::select("contig.id","family")
  tmp<-table(tmp)
  tmp1<-data.frame(colnames(tmp)[apply(tmp,1,which.max)])
  tmp<-uniprot%>%distinct(contig.id,orfs.id,.keep_all = T)%>%dplyr::select("contig.id","family")
  tmp<-data.frame(table(tmp))
  tmp2<-data.frame(tapply(tmp$Freq, tmp$contig.id, max))
  tmp3<-cbind(tmp1,tmp2)
  names(tmp3)<-c("uniprot.family","uniprot.family.orf.number")
  tmp3$contig.id<-rownames(tmp3)
  rownames(tmp3)<-NULL
  data<-left_join(data,tmp3,by="contig.id")
  data$uniprot.family<-as.character(data$uniprot.family)
  data$uniprot.family.orf.number<-as.numeric(data$uniprot.family.orf.number)
  data[is.na(data$uniprot.family.orf.number),"uniprot.family.orf.number"]<-0
  data<-data%>%mutate(uniprot.family=ifelse(uniprot.family.orf.number<=0,"Others",ifelse(uniprot.family.orf.number/uniprot.total.orf.number<=0.5,"Others",uniprot.family)))
  tmp<-uniprot%>%distinct(contig.id,orfs.id,.keep_all = T)%>%dplyr::select("contig.id","family","alignment.len")
  tmp<-tmp%>%group_by(contig.id,family)%>%
    dplyr::summarise(alignment=sum(alignment.len))%>%
    arrange(desc(alignment))%>%distinct(contig.id,.keep_all = T)%>%
    dplyr::rename(len.family=family)
  data<-data%>%left_join(tmp%>%dplyr::select(contig.id,len.family),by="contig.id")%>%
    mutate(final.family=ifelse(len.family==uniprot.family,as.character(len.family),"Others"))%>%
    mutate(final.family=ifelse(is.na(final.family),"Others",final.family))
  
  #data<-data%>%filter(uniprot.total.orf.number/orf.number > 0.5) #
  # all reads
  read.id<-list.files(x, full.names = T,pattern = ".read.txt")
  bam_read<-function(y){
    p = read.delim(y, header=F,
                   stringsAsFactors = FALSE)
    p<-p%>%mutate(rate=V3/(sum(V3)+sum(V4))*1000000)
    names(p)<-c("contig.id", 
                "contig.len",
                "mapped", 
                "unmapped",
                strsplit(strsplit(y,"\\/")[[1]][length(strsplit(y,"\\/")[[1]])],"\\.")[[1]][2])
    return(p[-nrow(p),c(1,5)])}
  l<-lapply(read.id,bam_read)
  tmp<-Reduce(function(x, y) left_join(x, y),l)
  data<-left_join(data,tmp)
  #tmp1<-data%>%arrange_(paste(paste0(id,"0")))%>%slice((nrow(data)-4):nrow(data))
  #tmp2<-data%>%arrange_(paste(paste0(id,"1")))%>%slice((nrow(data)-4):nrow(data))
  #tmp3<-data%>%arrange_(paste(paste0(id,"4")))%>%slice((nrow(data)-4):nrow(data))
  #data<-rbind(tmp1,tmp2,tmp3)%>%distinct(contig.id,.keep_all = T)
  return(data)
}
l<-c()
l<-lapply(name,read_cross)
data<-do.call(rbind,l)
l<-c()
j<-c()
m<-c()
i=1
name<-list.files(pa, full.names = F,pattern = "D")
for (x in name) {
  tmp<-data%>%filter(grepl(paste0(x,".contig"),contig.id))%>%
    dplyr::select("contig.id","length","circularity","orf.number","uniprot.total.orf.number","uniprot.species","uniprot.species.orf.number","uniprot.family.orf.number","uniprot.family","len.family","final.family",contains(x),contains("DB"),contains("DD"),contains("DT"))%>%
    mutate(max1=pmax(!!as.name(paste0(x,"0")),!!as.name(paste0(x,"1")),!!as.name(paste0(x,"4"))),max2=pmax(DB1,DB2,DB3,DB4,DB5,DD1,DD2,DDH,DDN1,DDN2,DT1,DT2,DT3,DT4,DT5))%>%
    mutate(source=ifelse(max1<=max2,"Contamination",ifelse(uniprot.total.orf.number/orf.number > 0.5,"Viruses","Nonviruses")))%>%
    left_join(host[,c("family","host.name")],by=c("final.family" = "family"))%>%
    left_join(contig.nt)%>%
    mutate(coverage=alignment.len/qlen)
  tmp$host.name<-as.character(tmp$host.name)
  tmp[is.na(tmp$host.name),"host.name"]<-"Others"
  tmp$len.family<-as.character(tmp$len.family)
  tmp[is.na(tmp$len.family),"len.family"]<-"Others"
  tmp$superkingdom<-as.character(tmp$superkingdom)
  tmp[is.na(tmp$superkingdom),"superkingdom"]<-"Others"
  tmp[is.na(tmp$pident),"pident"]<-0
  tmp[is.na(tmp$coverage),"coverage"]<-0
  
  tmp<-tmp%>%mutate(new.source=ifelse(source=="Contamination","Contamination",
                                      ifelse(host.name=="Bacteriophages",source,ifelse(pident>80&coverage>0.8,superkingdom,source))))%>%
    mutate(new.source=ifelse(new.source=="Others",source,new.source))
  
  
  assign(paste0(x,".data"),tmp)
  
  read<-tmp%>%
    dplyr::select(contains(x),new.source)%>%
    gather(key=library_id,value=rpm,-new.source)%>%
    group_by(library_id,new.source)%>%
    dplyr::summarise(total=sum(rpm)/10000)
  l[[i]]<-read
  
  rich<-tmp%>%filter(new.source=="Viruses")%>%
    dplyr::select(uniprot.species,contains(x),contains("DB"),contains("DD"),contains("DT"))%>%
    gather(key=library_id,value=rpm,-uniprot.species)%>%
    group_by(uniprot.species,library_id)%>%
    dplyr::summarise(total=sum(rpm))%>%
    filter(total>10)%>%
    group_by(library_id)%>%
    dplyr::summarise(n=n())
  j[[i]]<-rich
  
  family<-tmp%>%filter(new.source=="Viruses")%>%
    dplyr::select(final.family,contains(x),contains("DB"),contains("DD"),contains("DT"))%>%
    gather(key=library_id,value=rpm,-final.family)%>%
    group_by(final.family,library_id)%>%
    dplyr::summarise(total=sum(rpm))
  
  m[[i]]<-family
  i=i+1
}

rich1<-do.call(rbind,j)


# RNA contig
pa="input/wet.input/discovery.virome/my.analysis/cross.assembly"
rna=c("R345","R352","R357","R359","R365","R368","R375","R385","R396","R400","R401","R428")
name<-paste0(pa,"/",rna)
l<-c()
l<-lapply(name,read_cross)
data<-do.call(rbind,l)

# Based on family orf length
l<-c()
j<-c()
m<-c()
i=1
name<-rna
for (x in name) {
  tmp<-data%>%filter(grepl(paste0(x,".contig"),contig.id))%>%
    dplyr::select("contig.id","length","circularity","orf.number","uniprot.total.orf.number","uniprot.species","uniprot.species.orf.number","uniprot.family.orf.number","uniprot.family","len.family","final.family",contains(x),contains("RB"),contains("RD"),contains("RT"))%>%
    mutate(max1=pmax(!!as.name(paste0(x,"0")),!!as.name(paste0(x,"1")),!!as.name(paste0(x,"4"))),max2=pmax(RB1,RB2,RB3,RB4,RB5,RD1,RD2,RDH,RDN1,RDN2,RT1,RT2,RT3,RT4,RT5))%>%
    mutate(source=ifelse(max1<=max2,"Contamination",ifelse(uniprot.total.orf.number/orf.number > 0.5,"Viruses","Nonviruses")))%>%
    left_join(host[,c("family","host.name")],by=c("len.family" = "family"))%>%
    left_join(contig.nt)%>%
    mutate(coverage=alignment.len/qlen)
  tmp$host.name<-as.character(tmp$host.name)
  tmp[is.na(tmp$host.name),"host.name"]<-"Others"
  tmp$len.family<-as.character(tmp$len.family)
  tmp[is.na(tmp$len.family),"len.family"]<-"Others"
  tmp$superkingdom<-as.character(tmp$superkingdom)
  tmp[is.na(tmp$superkingdom),"superkingdom"]<-"Others"
  tmp[is.na(tmp$pident),"pident"]<-0
  tmp[is.na(tmp$coverage),"coverage"]<-0
  
  tmp<-tmp%>%mutate(new.source=ifelse(source=="Contamination","Contamination",
                                      ifelse(host.name=="Bacteriophages",source,ifelse(pident>80&coverage>0.8,superkingdom,source))))%>%
    mutate(new.source=ifelse(new.source=="Others",source,new.source))
  
  
  assign(paste0(x,".data"),tmp)
  
  read<-tmp%>%
    dplyr::select(contains(x),new.source)%>%
    gather(key=library_id,value=rpm,-new.source)%>%
    group_by(library_id,new.source)%>%
    dplyr::summarise(total=sum(rpm)/10000)
  l[[i]]<-read
  
  
  rich<-tmp%>%filter(new.source=="Viruses")%>%
    dplyr::select(uniprot.species,contains(x),contains("RB"),contains("RD"),contains("RT"))%>%
    gather(key=library_id,value=rpm,-uniprot.species)%>%
    group_by(uniprot.species,library_id)%>%
    dplyr::summarise(total=sum(rpm))%>%
    filter(total>10)%>%
    group_by(library_id)%>%
    dplyr::summarise(n=n())
  j[[i]]<-rich
  
  family<-tmp%>%filter(new.source=="Viruses")%>%
    dplyr::select(final.family,contains(x),contains("RB"),contains("RD"),contains("RT"))%>%
    gather(key=library_id,value=rpm,-final.family)%>%
    group_by(final.family,library_id)%>%
    dplyr::summarise(total=sum(rpm))
  
  m[[i]]<-family
  i=i+1
}

rich2<-do.call(rbind,j)
rich2$library_id<-gsub(rich2$library_id,pattern = "R",replacement = "D")
rich1<-rich1%>%arrange(desc(n))%>%distinct(library_id,.keep_all = T)
rich2<-rich2%>%arrange(desc(n))%>%distinct(library_id,.keep_all = T)
tmp<-rbind(rich1,rich2)%>%arrange(desc(n))%>%group_by(library_id)%>%dplyr::summarise(n=sum(n))%>%left_join(meta)
wilcox.test(tmp[tmp$Study_group=="Month 0",]$n,tmp[tmp$Study_group=="Month 1",]$n,paired = T)
wilcox.test(tmp[tmp$Study_group=="Month 0",]$n,tmp[tmp$Study_group=="Month 4",]$n,paired = T)
wilcox.test(tmp[tmp$Study_group=="Month 1",]$n,tmp[tmp$Study_group=="Month 4",]$n,paired = T)
wilcox.test(tmp[tmp$Study_group=="Month 1",]$n,tmp[tmp$Study_group=="Negative Control",]$n,paired = F)
wilcox.test(tmp[tmp$Study_group=="Negative Control",]$n,tmp[tmp$Study_group=="Month 4",]$n,paired = F)
wilcox.test(tmp[tmp$Study_group=="Negative Control",]$n,tmp[tmp$Study_group=="Month 0",]$n,paired = F)

g<-ggplot(tmp%>%filter(Study_group!="Negative Control"),aes(Study_group,n))+
  geom_violin(aes(fill=Study_group))+
  geom_beeswarm(cex=4,size=3,shape=21,colour="white",fill="black",alpha=0.8,groupOnX =T )+
  #scale_y_continuous(trans='log10', breaks=c(1,1e4,1e8,1e12),labels = c("0",expression(10^4),expression(10^8),expression(10^12)), expand = c(0.05,0))+
  expand_limits(y = c(1, 130))+
  scale_x_discrete(breaks=c("Month 0","Month 1","Month 4"),
                   labels=c("0", "1", "4"))+
  ylab("Viral richness")+
  xlab("Month")+
  scale_fill_manual(values = guan.color[1:3])+
  annotate("text", x=c(1.5,2.5,2), y=c(115,115,125),size=5, label= c(expression(paste(italic(P)," < 0.0001")),"NS",expression(paste(italic(P)," < 0.0001"))))+
  annotate("segment", x = c(1.05,1.05,2.05), xend = c(1.95,2.95,2.95), y = c(110,120,110), yend = c(110,120,110))+
  theme_classic()+
  theme(axis.title = element_text(size=35,colour="black"),
        axis.text = element_text(size=30,colour = "black"),
        legend.position = "none"
  )
pdf("output/Fig.1d.pdf", width = 6, height = 5)
g
dev.off()

g<-ggplot(tmp,aes(Study_group,n))+
  geom_violin(aes(fill=Study_group))+
  geom_beeswarm(cex=3,size=3,shape=21,colour="white",fill="black",alpha=0.8,groupOnX =T )+
  #scale_y_continuous(trans='log10', breaks=c(1,1e4,1e8,1e12),labels = c("0",expression(10^4),expression(10^8),expression(10^12)), expand = c(0.05,0))+
  expand_limits(y = c(1, 150))+
  #scale_x_discrete(breaks=c("Month 0","Month 1","Month 4","Negative Control"),
  #                 labels=c("0", "1", "4"))+
  ylab("Viral richness")+
  xlab("")+
  scale_fill_manual(values = guan.color[1:4])+
  annotate("text", x=c(1.5,2.5,2,3,3.5,2.5), y=c(115,115,125,135,115,145),size=5, label= c(expression(paste(italic(P)," < 0.0001")),expression(paste(italic(P)," = 0.05")),expression(paste(italic(P)," < 0.0001")),expression(paste(italic(P)," < 0.0001")),expression(paste(italic(P)," < 0.0001")),expression(paste(italic(P)," = 0.05"))))+
  annotate("segment", x = c(1.05,1.05,2.05,2.05,3.05,1.05), xend = c(1.95,2.95,2.95,3.95,3.95,3.95), y = c(110,120,110,130,110,140), yend = c(110,120,110,130,110,140))+
  theme_classic()+
  theme(axis.title = element_text(size=35,colour="black"),
        axis.text.y = element_text(size=30,colour = "black"),
        axis.text.x = element_text(size=20,colour = "black",angle = 0),
        legend.position = "none"
  )


#### End ####

#### Fig.1e discovery virome heatmap ####
meta<-read.delim(paste0(meta.pa,"meta.data.txt"),colClasses = c(rep("character",12),rep("numeric",4),"character",rep("numeric",2),rep("character",4)),stringsAsFactors = T)
# DNA contig
pa="input/wet.input/discovery.virome/my.analysis/cross.assembly"
name<-list.files(pa, full.names = T,pattern = "D")
l<-c()
l<-lapply(name,read_cross)
data<-do.call(rbind,l)

l<-c()
j<-c()
m<-c()
i=1
name<-list.files(pa, full.names = F,pattern = "D")
for (x in name) {
  tmp<-data%>%filter(grepl(paste0(x,".contig"),contig.id))%>%
    dplyr::select("contig.id","length","circularity","orf.number","uniprot.total.orf.number","uniprot.species","uniprot.species.orf.number","uniprot.family.orf.number","uniprot.family","len.family","final.family",contains(x),contains("DB"),contains("DD"),contains("DT"))%>%
    mutate(max1=pmax(!!as.name(paste0(x,"0")),!!as.name(paste0(x,"1")),!!as.name(paste0(x,"4"))),max2=pmax(DB1,DB2,DB3,DB4,DB5,DD1,DD2,DDH,DDN1,DDN2,DT1,DT2,DT3,DT4,DT5))%>%
    mutate(source=ifelse(max1<=max2,"Contamination",ifelse(uniprot.total.orf.number/orf.number > 0.5,"Viruses","Nonviruses")))%>%
    left_join(host[,c("family","host.name")],by=c("final.family" = "family"))%>%
    left_join(contig.nt)%>%
    mutate(coverage=alignment.len/qlen)
  tmp$host.name<-as.character(tmp$host.name)
  tmp[is.na(tmp$host.name),"host.name"]<-"Others"
  tmp$len.family<-as.character(tmp$len.family)
  tmp[is.na(tmp$len.family),"len.family"]<-"Others"
  tmp$superkingdom<-as.character(tmp$superkingdom)
  tmp[is.na(tmp$superkingdom),"superkingdom"]<-"Others"
  tmp[is.na(tmp$pident),"pident"]<-0
  tmp[is.na(tmp$coverage),"coverage"]<-0
  
  tmp<-tmp%>%mutate(new.source=ifelse(source=="Contamination","Contamination",
                                      ifelse(host.name=="Bacteriophages",source,ifelse(pident>80&coverage>0.8,superkingdom,source))))%>%
    mutate(new.source=ifelse(new.source=="Others",source,new.source))
  
  
  assign(paste0(x,".data"),tmp)
  
  read<-tmp%>%
    dplyr::select(contains(x),new.source)%>%
    gather(key=library_id,value=rpm,-new.source)%>%
    group_by(library_id,new.source)%>%
    dplyr::summarise(total=sum(rpm)/10000)
  l[[i]]<-read
  
  rich<-tmp%>%filter(new.source=="Viruses")%>%
    dplyr::select(uniprot.species,contains(x),contains("DB"),contains("DD"),contains("DT"))%>%
    gather(key=library_id,value=rpm,-uniprot.species)%>%
    group_by(uniprot.species,library_id)%>%
    dplyr::summarise(total=sum(rpm))%>%
    filter(total>10)%>%
    group_by(library_id)%>%
    dplyr::summarise(n=n())
  j[[i]]<-rich
  
  family<-tmp%>%filter(new.source=="Viruses")%>%
    dplyr::select(final.family,contains(x),contains("DB"),contains("DD"),contains("DT"))%>%
    gather(key=library_id,value=rpm,-final.family)%>%
    group_by(final.family,library_id)%>%
    dplyr::summarise(total=sum(rpm))
  
  m[[i]]<-family
  i=i+1
}

family<-do.call(rbind,m)
tmp<-family%>%left_join(meta)



tmp<-tmp%>%dplyr::select(final.family,library_id,total,Study_group)%>%
  arrange(desc(total))%>%
  distinct(final.family,library_id,.keep_all = T)%>%
  spread(key=final.family,value=total)%>%
  replace(is.na(.), 0)%>%
  gather(key=final.family,value=total,-library_id,-Study_group)


tmp1<-left_join(tmp,host[,c("family","host.name")],by=c("final.family" = "family"))
tmp1<-tmp1%>%filter(host.name%in%c("Animal DNA Viruses","Bacteriophages"))
#tmp$host.name<-factor(tmp$host.name,levels=c("Animal DNA Viruses","Animal RNA Viruses","Plant RNA Viruses","Bacteriophages"))


# RNA contig
pa="input/wet.input/discovery.virome/my.analysis/cross.assembly"
rna=c("R345","R352","R357","R359","R365","R368","R375","R385","R396","R400","R401","R428")
name<-paste0(pa,"/",rna)

l<-c()
l<-lapply(name,read_cross)
data<-do.call(rbind,l)

# Based on family orf length
l<-c()
j<-c()
m<-c()
i=1
name<-rna
for (x in name) {
  tmp<-data%>%filter(grepl(paste0(x,".contig"),contig.id))%>%
    dplyr::select("contig.id","length","circularity","orf.number","uniprot.total.orf.number","uniprot.species","uniprot.species.orf.number","uniprot.family.orf.number","uniprot.family","len.family","final.family",contains(x),contains("RB"),contains("RD"),contains("RT"))%>%
    mutate(max1=pmax(!!as.name(paste0(x,"0")),!!as.name(paste0(x,"1")),!!as.name(paste0(x,"4"))),max2=pmax(RB1,RB2,RB3,RB4,RB5,RD1,RD2,RDH,RDN1,RDN2,RT1,RT2,RT3,RT4,RT5))%>%
    mutate(source=ifelse(max1<=max2,"Contamination",ifelse(uniprot.total.orf.number/orf.number > 0.5,"Viruses","Nonviruses")))%>%
    left_join(host[,c("family","host.name")],by=c("len.family" = "family"))%>%
    left_join(contig.nt)%>%
    mutate(coverage=alignment.len/qlen)
  tmp$host.name<-as.character(tmp$host.name)
  tmp[is.na(tmp$host.name),"host.name"]<-"Others"
  tmp$len.family<-as.character(tmp$len.family)
  tmp[is.na(tmp$len.family),"len.family"]<-"Others"
  tmp$superkingdom<-as.character(tmp$superkingdom)
  tmp[is.na(tmp$superkingdom),"superkingdom"]<-"Others"
  tmp[is.na(tmp$pident),"pident"]<-0
  tmp[is.na(tmp$coverage),"coverage"]<-0
  
  tmp<-tmp%>%mutate(new.source=ifelse(source=="Contamination","Contamination",
                                      ifelse(host.name=="Bacteriophages",source,ifelse(pident>80&coverage>0.8,superkingdom,source))))%>%
    mutate(new.source=ifelse(new.source=="Others",source,new.source))
  
  
  assign(paste0(x,".data"),tmp)
  
  read<-tmp%>%
    dplyr::select(contains(x),new.source)%>%
    gather(key=library_id,value=rpm,-new.source)%>%
    group_by(library_id,new.source)%>%
    dplyr::summarise(total=sum(rpm)/10000)
  l[[i]]<-read
  
  
  rich<-tmp%>%filter(new.source=="Viruses")%>%
    dplyr::select(uniprot.species,contains(x),contains("RB"),contains("RD"),contains("RT"))%>%
    gather(key=library_id,value=rpm,-uniprot.species)%>%
    group_by(uniprot.species,library_id)%>%
    dplyr::summarise(total=sum(rpm))%>%
    filter(total>10)%>%
    group_by(library_id)%>%
    dplyr::summarise(n=n())
  j[[i]]<-rich
  
  family<-tmp%>%filter(new.source=="Viruses")%>%
    dplyr::select(final.family,contains(x),contains("RB"),contains("RD"),contains("RT"))%>%
    gather(key=library_id,value=rpm,-final.family)%>%
    group_by(final.family,library_id)%>%
    dplyr::summarise(total=sum(rpm))
  
  m[[i]]<-family
  i=i+1
}


family<-do.call(rbind,m)
tmp<-family%>%left_join(meta)



tmp<-tmp%>%dplyr::select(final.family,library_id,total,Study_group)%>%
  arrange(desc(total))%>%
  distinct(final.family,library_id,.keep_all = T)%>%
  spread(key=final.family,value=total)%>%
  replace(is.na(.), 0)%>%
  gather(key=final.family,value=total,-library_id,-Study_group)


tmp2<-left_join(tmp,host[,c("family","host.name")],by=c("final.family" = "family"))
tmp2<-tmp2%>%filter(host.name%in%c("Animal RNA Viruses"))
#tmp$host.name<-factor(tmp$host.name,levels=c("Animal DNA Viruses","Animal RNA Viruses","Plant RNA Viruses","Bacteriophages"))
tmp2$library_id<-gsub(tmp2$library_id,pattern = "R",replacement = "D")


tmp<-rbind(tmp1[,-5],tmp2[,-5])
tmp<-tmp%>%
  spread(key=final.family,value=total)%>%
  replace(is.na(.), 0)%>%
  gather(key=final.family,value=total,-library_id,-Study_group)
tmp<-left_join(tmp,host[,c("family","host.name")],by=c("final.family" = "family"))
tmp<-tmp%>%group_by(final.family)%>%filter(max(total)>10)


g1<-ggplot(tmp%>%filter(Study_group!="Negative Control"), aes(x = library_id, y = final.family, fill = log(total+1,10))) +
  geom_tile(color="grey", size=0.1) +
  labs(x="",y="")+
  theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        axis.ticks =element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=15,colour="black"),
        strip.background = element_rect(linetype = 0),
        strip.text.x  = element_text(size=20,colour="black"),
        strip.text.y  = element_text(size=15,colour="black",angle = 360),
        legend.position="bottom",
        legend.title = element_text(size=15,margin = margin(t = -20, unit = "pt")),
        legend.title.align = 0.5,
        legend.text = element_text(size=15),
        legend.key=element_rect(color='black'),
        legend.key.height = unit(0.1, "inch"),
        legend.key.width=unit(0.5, "inch")
  )+
  facet_grid(host.name~Study_group,scales = "free", space = "free")+
  scale_fill_continuous(name="Reads per million total reads",na.value="white",guide=guide_colorbar(title.position = "top"),
                        limits = c(log(11,10),log(max(tmp$total)+1,10)),
                        low=brewer.pal(9, name = 'OrRd')[2],high=brewer.pal(9, name = 'OrRd')[8],
                        breaks = c(log(11,10),log(1001,10),log(100001,10)),
                        labels=c(expression("10"),expression("10"^3),expression("10"^5)))+
  
  theme(aspect.ratio = 1)

g <- ggplot_gtable(ggplot_build(g1))
stript <- which(grepl('strip-t', g$layout$name))
fills <- guan.color[1:3]
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

pdf("output/Fig.1e.pdf", width = 20, height = 10)
grid.draw(g)
dev.off()

#### End ####

#### Fig.3a,b animal cell virus Infection ####
#binom.test(1, 79, 0.5,alternative="two.sided",conf.level=0.95)
tmp<-read.delim(paste0(data.pa,"infection.ratio.all.txt"))
tmp<-tmp%>%gather(key=virus,value=rr,-technique,-cohort,-group,-Value)
tmp$virus<-factor(tmp$virus,levels=rev(unique(as.character(tmp$virus))))
tmp<-tmp%>%spread(key=Value,value=rr)
tmp$cohort<-gsub(tmp$cohort,pattern="Urban US 1",replacement = "Discovery cohort\nUS Urban\nn = 20")
tmp$cohort<-gsub(tmp$cohort,pattern="Urban US 2",replacement = "Validation cohort\nUS Urban\nn = 125")
tmp$cohort<-gsub(tmp$cohort,pattern="Africa",replacement = "Validation cohort\nAfrica Rural\nn = 100")
tmp$cohort<-factor(tmp$cohort,levels=c("Discovery cohort\nUS Urban\nn = 20","Validation cohort\nUS Urban\nn = 125","Validation cohort\nAfrica Rural\nn = 100"))
tmp$group<-gsub(tmp$group,pattern = "Breastmilk",replacement = "Breastmilk + Mixed")
tmp$group<-factor(tmp$group,levels=c("Formula","Breastmilk + Mixed"))
g<-ggplot(tmp%>%filter(technique=="Seq",!virus%in%c("size")),aes(x=virus,y=Ratio*100,fill=group))+
  geom_errorbar(aes(ymin=Ratio*100, ymax=CI*100), width=0.2,size=0.5,
                position=position_dodge(width = 0.5))+
  geom_bar(stat="identity",width=0.4,position=position_dodge(width=0.5),color="black")+
  scale_fill_manual(values =guan.color[c(10,11)])+
  coord_flip()+
  labs(x=c(""),y=c("% of positive subject\n(VLP metagenomic sequencing)"))+
  geom_vline(xintercept = 1.5)+
  facet_grid(.~cohort,scales = "free_x")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0, "lines"),
        axis.title = element_text(size=25,color="black"),
        axis.text.y = element_text(size=25,color="black"),
        axis.text.x = element_text(size=10,color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=25,color="black"),
        legend.position = "none"
  )
pdf("output/Fig.3a.pdf",width = 9,height = 12)
g
dev.off()

g<-ggplot(tmp%>%filter(technique=="qPCR",!virus%in%c("size"),virus!="Astrovirus",virus!="Parvovirus"),aes(x=virus,y=Ratio*100,fill=group))+
  geom_errorbar(aes(ymin=Ratio*100, ymax=CI*100), width=0.2,size=0.5,
                position=position_dodge(width = 0.5))+
  geom_bar(stat="identity", width = 0.4,position=position_dodge(width=0.5),color="black")+
  scale_fill_manual(values =guan.color[c(10,11)])+
  coord_flip()+
  labs(x=c(""),y=c("% of positive subject\n(qPCR)"))+
  geom_vline(xintercept = 1.5)+
  facet_grid(.~cohort)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0, "lines"),
        axis.title = element_text(size=25,color="black"),
        axis.text.y = element_text(size=25,color="black"),
        axis.text.x = element_text(size=10,color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=25,color="black"),
        legend.position = "none"
  )
pdf("output/Fig.3b.pdf",width = 12,height = 8)
g
dev.off()
#### End ####

#### Extended Data Fig.2a,b longitudinal cross-assembly ####
pa="input/wet.input/discovery.virome/my.analysis/cross.assembly"
name<-list.files(pa, full.names = T,pattern = "D")
read_cross<-function(x){
  # filter viral contig
  id<-strsplit(x,"\\/")[[1]][length(strsplit(x,"\\/")[[1]])]
  # circularity
  data<-read.delim(paste(x,"/",id,".contig.3k.cir",sep=""),sep="\t",header=T)
  # total orf number
  tmp<-read.delim(paste(x,"/",id,".orf.number",sep=""),sep="\t",header=F)
  tmp<-separate(data = tmp, col = V1, into = c("contig.id"), sep = "\\_")
  tmp<-data.frame(table(tmp$contig))
  names(tmp)<-c("contig.id","orf.number")
  data<-left_join(data,tmp)
  # uniprot
  uniprot<-read.delim(paste(x,"/",id,".3k.uniprot",sep=""),header=F)
  names(uniprot)<-c("contig.id", "uniprot.hit.id", "pident", "qlen","slen","alignment.len" ,"mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  uniprot<- uniprot  %>% 
    separate(col = contig.id, into = c("contig.id","orfs.id"), sep = "\\_") 
  uniprot<-inner_join(uniprot,uniprottax,by="uniprot.hit.id")
  uniprot<-cbind(uniprot,getTaxonomy(uniprot$taxid,"input/database/ac2tax/accessionTaxa.sql"))
  # uniprot orf number >50% belong to virus
  tmp<-uniprot%>%group_by(contig.id)%>%distinct(contig.id,orfs.id,.keep_all= TRUE)%>%
    dplyr::summarise(uniprot.total.orf.number=n())
  data<-left_join(data,tmp)
  data[is.na(data$uniprot.total.orf.number),"uniprot.total.orf.number"]<-0
  
  tmp<-uniprot%>%distinct(contig.id,orfs.id,.keep_all = T)%>%dplyr::select("contig.id","species")
  tmp<-table(tmp)
  tmp1<-data.frame(colnames(tmp)[apply(tmp,1,which.max)])
  tmp<-uniprot%>%distinct(contig.id,orfs.id,.keep_all = T)%>%dplyr::select("contig.id","species")
  tmp<-data.frame(table(tmp))
  tmp2<-data.frame(tapply(tmp$Freq, tmp$contig.id, max))
  tmp3<-cbind(tmp1,tmp2)
  names(tmp3)<-c("uniprot.species","uniprot.species.orf.number")
  tmp3$contig.id<-rownames(tmp3)
  rownames(tmp3)<-NULL
  data<-left_join(data,tmp3,by="contig.id")
  data$uniprot.species<-as.character(data$uniprot.species)
  data$uniprot.species.orf.number<-as.numeric(data$uniprot.species.orf.number)
  data[is.na(data$uniprot.species.orf.number),"uniprot.species.orf.number"]<-0
  data<-data%>%mutate(uniprot.species=ifelse(uniprot.species.orf.number==0,"Others",uniprot.species))
  
  
  tmp<-uniprot%>%distinct(contig.id,orfs.id,.keep_all = T)%>%dplyr::select("contig.id","family")
  tmp<-table(tmp)
  tmp1<-data.frame(colnames(tmp)[apply(tmp,1,which.max)])
  tmp<-uniprot%>%distinct(contig.id,orfs.id,.keep_all = T)%>%dplyr::select("contig.id","family")
  tmp<-data.frame(table(tmp))
  tmp2<-data.frame(tapply(tmp$Freq, tmp$contig.id, max))
  tmp3<-cbind(tmp1,tmp2)
  names(tmp3)<-c("uniprot.family","uniprot.family.orf.number")
  tmp3$contig.id<-rownames(tmp3)
  rownames(tmp3)<-NULL
  data<-left_join(data,tmp3,by="contig.id")
  data$uniprot.family<-as.character(data$uniprot.family)
  data$uniprot.family.orf.number<-as.numeric(data$uniprot.family.orf.number)
  data[is.na(data$uniprot.family.orf.number),"uniprot.family.orf.number"]<-0
  data<-data%>%mutate(uniprot.family=ifelse(uniprot.family.orf.number<=0,"Others",ifelse(uniprot.family.orf.number/uniprot.total.orf.number<=0.5,"Others",uniprot.family)))
  tmp<-uniprot%>%distinct(contig.id,orfs.id,.keep_all = T)%>%dplyr::select("contig.id","family","alignment.len")
  tmp<-tmp%>%group_by(contig.id,family)%>%
    dplyr::summarise(alignment=sum(alignment.len))%>%
    arrange(desc(alignment))%>%distinct(contig.id,.keep_all = T)%>%
    dplyr::rename(len.family=family)
  data<-data%>%left_join(tmp%>%dplyr::select(contig.id,len.family),by="contig.id")%>%
    mutate(final.family=ifelse(len.family==uniprot.family,as.character(len.family),"Others"))%>%
    mutate(final.family=ifelse(is.na(final.family),"Others",final.family))
  
  #data<-data%>%filter(uniprot.total.orf.number/orf.number > 0.5) #
  # all reads
  read.id<-list.files(x, full.names = T,pattern = ".read.txt")
  bam_read<-function(y){
    p = read.delim(y, header=F,
                   stringsAsFactors = FALSE)
    p<-p%>%mutate(rate=V3/(sum(V3)+sum(V4))*1000000)
    names(p)<-c("contig.id", 
                "contig.len",
                "mapped", 
                "unmapped",
                strsplit(strsplit(y,"\\/")[[1]][length(strsplit(y,"\\/")[[1]])],"\\.")[[1]][2])
    return(p[-nrow(p),c(1,5)])}
  l<-lapply(read.id,bam_read)
  tmp<-Reduce(function(x, y) left_join(x, y),l)
  data<-left_join(data,tmp)
  #tmp1<-data%>%arrange_(paste(paste0(id,"0")))%>%slice((nrow(data)-4):nrow(data))
  #tmp2<-data%>%arrange_(paste(paste0(id,"1")))%>%slice((nrow(data)-4):nrow(data))
  #tmp3<-data%>%arrange_(paste(paste0(id,"4")))%>%slice((nrow(data)-4):nrow(data))
  #data<-rbind(tmp1,tmp2,tmp3)%>%distinct(contig.id,.keep_all = T)
  return(data)
}
l<-c()
l<-lapply(name,read_cross)
data<-do.call(rbind,l)

name<-list.files(pa, full.names = F,pattern = "D")
read_virus<-function(x){
  tmp<-data%>%filter(grepl(paste0(x,".contig"),contig.id))%>%
    mutate(max1=pmax(!!as.name(paste0(x,"0")),!!as.name(paste0(x,"1")),!!as.name(paste0(x,"4"))),max2=pmax(DB1,DB2,DB3,DB4,DB5,DD1,DD2,DDH,DDN1,DDN2,DT1,DT2,DT3,DT4,DT5))%>%
    mutate(source=ifelse(max1<=max2,"Contamination",ifelse(uniprot.total.orf.number/orf.number > 0.5,"Viruses","Nonviruses")))%>%
    left_join(host[,c("family","host.name")],by=c("final.family" = "family"))%>%
    left_join(contig.nt)%>%
    mutate(coverage=alignment.len/qlen)%>%
    mutate(new.source=ifelse(source=="Contamination","Contamination",
                             ifelse(!is.na(host.name),ifelse(host.name=="Bacteriophages",source,
                                                             ifelse(pident>80&coverage>0.8,superkingdom,source)),
                                    ifelse(pident>80&coverage>0.8,superkingdom,source))))
  tmp<-tmp%>%filter(new.source=="Viruses")
  tmp1<-tmp%>%arrange_(paste(paste0(x,"0")))%>%slice((nrow(tmp)-4):nrow(tmp))
  tmp2<-tmp%>%arrange_(paste(paste0(x,"1")))%>%slice((nrow(tmp)-4):nrow(tmp))
  tmp3<-tmp%>%arrange_(paste(paste0(x,"4")))%>%slice((nrow(tmp)-4):nrow(tmp))
  data<-rbind(tmp1,tmp2,tmp3)%>%distinct(contig.id,.keep_all = T)
  return(data)
}
l<-c()
l<-lapply(name,read_virus)
data<-do.call(rbind,l)


tmp<-data%>%separate(col = contig.id,sep = "\\.",into = c("subject.id","contig.id"))%>%
  gather(key=library_id,value=rpm,-subject.id,-contig.id,-length,-circularity,-orf.number,-uniprot.total.orf.number,-uniprot.species,-uniprot.species.orf.number,-uniprot.family,-uniprot.family.orf.number,-len.family)%>%
  mutate(group=ifelse(subject.id==substr(library_id, 1, 4),"Within Subjects","Between Subjects"))%>%left_join(meta)%>%
  filter(Study_group=="Month 0"|Study_group=="Month 1"|Study_group=="Month 4")%>%
  mutate(rpm=as.numeric(rpm))


g<-ggplot(tmp,aes(x=group,y=log(as.numeric(rpm)+1,10)))+
  geom_boxplot(outlier.shape = NA,fill="grey")+
  ylab("Log10 reads\nper million total reads")+
  xlab("")+
  theme_classic()+
  annotate("text", x = c(1.5), y = c(6), size=5, label= c(expression(paste(italic(P)," < 0.0001", sep=""))))+
  annotate("segment", x = c(1.05), xend = c(1.95), y = c(5.8), yend = c(5.8))+
  theme(axis.text = element_text(size=20,color="black"),
        axis.title  = element_text(size=25,color="black")
  )
pdf("output/extended.fig.2b.pdf",width=6,height = 5)
g
dev.off()


wilcox.test((tmp%>%filter(group=="Within Subjects"))$rpm,(tmp%>%filter(group=="Between Subjects"))$rpm)

tmp<-log(data[,12:86]+1,10)
rownames(tmp)<-data[,1]
rowann<-data%>%dplyr::select(contig.id,circularity)
rownames(rowann)<-data[,1]
rowann<-rowann%>%separate(col = contig.id,sep = "\\.",into = c("subject.id","contig.id"))
g<-pheatmap(tmp, cluster_rows = F, cluster_cols = F,
            gaps_col = c(seq(0,60,3)), gaps_row = c(15,27,41,55,70,85,100,114,128,142,157,169,183,198,209,224,239,254,269),
            show_rownames = F,cellwidth =8,cellheight = 2,
            color = colorRampPalette(c("white","yellow","red"))(100),
            annotation_row = rowann[,c(1,3)],border_color = NA
)
pdf("output/extended.fig.2a.pdf",width = 15,height = 12)
g
dev.off()
#### End ####

#### Extended Data Fig.2c-h discovery reads proportion ####
# DNA contig
pa="input/wet.input/discovery.virome/my.analysis/cross.assembly"
name<-list.files(pa, full.names = T,pattern = "D")
l<-c()
l<-lapply(name,read_cross)
data<-do.call(rbind,l)
l<-c()
j<-c()
m<-c()
i=1
name<-list.files(pa, full.names = F,pattern = "D")
for (x in name) {
  tmp<-data%>%filter(grepl(paste0(x,".contig"),contig.id))%>%
    dplyr::select("contig.id","length","circularity","orf.number","uniprot.total.orf.number","uniprot.species","uniprot.species.orf.number","uniprot.family.orf.number","uniprot.family","len.family","final.family",contains(x),contains("DB"),contains("DD"),contains("DT"))%>%
    mutate(max1=pmax(!!as.name(paste0(x,"0")),!!as.name(paste0(x,"1")),!!as.name(paste0(x,"4"))),max2=pmax(DB1,DB2,DB3,DB4,DB5,DD1,DD2,DDH,DDN1,DDN2,DT1,DT2,DT3,DT4,DT5))%>%
    mutate(source=ifelse(max1<=max2,"Contamination",ifelse(uniprot.total.orf.number/orf.number > 0.5,"Viruses","Nonviruses")))%>%
    left_join(host[,c("family","host.name")],by=c("final.family" = "family"))%>%
    left_join(contig.nt)%>%
    mutate(coverage=alignment.len/qlen)
  tmp$host.name<-as.character(tmp$host.name)
  tmp[is.na(tmp$host.name),"host.name"]<-"Others"
  tmp$len.family<-as.character(tmp$len.family)
  tmp[is.na(tmp$len.family),"len.family"]<-"Others"
  tmp$superkingdom<-as.character(tmp$superkingdom)
  tmp[is.na(tmp$superkingdom),"superkingdom"]<-"Others"
  tmp[is.na(tmp$pident),"pident"]<-0
  tmp[is.na(tmp$coverage),"coverage"]<-0
  
  tmp<-tmp%>%mutate(new.source=ifelse(source=="Contamination","Contamination",
                                      ifelse(host.name=="Bacteriophages",source,ifelse(pident>80&coverage>0.8,superkingdom,source))))%>%
    mutate(new.source=ifelse(new.source=="Others",source,new.source))
  
  
  assign(paste0(x,".data"),tmp)
  
  read<-tmp%>%
    dplyr::select(contains(x),new.source)%>%
    gather(key=library_id,value=rpm,-new.source)%>%
    group_by(library_id,new.source)%>%
    dplyr::summarise(total=sum(rpm)/10000)
  l[[i]]<-read
  
  rich<-tmp%>%filter(new.source=="Viruses")%>%
    dplyr::select(uniprot.species,contains(x),contains("DB"),contains("DD"),contains("DT"))%>%
    gather(key=library_id,value=rpm,-uniprot.species)%>%
    group_by(uniprot.species,library_id)%>%
    dplyr::summarise(total=sum(rpm))%>%
    filter(total>10)%>%
    group_by(library_id)%>%
    dplyr::summarise(n=n())
  j[[i]]<-rich
  
  family<-tmp%>%filter(new.source=="Viruses")%>%
    dplyr::select(final.family,contains(x),contains("DB"),contains("DD"),contains("DT"))%>%
    gather(key=library_id,value=rpm,-final.family)%>%
    group_by(final.family,library_id)%>%
    dplyr::summarise(total=sum(rpm))
  
  m[[i]]<-family
  i=i+1
}

tmp<-do.call(rbind,l)
tmp$new.source<-gsub(tmp$new.source,pattern = "Eukaryota",replacement = "Nonviruses")
tmp$new.source<-gsub(tmp$new.source,pattern = "Bacteria",replacement = "Nonviruses")

tmp<-tmp%>%group_by(library_id,new.source)%>%dplyr::summarise(total=sum(total))%>%
  spread(key=new.source,value=total)%>%mutate(Unassigned=100-Viruses-Contamination)%>%
  gather(key=new.source,value=total,-library_id)%>%
  left_join(meta)

g<-ggplot(tmp%>%filter(new.source=="Viruses"),aes(x=Study_group,y=total,fill=Study_group))+
  geom_boxplot()+
  scale_fill_manual(values = guan.color[1:3])+
  ggtitle("DNA virome") +
  xlab("")+
  ylab("Viral read proportion")+
  theme_classic()+
  theme(legend.position = "none",
        plot.title = element_text(color="black", size=20),
        axis.title = element_text(size=20,color="black"),
        axis.text = element_text(size=15,color="black"))
pdf("output/extended.fig2c.pdf",width = 5,height = 5)
g
dev.off()


g<-ggplot(tmp%>%filter(new.source=="Contamination"),aes(x=Study_group,y=total,fill=Study_group))+
  geom_boxplot()+
  scale_fill_manual(values = guan.color[1:3])+
  ggtitle("DNA virome") +
  xlab("")+
  ylab("Contamination read proportion")+
  theme_classic()+
  theme(legend.position = "none",
        plot.title = element_text(color="black", size=20),
        axis.title = element_text(size=20,color="black"),
        axis.text = element_text(size=15,color="black"))
pdf("output/extended.fig2e.pdf",width = 5,height = 5)
g
dev.off()


g<-ggplot(tmp%>%filter(new.source=="Nonviruses"),aes(x=Study_group,y=total,fill=Study_group))+
  geom_boxplot()+
  scale_fill_manual(values = guan.color[1:3])+
  ggtitle("DNA virome") +
  xlab("")+
  ylab("Non-virus microorganism\nread proportion")+
  theme_classic()+
  theme(legend.position = "none",
        plot.title = element_text(color="black", size=20),
        axis.title = element_text(size=20,color="black"),
        axis.text = element_text(size=15,color="black"))
#pdf("output/extended.fig2e.pdf",width = 5,height = 5)
#g
#dev.off()


g<-ggplot(tmp%>%filter(new.source=="Unassigned"),aes(x=Study_group,y=total,fill=Study_group))+
  geom_boxplot()+
  scale_fill_manual(values = guan.color[1:3])+
  ggtitle("DNA virome") +
  xlab("")+
  ylab("Unassigned read proportion")+
  theme_classic()+
  theme(legend.position = "none",
        plot.title = element_text(color="black", size=20),
        axis.title = element_text(size=20,color="black"),
        axis.text = element_text(size=15,color="black"))
pdf("output/extended.fig2d.pdf",width = 5,height = 5)
g
dev.off()

meanandsd<-function(x,y){
  a<-mean((tmp%>%filter(new.source==x,Study_group==y))$total)
  b<-sd((tmp%>%filter(new.source==x,Study_group==y))$total)/sqrt(length((tmp%>%filter(new.source==x,Study_group==y))$total))
  c<-range((tmp%>%filter(new.source==x,Study_group==y))$total)
  l=c(a,b,c)
  return(l)
}

meanandsd("Viruses","Month 0")
meanandsd("Viruses","Month 1")
meanandsd("Viruses","Month 4")
meanandsd("Contamination","Month 0")
meanandsd("Contamination","Month 1")
meanandsd("Contamination","Month 4")
meanandsd("Unassigned","Month 0")
meanandsd("Unassigned","Month 1")
meanandsd("Unassigned","Month 4")





# RNA contig
pa="input/wet.input/discovery.virome/my.analysis/cross.assembly"
rna=c("R345","R352","R357","R359","R365","R368","R375","R385","R396","R400","R401","R428")
name<-paste0(pa,"/",rna)
l<-c()
l<-lapply(name,read_cross)
data<-do.call(rbind,l)
# Based on family orf length
l<-c()
j<-c()
m<-c()
i=1
name<-rna
for (x in name) {
  tmp<-data%>%filter(grepl(paste0(x,".contig"),contig.id))%>%
    dplyr::select("contig.id","length","circularity","orf.number","uniprot.total.orf.number","uniprot.species","uniprot.species.orf.number","uniprot.family.orf.number","uniprot.family","len.family","final.family",contains(x),contains("RB"),contains("RD"),contains("RT"))%>%
    mutate(max1=pmax(!!as.name(paste0(x,"0")),!!as.name(paste0(x,"1")),!!as.name(paste0(x,"4"))),max2=pmax(RB1,RB2,RB3,RB4,RB5,RD1,RD2,RDH,RDN1,RDN2,RT1,RT2,RT3,RT4,RT5))%>%
    mutate(source=ifelse(max1<=max2,"Contamination",ifelse(uniprot.total.orf.number/orf.number > 0.5,"Viruses","Nonviruses")))%>%
    left_join(host[,c("family","host.name")],by=c("len.family" = "family"))%>%
    left_join(contig.nt)%>%
    mutate(coverage=alignment.len/qlen)
  tmp$host.name<-as.character(tmp$host.name)
  tmp[is.na(tmp$host.name),"host.name"]<-"Others"
  tmp$len.family<-as.character(tmp$len.family)
  tmp[is.na(tmp$len.family),"len.family"]<-"Others"
  tmp$superkingdom<-as.character(tmp$superkingdom)
  tmp[is.na(tmp$superkingdom),"superkingdom"]<-"Others"
  tmp[is.na(tmp$pident),"pident"]<-0
  tmp[is.na(tmp$coverage),"coverage"]<-0
  
  tmp<-tmp%>%mutate(new.source=ifelse(source=="Contamination","Contamination",
                                      ifelse(host.name=="Bacteriophages",source,ifelse(pident>80&coverage>0.8,superkingdom,source))))%>%
    mutate(new.source=ifelse(new.source=="Others",source,new.source))
  
  
  assign(paste0(x,".data"),tmp)
  
  read<-tmp%>%
    dplyr::select(contains(x),new.source)%>%
    gather(key=library_id,value=rpm,-new.source)%>%
    group_by(library_id,new.source)%>%
    dplyr::summarise(total=sum(rpm)/10000)
  l[[i]]<-read
  
  
  rich<-tmp%>%filter(new.source=="Viruses")%>%
    dplyr::select(uniprot.species,contains(x),contains("RB"),contains("RD"),contains("RT"))%>%
    gather(key=library_id,value=rpm,-uniprot.species)%>%
    group_by(uniprot.species,library_id)%>%
    dplyr::summarise(total=sum(rpm))%>%
    filter(total>10)%>%
    group_by(library_id)%>%
    dplyr::summarise(n=n())
  j[[i]]<-rich
  
  family<-tmp%>%filter(new.source=="Viruses")%>%
    dplyr::select(final.family,contains(x),contains("RB"),contains("RD"),contains("RT"))%>%
    gather(key=library_id,value=rpm,-final.family)%>%
    group_by(final.family,library_id)%>%
    dplyr::summarise(total=sum(rpm))
  
  m[[i]]<-family
  i=i+1
}
tmp<-do.call(rbind,l)
tmp$new.source<-gsub(tmp$new.source,pattern = "Eukaryota",replacement = "Nonviruses")
tmp$new.source<-gsub(tmp$new.source,pattern = "Bacteria",replacement = "Nonviruses")


tmp<-tmp%>%group_by(library_id,new.source)%>%dplyr::summarise(total=sum(total))%>%
  spread(key=new.source,value=total)%>%mutate(Unassigned=100-Viruses-Contamination)%>%
  gather(key=new.source,value=total,-library_id)%>%
  left_join(meta)


g<-ggplot(tmp%>%filter(new.source=="Viruses"),aes(x=Study_group,y=total,fill=Study_group))+
  geom_boxplot()+
  scale_fill_manual(values = guan.color[1:3])+
  ggtitle("RNA virome") +
  xlab("")+
  ylab("Viral read proportion")+
  theme_classic()+
  theme(legend.position = "none",
        plot.title = element_text(color="black", size=20),
        axis.title = element_text(size=20,color="black"),
        axis.text = element_text(size=15,color="black"))
pdf("output/extended.fig2f.pdf",width = 5,height = 5)
g
dev.off()

g<-ggplot(tmp%>%filter(new.source=="Contamination"),aes(x=Study_group,y=total,fill=Study_group))+
  geom_boxplot()+
  scale_fill_manual(values = guan.color[1:3])+
  ggtitle("RNA virome") +
  xlab("")+
  ylab("Contamination read proportion")+
  theme_classic()+
  theme(legend.position = "none",
        plot.title = element_text(color="black", size=20),
        axis.title = element_text(size=20,color="black"),
        axis.text = element_text(size=15,color="black"))
pdf("output/extended.fig2h.pdf",width = 5,height = 5)
g
dev.off()

# g<-ggplot(tmp%>%filter(new.source!="Contamination",new.source!="Viruses"),aes(x=Study_group,y=total,fill=Study_group))+
#   geom_boxplot()+
#   scale_fill_manual(values = guan.color[1:3])+
#   ggtitle("RNA virome") +
#   xlab("")+
#   ylab("Non-virus microorganism\nread proportion")+
#   theme_classic()+
#   theme(legend.position = "none",
#         plot.title = element_text(color="black", size=20),
#         axis.title = element_text(size=20,color="black"),
#         axis.text = element_text(size=15,color="black"))
# pdf("output/nonviruas.pdf",width = 5,height = 5)
# g
# dev.off()


g<-ggplot(tmp%>%filter(new.source=="Unassigned"),aes(x=Study_group,y=total,fill=Study_group))+
  geom_boxplot()+
  scale_fill_manual(values = guan.color[1:3])+
  ggtitle("RNA virome") +
  xlab("")+
  ylab("Unassigned read proportion")+
  theme_classic()+
  theme(legend.position = "none",
        plot.title = element_text(color="black", size=20),
        axis.title = element_text(size=20,color="black"),
        axis.text = element_text(size=15,color="black"))
pdf("output/extended.fig2g.pdf",width = 5,height = 5)
g
dev.off()

tmp<-tmp%>%filter(!is.na(total))
meanandsd("Viruses","Month 0")
meanandsd("Viruses","Month 1")
meanandsd("Viruses","Month 4")
meanandsd("Contamination","Month 0")
meanandsd("Contamination","Month 1")
meanandsd("Contamination","Month 4")
meanandsd("Unassigned","Month 0")
meanandsd("Unassigned","Month 1")
meanandsd("Unassigned","Month 4")


#### End ####

#### Extended Data Fig.5 Longitudinal cross-assembly contig life cycle prediction ####

pa="input/wet.input/discovery.virome/my.analysis/cross.assembly"
name<-list.files(pa, full.names = T,pattern = "D")
read_cross<-function(x){
  # filter viral contig
  id<-strsplit(x,"\\/")[[1]][length(strsplit(x,"\\/")[[1]])]
  # circularity
  data<-read.delim(paste(x,"/",id,".contig.3k.cir",sep=""),sep="\t",header=T)
  # total orf number
  tmp<-read.delim(paste(x,"/",id,".orf.number",sep=""),sep="\t",header=F)
  tmp<-separate(data = tmp, col = V1, into = c("contig.id"), sep = "\\_")
  tmp<-data.frame(table(tmp$contig))
  names(tmp)<-c("contig.id","orf.number")
  data<-left_join(data,tmp)
  # uniprot
  uniprot<-read.delim(paste(x,"/",id,".3k.uniprot",sep=""),header=F)
  names(uniprot)<-c("contig.id", "uniprot.hit.id", "pident", "qlen","slen","alignment.len" ,"mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  uniprot<- uniprot  %>% 
    separate(col = contig.id, into = c("contig.id","orfs.id"), sep = "\\_") 
  uniprot<-inner_join(uniprot,uniprottax,by="uniprot.hit.id")
  uniprot<-cbind(uniprot,getTaxonomy(uniprot$taxid,"/media/lorax/users/guanxiang/1_igram/liang2019/input/database/ac2tax/accessionTaxa.sql"))
  # uniprot orf number >50% belong to virus
  tmp<-uniprot%>%group_by(contig.id)%>%distinct(contig.id,orfs.id,.keep_all= TRUE)%>%
    dplyr::summarise(uniprot.total.orf.number=n())
  data<-left_join(data,tmp)
  data[is.na(data$uniprot.total.orf.number),"uniprot.total.orf.number"]<-0
  
  tmp<-uniprot%>%distinct(contig.id,orfs.id,.keep_all = T)%>%dplyr::select("contig.id","species")
  tmp<-table(tmp)
  tmp1<-data.frame(colnames(tmp)[apply(tmp,1,which.max)])
  tmp<-uniprot%>%distinct(contig.id,orfs.id,.keep_all = T)%>%dplyr::select("contig.id","species")
  tmp<-data.frame(table(tmp))
  tmp2<-data.frame(tapply(tmp$Freq, tmp$contig.id, max))
  tmp3<-cbind(tmp1,tmp2)
  names(tmp3)<-c("uniprot.species","uniprot.species.orf.number")
  tmp3$contig.id<-rownames(tmp3)
  rownames(tmp3)<-NULL
  data<-left_join(data,tmp3,by="contig.id")
  data$uniprot.species<-as.character(data$uniprot.species)
  data$uniprot.species.orf.number<-as.numeric(data$uniprot.species.orf.number)
  data[is.na(data$uniprot.species.orf.number),"uniprot.species.orf.number"]<-0
  data<-data%>%mutate(uniprot.species=ifelse(uniprot.species.orf.number==0,"Others",uniprot.species))
  
  
  tmp<-uniprot%>%distinct(contig.id,orfs.id,.keep_all = T)%>%dplyr::select("contig.id","family")
  tmp<-table(tmp)
  tmp1<-data.frame(colnames(tmp)[apply(tmp,1,which.max)])
  tmp<-uniprot%>%distinct(contig.id,orfs.id,.keep_all = T)%>%dplyr::select("contig.id","family")
  tmp<-data.frame(table(tmp))
  tmp2<-data.frame(tapply(tmp$Freq, tmp$contig.id, max))
  tmp3<-cbind(tmp1,tmp2)
  names(tmp3)<-c("uniprot.family","uniprot.family.orf.number")
  tmp3$contig.id<-rownames(tmp3)
  rownames(tmp3)<-NULL
  data<-left_join(data,tmp3,by="contig.id")
  data$uniprot.family<-as.character(data$uniprot.family)
  data$uniprot.family.orf.number<-as.numeric(data$uniprot.family.orf.number)
  data[is.na(data$uniprot.family.orf.number),"uniprot.family.orf.number"]<-0
  data<-data%>%mutate(uniprot.family=ifelse(uniprot.family.orf.number<=0,"Others",ifelse(uniprot.family.orf.number/uniprot.total.orf.number<=0.5,"Others",uniprot.family)))
  tmp<-uniprot%>%distinct(contig.id,orfs.id,.keep_all = T)%>%dplyr::select("contig.id","family","alignment.len")
  tmp<-tmp%>%group_by(contig.id,family)%>%
    dplyr::summarise(alignment=sum(alignment.len))%>%
    arrange(desc(alignment))%>%distinct(contig.id,.keep_all = T)%>%
    dplyr::rename(len.family=family)
  data<-data%>%left_join(tmp%>%dplyr::select(contig.id,len.family),by="contig.id")%>%
    mutate(final.family=ifelse(len.family==uniprot.family,as.character(len.family),"Others"))%>%
    mutate(final.family=ifelse(is.na(final.family),"Others",final.family))
  
  #data<-data%>%filter(uniprot.total.orf.number/orf.number > 0.5) #
  # all reads
  read.id<-list.files(x, full.names = T,pattern = ".read.txt")
  bam_read<-function(y){
    p = read.delim(y, header=F,
                   stringsAsFactors = FALSE)
    p<-p%>%mutate(rate=V3/(sum(V3)+sum(V4))*1000000)
    names(p)<-c("contig.id", 
                "contig.len",
                "mapped",
                "unmapped",
                strsplit(strsplit(y,"\\/")[[1]][length(strsplit(y,"\\/")[[1]])],"\\.")[[1]][2])
    return(p[-nrow(p),c(1,5)])}
  l<-lapply(read.id,bam_read)
  tmp<-Reduce(function(x, y) left_join(x, y),l)
  data<-left_join(data,tmp)
  #tmp1<-data%>%arrange_(paste(paste0(id,"0")))%>%slice((nrow(data)-4):nrow(data))
  #tmp2<-data%>%arrange_(paste(paste0(id,"1")))%>%slice((nrow(data)-4):nrow(data))
  #tmp3<-data%>%arrange_(paste(paste0(id,"4")))%>%slice((nrow(data)-4):nrow(data))
  #data<-rbind(tmp1,tmp2,tmp3)%>%distinct(contig.id,.keep_all = T)
  return(data)
}
l<-c()
l<-lapply(name,read_cross)
data<-do.call(rbind,l)

name<-list.files(pa, full.names = F,pattern = "D")
read_virus<-function(x){
  tmp<-data%>%filter(grepl(paste0(x,".contig"),contig.id))%>%
    mutate(max1=pmax(!!as.name(paste0(x,"0")),!!as.name(paste0(x,"1")),!!as.name(paste0(x,"4"))),max2=pmax(DB1,DB2,DB3,DB4,DB5,DD1,DD2,DDH,DDN1,DDN2,DT1,DT2,DT3,DT4,DT5))%>%
    mutate(source=ifelse(max1<=max2,"Contamination",ifelse(uniprot.total.orf.number/orf.number > 0.5,"Viruses","Nonviruses")))%>%
    left_join(host[,c("family","host.name")],by=c("final.family" = "family"))%>%
    left_join(contig.nt)%>%
    mutate(coverage=alignment.len/qlen)%>%
    mutate(new.source=ifelse(source=="Contamination","Contamination",
                             ifelse(!is.na(host.name),ifelse(host.name=="Bacteriophages",source,
                                                             ifelse(pident>80&coverage>0.8,superkingdom,source)),
                                    ifelse(pident>80&coverage>0.8,superkingdom,source))))
  tmp<-tmp%>%filter(new.source=="Viruses")
  return(tmp)
}
l<-c()
l<-lapply(name,read_virus)
viral.contig.all<-do.call(rbind,l)


viral.contig.10orf<-viral.contig.all%>%filter(orf.number>=10,host.name=="Bacteriophages")
#write.table(viral.contig.10orf[,1:2],"/media/lorax/users/guanxiang/1_igram/1_virome/analysis/1_vlp_analysis/my.analysis/lifecycle/viral.contig.txt",quote = F,row.names = F)
#bash
#tail -n +2 viral.contig.txt | awk -F" " '{print $1}' |while read p; do awk "/^>/{x = /${p}_/;}(x)" orf.fa > ${p}.fa; perl ~/software/PHACTS/phacts.pl -f ${p}.fa --classes ~/software/PHACTS/classes_lifestyle -r 100 > ${p}.lifestyle; remove ${p}.fa;done

name<-list.files("input/wet.input/discovery.virome/my.analysis/life",full.names = T,pattern = "*style")
l<-c()
j=1
for (i in name){
  tmp<-read.delim(i,skip=3,header=F)[,c("V1","V2")]%>%spread(key=V1,value=V2)
  contig.id<-paste(unlist(strsplit(unlist(strsplit(i, "\\/"))[6],"\\."))[1],unlist(strsplit(unlist(strsplit(i, "\\/"))[6],"\\."))[2],sep=".")
  tmp$contig.id=contig.id
  l[[j]]=tmp
  j=j+1
}
data.life<-do.call("rbind", l)%>%mutate(life=Temperate-0.5)
name<-list.files(pa, full.names = F,pattern = "D")
l<-c()
j=1
for (x in name) {
  tmp<-viral.contig.10orf%>%left_join(data.life)%>%
    #mutate(life=ifelse(is.na(coverage),life,ifelse(pident>0.8&coverage>0.8&superkingdom=="Bacteria",max(life),life)))%>%
    #mutate(life=ifelse(is.na(life),0,life))%>%
    filter(grepl(paste0(x,".contig"),contig.id))%>%
    dplyr::select(life,contig.id,contains(x),contains("DD"),contains("DB"),contains("DT"))%>%
    gather(key=library_id,value=rpm,-life,-contig.id)
  l[[j]]=tmp
  j=j+1
}

data<-do.call(rbind, l)%>%left_join(meta)%>%dplyr::select(Study_group,life,rpm,contig.id)



sum((data%>%filter(life<0,Study_group=="Month 0"))$rpm)/(sum((data%>%filter(life>0,Study_group=="Month 0"))$rpm)+sum((data%>%filter(life<0,Study_group=="Month 0"))$rpm))
sum((data%>%filter(life<0,Study_group=="Month 1"))$rpm)/(sum((data%>%filter(life>0,Study_group=="Month 1"))$rpm)+sum((data%>%filter(life<0,Study_group=="Month 1"))$rpm))
sum((data%>%filter(life<0,Study_group=="Month 4"))$rpm)/(sum((data%>%filter(life>0,Study_group=="Month 4"))$rpm)+sum((data%>%filter(life<0,Study_group=="Month 4"))$rpm))


data%>%filter(life<0,Study_group=="Month 0",rpm>1)%>%distinct(contig.id)%>%nrow()
data%>%filter(life<0,Study_group=="Month 1",rpm>1)%>%distinct(contig.id)%>%nrow()
data%>%filter(life<0,Study_group=="Month 4",rpm>1)%>%distinct(contig.id)%>%nrow()



pdf("output/extended.data.fig5.pdf",width=10,height=5)
ggplot(data%>%group_by(contig.id)%>%distinct(contig.id,.keep_all = T), aes(x=as.numeric(life))) +
  #geom_point()+
  #geom_density(fill=guan.color[5])+
  geom_histogram(color="white", position="identity",binwidth=0.01,boundary = -0.05,fill=guan.color[5])+
  #geom_freqpoly( position="identity",binwidth=0.005,boundary = -0.05)+
  #geom_density(aes(y=..count..))+
  # geom_freqpoly(data=tmp2,aes(x=life))+
  #scale_fill_manual(values = guan.color[5])+
  coord_cartesian(xlim = c(-0.15, 0.15)) +
  #geom_rect(aes(xmin = -0.8, xmax = 0, ymin = -Inf, ymax = Inf),
  #          fill = "pink", alpha = 0.01)+
  labs(x="",y="Viral contig number")+
  geom_vline(xintercept = 0)+
  scale_x_continuous(breaks=c(-0.15,-0.075,0,0.075,0.15),
                     labels=c("1","0.75","0.5","0.75","1"))+
  #facet_grid(Study_group ~ .)+
  #  annotate("text", x=0.4, y=0.1, label= paste("P = ",,sep="")+
  theme_classic()+
  theme(
    axis.title = element_text(size=20,colour="black"),
    axis.text = element_text(size=15,colour="black"),
    strip.background = element_blank(),
    strip.text=element_blank(),
    legend.position="none"
  )

# grid.text(percent(r1), x = unit(0.45, "npc"), y = unit(0.90, "npc"),gp=gpar(fontsize=20, col="black"))
# grid.text(percent(1-r1), x = unit(0.62, "npc"), y = unit(0.90, "npc"),gp=gpar(fontsize=20, col="black"))
# grid.text(percent(r2), x = unit(0.45, "npc"), y = unit(0.60, "npc"),gp=gpar(fontsize=20, col="black"))
# grid.text(percent(1-r2), x = unit(0.62, "npc"), y = unit(0.60, "npc"),gp=gpar(fontsize=20, col="black"))
# grid.text(percent(r3), x = unit(0.45, "npc"), y = unit(0.3, "npc"),gp=gpar(fontsize=20, col="black"))
# grid.text(percent(1-r3), x = unit(0.62, "npc"), y = unit(0.3, "npc"),gp=gpar(fontsize=20, col="black"))
# grid.text("Month 0", x = unit(0.18, "npc"), y = unit(0.95, "npc"),gp=gpar(fontsize=20, col="black"))
# grid.text("Month 1", x = unit(0.18, "npc"), y = unit(0.65, "npc"),gp=gpar(fontsize=20, col="black"))
# grid.text("Month 4", x = unit(0.18, "npc"), y = unit(0.35, "npc"),gp=gpar(fontsize=20, col="black"))
grid.text("Lytic", x = unit(0.1, "npc"), y = unit(0.04, "npc"),gp=gpar(fontsize=20, col="black"))
grid.text("Temperate", x = unit(0.92, "npc"), y = unit(0.04, "npc"),gp=gpar(fontsize=20, col="black"))
grid.text("Probability", x = unit(0.53, "npc"), y = unit(0.04, "npc"),gp=gpar(fontsize=20, col="black"))


grid.lines(x = unit(c(0.25, 0.45), "npc"),
           y = unit(c(0.04, 0.04), "npc"),
           gp = gpar(fill="black"),
           arrow = arrow(length = unit(0.2, "inches"), angle = 5,
                         ends="first", type="closed"))
grid.lines(x = unit(c(0.6, 0.8), "npc"),
           y = unit(c(0.04, 0.04), "npc"),
           gp = gpar(fill="black"),
           arrow = arrow(length = unit(0.2, "inches"), angle = 5,
                         ends="last", type="closed"))
dev.off()



#### End ####

#### Extended Data Fig.8 Animal cell virus coverage and prevalence ####

meta<-read.delim(paste0(meta.pa,"meta.data.txt"),colClasses = c(rep("character",12),rep("numeric",4),"character",rep("numeric",2),rep("character",4)),stringsAsFactors = T)
meta$Infant_feeding_type<-gsub(meta$Infant_feeding_type,pattern = "Mixed",replacement = "Breastmilk")
meta$Infant_feeding_type<-gsub(meta$Infant_feeding_type,pattern = "Breastmilk",replacement = "Breastmilk + Mixed")
meta$Infant_feeding_type<-factor(meta$Infant_feeding_type,levels=c("Formula","Breastmilk + Mixed"))
meta$Infant_delivery_type<-gsub(meta$Infant_delivery_type,pattern = "C-Section without labor",replacement ="C-Section")
meta$Infant_delivery_type<-gsub(meta$Infant_delivery_type,pattern = "C-Section with labor",replacement ="C-Section")

# Fuctions  
read_bwa_coverage <- function(filename1,filename2) {
  p<-try(read.delim(filename1, header=F,stringsAsFactors = FALSE))
  if (!inherits(p, 'try-error')) {
    names(p)<-c("family", 
                "read",
                "align.len","genome.len", "cov")
    p<-p[,c(1,2,5)]
    x<- p %>% filter(read>0)
    if (nrow(x) == 0) {
      x<-p[1,c("family","cov")]
      x$cov=0
    } else {
      x<- aggregate(.~family,x[,c("family","cov")],sum)
    }
    y = read.delim(filename2, header=F,
                   stringsAsFactors = FALSE)
    names(y)<-c("family", 
                "genome.len",
                "mapped","unmapped")
    y<-y %>% mutate (total.read=sum(mapped)+sum(unmapped))
    y<-y[-nrow(y),c(1,2,3,5)]
    z<-merge(x,y,by="family")
    z<-z %>% mutate (rpkm=mapped/(genome.len/1000*total.read/1000000))
    #z<-z %>% mutate (rpm=(mapped/total.read)*1000000)
    z<- z %>% arrange(desc(cov))
    z<-z[1,c(1,2,4,5)]
    z$sample_id<-strsplit(strsplit(filename1,"\\/")[[1]][12],"\\.")[[1]][1]
    names(z)<-c("virus_id","cov","reads","total.read","library_id")
    return(z)
  } else {
    #p=data.frame(family=NA,read=0,align.len=0,genome.len=0,cov=0)
    #x<- p %>% filter(read>0)
    y = read.delim(filename2, header=F,
                   stringsAsFactors = FALSE)
    names(y)<-c("family", 
                "genome.len",
                "mapped","unmapped")
    y<-y %>% mutate (total.read=sum(mapped)+sum(unmapped))
    y<-y[,c(1,2,3,5)]
    y$cov<-0
    z<-y[1,c(1,5,3,4)]
    z$sample_id<-strsplit(strsplit(filename1,"\\/")[[1]][12],"\\.")[[1]][1]
    names(z)<-c("virus_id","cov","reads","total.read","library_id")
    return(z)
  }
  #n<-data.frame(virus_id=NA,cov=0,reads=0,total.read=0,sample_id=strsplit(strsplit(filename1,"\\/")[[1]][10],"\\.")[[1]][1])
  
}
table_bwa_coverage<- function(viralname){
  path1 = paste("/media/lorax/users/guanxiang/1_igram/liang2019/input/wet.input/coverage","/",viralname,sep="") 
  filename.cov <- list.files(path = path1,full.names = TRUE,pattern = ".coverage")
  filename.reads <- list.files(path = path1,full.names = TRUE,pattern = ".read")
  bwa.cov.data <- mapply(read_bwa_coverage,filename.cov,filename.reads,SIMPLIFY = FALSE) # mapply do not give list default as lapply, need simplify option
  bwa.cov<-do.call(rbind,bwa.cov.data)
  rownames(bwa.cov)<-NULL
  return(bwa.cov)
}

# Each family
adeno<-table_bwa_coverage("adeno")
adeno$family<-"adeno"
anello<-table_bwa_coverage("anello")
anello$family<-"anello"
astro<-table_bwa_coverage("astro")
astro$family<-"astro"
calici<-table_bwa_coverage("calici")
calici$family<-"calici"
parvo<-table_bwa_coverage("parvo")
parvo$family<-"parvo"
picona<-table_bwa_coverage("picona")
picona$family<-"picona"
polymavi<-table_bwa_coverage("polymavi")
polymavi$family<-"polymavi"
data<-rbind(adeno,anello,astro,calici,parvo,picona,polymavi)
data$library_id<-gsub(data$library_id,pattern = "R",replacement = "D")
data<-left_join(data,meta)%>%mutate(RPM=reads/total.read*1000000)

# Discovery cohort
df<-data%>%filter(Study_group=="Month 4",Cohort=="Discovery")%>%arrange(desc(cov))
df<-df[!duplicated(df$library_id),]
l<-c()
j=1
for (i in seq(0,1,by=0.001)){
  x=nrow(df %>%filter(Infant_feeding_type=="Breastmilk + Mixed",cov>=get("i")))/nrow(df %>%filter(Infant_feeding_type=="Breastmilk + Mixed"))
  y=nrow(df %>%filter(Infant_feeding_type=="Formula",cov>=get("i")))/nrow(df %>%filter(Infant_feeding_type=="Formula"))
  l[[j]]<-data.frame(cov=i,Breastmilk=x,Formula=y)
  j=j+1
}
cov<-do.call(rbind,l)
mcov<-melt(cov,id.vars = "cov") 
mcov$cov<-as.numeric(mcov$cov)
mcov$variable<-gsub(mcov$variable,pattern = "Breastmilk",replacement = "Breastmilk + Mixed")
mcov$variable<-factor(mcov$variable,levels=c("Formula","Breastmilk + Mixed"))

g<-ggplot(mcov,aes(x=cov,y=value*100,color=variable))+
  geom_line(size=2)+
  coord_cartesian(ylim = c(0, 100))+
  scale_color_manual(values=guan.color[10:11])+
  ggtitle("Discovery cohort",subtitle = "Feeding type")+
  xlab("Genome coverage cutoff")+
  ylab("Percentage of\npositive subjects (%)")+
  theme_classic()+
  theme(legend.title = element_blank(),
        legend.position = c(0.6,0.9),
        legend.text = element_text(size=15,color="black"),
        plot.title = element_text(size=20,color="black"),
        plot.subtitle = element_text(size=15,color="black"),
        axis.text = element_text(size=10,color="black"),
        axis.title = element_text(size=20,color="black"))
pdf("output/extended.fig8a.pdf",width = 5,height = 5)
g
dev.off()
#  pvalue
df<-data%>%filter(Study_group=="Month 4",Cohort=="Discovery")%>%arrange(desc(cov))
df<-df[!duplicated(df$library_id),]
l<-c()
j=1
for (i in seq(0,1,by=0.001)){
  df$Infection<-"Uninfected"
  df[df$cov>=i,]$Infection<-"Infected"
  fisher=table(df%>%dplyr::select(Infant_feeding_type,Infection))
  if (ncol(fisher)!=1){
    p<-fisher.test(fisher)[1][[1]]} else{
      p=1
    }
  l[[j]]<-data.frame(cov=i,p=p)
  j=j+1
}
p<-do.call(rbind,l)

g<-ggplot(p,aes(x=cov,y=p))+
  geom_line()+
  geom_hline(yintercept = 0.05,color="red")+
  coord_cartesian(ylim = c(0, 1))+
  ggtitle("Discovery cohort",subtitle = "Feeding type")+
  xlab("Genome coverage cutoff")+
  ylab("P value")+
  theme_classic()+
  theme(legend.title = element_blank(),
        legend.position = c(0.6,0.9),
        legend.text = element_text(size=15,color="black"),
        plot.title = element_text(size=20,color="black"),
        plot.subtitle = element_text(size=15,color="black"),
        axis.text = element_text(size=10,color="black"),
        axis.title = element_text(size=20,color="black"))
pdf("output/extended.fig8b.pdf",width = 5,height = 5)
g
dev.off()

# Validation cohort
df<-data%>%filter(Study_group=="Month 4"|Study_group=="Month 3",Cohort=="Validation")%>%arrange(desc(cov))
df<-df[!duplicated(df$library_id),]
l<-c()
j=1
for (i in seq(0,1,by=0.001)){
  x=nrow(df %>%filter(Infant_feeding_type=="Breastmilk + Mixed",cov>=get("i")))/nrow(df %>%filter(Infant_feeding_type=="Breastmilk + Mixed"))
  y=nrow(df %>%filter(Infant_feeding_type=="Formula",cov>=get("i")))/nrow(df %>%filter(Infant_feeding_type=="Formula"))
  l[[j]]<-data.frame(cov=i,Breastmilk=x,Formula=y)
  j=j+1
}
cov<-do.call(rbind,l)
mcov<-melt(cov,id.vars = "cov") 
mcov$cov<-as.numeric(mcov$cov)
mcov$variable<-gsub(mcov$variable,pattern = "Breastmilk",replacement = "Breastmilk + Mixed")
mcov$variable<-factor(mcov$variable,levels=c("Formula","Breastmilk + Mixed"))

g<-ggplot(mcov,aes(x=cov,y=value*100,color=variable))+
  geom_line(size=2)+
  coord_cartesian(ylim = c(0, 100))+
  scale_color_manual(values=guan.color[10:11])+
  ggtitle("Validation cohort",subtitle = "Feeding type")+
  xlab("Genome coverage cutoff")+
  ylab("Percentage of\npositive subjects (%)")+
  theme_classic()+
  theme(legend.title = element_blank(),
        legend.position = c(0.6,0.9),
        legend.text = element_text(size=15,color="black"),
        plot.title = element_text(size=20,color="black"),
        plot.subtitle = element_text(size=15,color="black"),
        axis.text = element_text(size=10,color="black"),
        axis.title = element_text(size=20,color="black"))
pdf("output/extended.fig8c.pdf",width = 5,height = 5)
g
dev.off()

df<-data%>%filter(Study_group=="Month 4"|Study_group=="Month 3",Cohort=="Validation")%>%arrange(desc(cov))
df<-df[!duplicated(df$library_id),]
l<-c()
j=1
for (i in seq(0,1,by=0.001)){
  df$Infection<-"Uninfected"
  df[df$cov>=i,]$Infection<-"Infected"
  fisher=table(df%>%dplyr::select(Infant_feeding_type,Infection))
  if (ncol(fisher)!=1){
    p<-fisher.test(fisher)[1][[1]]} else{
      p=1
    }
  l[[j]]<-data.frame(cov=i,p=p)
  j=j+1
}
p<-do.call(rbind,l)

g<-ggplot(p,aes(x=cov,y=p))+
  geom_line()+
  geom_hline(yintercept = 0.05,color="red")+
  coord_cartesian(ylim = c(0, 1))+
  ggtitle("Validation cohort",subtitle = "Feeding type")+
  xlab("Genome coverage cutoff")+
  ylab("P value")+
  theme_classic()+
  theme(legend.title = element_blank(),
        legend.position = c(0.6,0.9),
        legend.text = element_text(size=15,color="black"),
        plot.title = element_text(size=20,color="black"),
        plot.subtitle = element_text(size=15,color="black"),
        axis.text = element_text(size=10,color="black"),
        axis.title = element_text(size=20,color="black"))
pdf("output/extended.fig8d.pdf",width = 5,height = 5)
g
dev.off()

# negative control
df<-data%>%filter(Study_group=="Negative Control")%>%arrange(desc(cov))%>%dplyr::distinct(library_id,.keep_all = T)
df$library_id<-factor(df$library_id,levels=unique(as.character(df$library_id)))
g<-ggplot(df,aes(x=library_id,y=cov))+
  geom_bar(stat="identity")+
  coord_cartesian(ylim = c(0, 0.2))+
  theme_classic()+
  ylab("Genome coverage")+
  xlab("Negative control samples")+
  theme(legend.title = element_blank(),
        legend.position = c(0.6,0.9),
        legend.text = element_text(size=15,color="black"),
        plot.title = element_text(size=20,color="black"),
        plot.subtitle = element_text(size=15,color="black"),
        axis.text.y = element_text(size=10,color="black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=20,color="black"))
pdf("output/extended.fig8k.pdf",width = 5,height = 5)
g
dev.off()

df<-data%>%filter(Study_group=="Month 4"|Study_group=="Month 3",Cohort=="Validation")
wilcox.test(RPM~Infant_feeding_type,data=df)
g<-ggplot(df,aes(x=Infant_feeding_type,y=log(RPM+1,2),fill=Infant_feeding_type))+
  geom_boxplot()+
  coord_cartesian(ylim = c(0, 5))+
  annotate("text", x = c(1.5), y = c(4.5), size=7, label= c(expression(paste(italic(P)," < 0.0001 ", sep=""))))+
  annotate("segment", x = c(1.05), xend = c(1.95), y = c(4), yend = c(4))+
  ggtitle("Validation cohort",subtitle = "Feeding type")+
  xlab("")+
  ylab("Log2(reads per million total reads)")+
  scale_fill_manual(values = guan.color[10:11])+
  theme_classic()+
  theme(axis.title = element_text(size=20,color = "black"),
        axis.text.y = element_text(size=10,color = "black"),
        axis.text.x = element_text(size=20,color = "black"),
        plot.title = element_text(size=20,color="black"),
        plot.subtitle = element_text(size=15,color="black"),
        legend.position = "none")
pdf("output/extended.fig8e.pdf",width = 5,height = 5)
g
dev.off()

## Delivery type
# Discovery cohort
df<-data%>%filter(Study_group=="Month 4",Cohort=="Discovery")%>%arrange(desc(cov))
df<-df[!duplicated(df$library_id),]
l<-c()
j=1
for (i in seq(0,1,by=0.001)){
  x=nrow(df %>%filter(Infant_delivery_type=="C-Section",cov>=get("i")))/nrow(df %>%filter(Infant_delivery_type=="C-Section"))
  y=nrow(df %>%filter(Infant_delivery_type=="Spontaneous vaginal delivery",cov>=get("i")))/nrow(df %>%filter(Infant_delivery_type=="Spontaneous vaginal delivery"))
  l[[j]]<-data.frame(cov=i,CSection=x,SVD=y)
  j=j+1
}
cov<-do.call(rbind,l)
mcov<-melt(cov,id.vars = "cov") 
mcov$cov<-as.numeric(mcov$cov)
mcov$variable<-gsub(mcov$variable,pattern = "CSection",replacement = "C-Section")
mcov$variable<-gsub(mcov$variable,pattern = "SVD",replacement = "Spontaneous vaginal delivery")


g<-ggplot(mcov,aes(x=cov,y=value*100,color=variable))+
  geom_line(size=2)+
  coord_cartesian(ylim = c(0, 100))+
  scale_color_manual(values=guan.color[10:11])+
  ggtitle("Discovery cohort",subtitle = "Delivery type")+
  xlab("Genome coverage cutoff")+
  ylab("Percentage of\npositive subjects (%)")+
  theme_classic()+
  theme(legend.title = element_blank(),
        legend.position = c(0.6,0.9),
        legend.text = element_text(size=15,color="black"),
        plot.title = element_text(size=20,color="black"),
        plot.subtitle = element_text(size=15,color="black"),
        axis.text = element_text(size=10,color="black"),
        axis.title = element_text(size=20,color="black"))
pdf("output/extended.fig8f.pdf",width = 5,height = 5)
g
dev.off()
#  pvalue
df<-data%>%filter(Study_group=="Month 4",Cohort=="Discovery")%>%arrange(desc(cov))
df<-df[!duplicated(df$library_id),]
l<-c()
j=1
for (i in seq(0,1,by=0.001)){
  df$Infection<-"Uninfected"
  df[df$cov>=i,]$Infection<-"Infected"
  fisher=table(df%>%dplyr::select(Infant_delivery_type,Infection))
  if (ncol(fisher)!=1){
    p<-fisher.test(fisher)[1][[1]]} else{
      p=1
    }
  l[[j]]<-data.frame(cov=i,p=p)
  j=j+1
}
p<-do.call(rbind,l)

g<-ggplot(p,aes(x=cov,y=p))+
  geom_line()+
  geom_hline(yintercept = 0.05,color="red")+
  coord_cartesian(ylim = c(0, 1))+
  ggtitle("Discovery cohort",subtitle = "Delivery type")+
  xlab("Genome coverage cutoff")+
  ylab("P value")+
  theme_classic()+
  theme(legend.title = element_blank(),
        legend.position = c(0.6,0.9),
        legend.text = element_text(size=15,color="black"),
        plot.title = element_text(size=20,color="black"),
        plot.subtitle = element_text(size=15,color="black"),
        axis.text = element_text(size=10,color="black"),
        axis.title = element_text(size=20,color="black"))
pdf("output/extended.fig8g.pdf",width = 5,height = 5)
g
dev.off()

# Validation cohort
df<-data%>%filter(Study_group=="Month 4"|Study_group=="Month 3",Cohort=="Validation")%>%arrange(desc(cov))
df<-df[!duplicated(df$library_id),]
l<-c()
j=1
for (i in seq(0,1,by=0.001)){
  x=nrow(df %>%filter(Infant_delivery_type=="C-Section",cov>=get("i")))/nrow(df %>%filter(Infant_delivery_type=="C-Section"))
  y=nrow(df %>%filter(Infant_delivery_type=="Spontaneous vaginal delivery",cov>=get("i")))/nrow(df %>%filter(Infant_delivery_type=="Spontaneous vaginal delivery"))
  l[[j]]<-data.frame(cov=i,CSection=x,SVD=y)
  j=j+1
}
cov<-do.call(rbind,l)
mcov<-melt(cov,id.vars = "cov") 
mcov$cov<-as.numeric(mcov$cov)
mcov$variable<-gsub(mcov$variable,pattern = "CSection",replacement = "C-Section")
mcov$variable<-gsub(mcov$variable,pattern = "SVD",replacement = "Spontaneous vaginal delivery")


g<-ggplot(mcov,aes(x=cov,y=value*100,color=variable))+
  geom_line(size=2)+
  coord_cartesian(ylim = c(0, 100))+
  scale_color_manual(values=guan.color[10:11])+
  ggtitle("Validation cohort",subtitle = "Delivery type")+
  xlab("Genome coverage cutoff")+
  ylab("Percentage of\npositive subjects (%)")+
  theme_classic()+
  theme(legend.title = element_blank(),
        legend.position = c(0.6,0.9),
        legend.text = element_text(size=15,color="black"),
        plot.title = element_text(size=20,color="black"),
        plot.subtitle = element_text(size=15,color="black"),
        axis.text = element_text(size=10,color="black"),
        axis.title = element_text(size=20,color="black"))
pdf("output/extended.fig8h.pdf",width = 5,height = 5)
g
dev.off()

df<-data%>%filter(Study_group=="Month 4"|Study_group=="Month 3",Cohort=="Validation")%>%arrange(desc(cov))
df<-df[!duplicated(df$library_id),]
l<-c()
j=1
for (i in seq(0,1,by=0.001)){
  df$Infection<-"Uninfected"
  df[df$cov>=i,]$Infection<-"Infected"
  fisher=table(df%>%dplyr::select(Infant_delivery_type,Infection))
  if (ncol(fisher)!=1){
    p<-fisher.test(fisher)[1][[1]]} else{
      p=1
    }
  l[[j]]<-data.frame(cov=i,p=p)
  j=j+1
}
p<-do.call(rbind,l)

g<-ggplot(p,aes(x=cov,y=p))+
  geom_line()+
  geom_hline(yintercept = 0.05,color="red")+
  coord_cartesian(ylim = c(0, 1))+
  ggtitle("Validation cohort",subtitle = "Delivery type")+
  xlab("Genome coverage cutoff")+
  ylab("P value")+
  theme_classic()+
  theme(legend.title = element_blank(),
        legend.position = c(0.6,0.9),
        legend.text = element_text(size=15,color="black"),
        plot.title = element_text(size=20,color="black"),
        plot.subtitle = element_text(size=15,color="black"),
        axis.text = element_text(size=10,color="black"),
        axis.title = element_text(size=20,color="black"))
pdf("output/extended.fig8i.pdf",width = 5,height = 5)
g
dev.off()



df<-data%>%filter(Study_group=="Month 4"|Study_group=="Month 3",Cohort=="Validation")
wilcox.test(RPM~Infant_delivery_type,data=df)
g<-ggplot(df,aes(x=Infant_delivery_type,y=log(RPM+1,2),fill=Infant_delivery_type))+
  geom_boxplot()+
  coord_cartesian(ylim = c(0, 5))+
  scale_x_discrete(breaks=c("C-Section","Spontaneous vaginal delivery"),
                   labels=c("C-Section", "Spontaneous\nvaginal delivery"))+
  annotate("text", x = c(1.5), y = c(4.5), size=5, label= c(expression(paste(italic(P)," = 0.11 ", sep=""))))+
  annotate("segment", x = c(1.05), xend = c(1.95), y = c(4), yend = c(4))+
  ggtitle("Validation cohort",subtitle = "Delivery type")+
  xlab("")+
  ylab("Log2(reads per million total reads)")+
  scale_fill_manual(values = guan.color[10:11])+
  theme_classic()+
  theme(axis.title = element_text(size=20,color = "black"),
        axis.text.y = element_text(size=10,color = "black"),
        axis.text.x = element_text(size=20,color = "black"),
        plot.title = element_text(size=20,color="black"),
        plot.subtitle = element_text(size=15,color="black"),
        legend.position = "none")
pdf("output/extended.fig8j.pdf",width = 5,height = 5)
g
dev.off()

#### End ####

#### Extended Data Fig.7 crAssphage ####
all.crass<-table_bwa_coverage("crassphage")

data<-all.crass%>%left_join(meta)
data<-data%>%arrange(desc(cov))
data<-data[!duplicated(data$library_id),]
data$Infection<-"Uninfected"
data[data$cov>=(1/3),"Infection"]<-"Infected"
data$Infection<-factor(data$Infection,levels=c("Infected","Uninfected"))
data<-data%>%filter(library_type=="DNA")



df<-data%>%filter(Study_group=="Month 0"|Study_group=="Month 1"|Study_group=="Month 4"|Study_group=="Year 2~5"|Study_group=="Month 3")%>%dplyr::select(Study_group,Infection)%>%group_by(Study_group,Infection)%>%
  dplyr::summarise(n=n())%>%
  mutate(freq = n/sum(n))%>%
  ungroup() %>%
  complete(Study_group, Infection,
           fill = list(n = 0, freq = 0))


g<-ggplot(df%>%filter(Infection=="Infected")) +
  geom_bar(aes(y = freq*100, x = Study_group),fill="grey", color="black",stat="identity",width = 0.5,size=0.5)+
  expand_limits(y = c(0, 50))+
  labs(title="",x ="", y = "Percentage of subjects (%)\nwith crAssphage")+
  theme_classic()+
  theme(axis.title = element_text(size=25,colour="black"),
        axis.text.x = element_text(size=20,colour = "black"),
        axis.text.y = element_text(size=20,colour = "black",margin = margin(t = -20, unit = "pt")),
        strip.text = element_text(size=25),
        legend.position = "none"
  )

pdf("output/extended.fig7.pdf", width = 10, height = 6)
g
dev.off()
#### End ####

#### Extended Data Fig.9 Pfam phage pupulation ####
pfam.nc<-read.delim("input/database/pfamA_ncbi.txt",header=F)
names(pfam.nc)<-c("acc","un","tax.id")

# Discovery
path<-"input/wet.input/discovery.virome/my.analysis/orfs/"
name<-list.files(path,pattern = "^D")
l1=c()
l.tmp=c()
i=1
for (x in name) {
  y<-read.delim(paste(path,x,"/",x,".hmm.out",sep=""))
  y<-y %>% 
    filter(evalue<1e-5) %>%
    separate(col=acc,into = c("acc","acc.number"),sep = "\\.")
  y<-y[!duplicated(y$query),]
  pfam.ta<-read.delim("input/database/pfamA_tax_depth.txt",header=F)
  names(pfam.ta)<-c("acc","taxa","V3","V4","V5")
  a<-merge(y[,c("acc","query")],pfam.ta[,c("acc","taxa")]) %>% filter(taxa=="Viruses")%>%distinct(acc)%>%nrow()
  b<-merge(y[,c("acc","query")],pfam.ta[,c("acc","taxa")]) %>% filter(taxa=="Bacteria")%>%distinct(acc)%>%nrow()
  c<-merge(y[,c("acc","query")],pfam.ta[,c("acc","taxa")]) %>% filter(taxa=="Eukaryota")%>%distinct(acc)%>%nrow()
  d<-merge(y[,c("acc","query")],pfam.ta[,c("acc","taxa")]) %>% filter(taxa=="Archaea")%>%distinct(acc)%>%nrow()
  e<-merge(y[,c("acc","query")],pfam.ta[,c("acc","taxa")]) %>% filter(taxa=="Viruses"|taxa=="Bacteria")%>%distinct(acc,taxa)%>%group_by(acc)%>%filter(n()>1)%>%nrow()
  f<-merge(y[,c("acc","query")],pfam.ta[,c("acc","taxa")]) %>% filter(taxa=="Viruses"|taxa=="Eukaryota")%>%distinct(acc,taxa)%>%group_by(acc)%>%filter(n()>1)%>%nrow()
  g<-merge(y[,c("acc","query")],pfam.ta[,c("acc","taxa")]) %>% filter(taxa=="Viruses"|taxa=="Archaea")%>%distinct(acc,taxa)%>%group_by(acc)%>%filter(n()>1)%>%nrow()
  u<-merge(y[,c("acc","query")],pfam.ta[,c("acc","taxa")])%>%distinct(acc,taxa)%>%group_by(acc)%>%filter(n()==1,taxa=="Viruses")%>%nrow
  t<-merge(y[,c("acc","query")],pfam.ta[,c("acc","taxa")])%>%distinct(acc)%>%nrow
  h<-data.frame(virus=a,bacteria=b,eukaryota=c,archaea=d,vir.bac=e/2,vir.euk=f/2,vir.arc=g/2,vir.u=u,total=t,library_id=x)
  y<-merge(y[,c("acc","query")],pfam.ta[,c("acc","taxa")]) %>% filter(taxa=="Viruses")
  z<-read.delim(paste(path,x,"/",x,".orf.count",sep=""),skip = 1)
  names(z)[1]<-"query"
  names(z)[7]<-"read"
  y<-merge(y,z,all.x=T,by="query")
  z<-read.delim(paste(path,x,"/",x,".orf.count.summary",sep=""),skip = 1,header=F)
  y$total<-z[1,2]
  y<-y%>%mutate(ratio=read/total)
  if (nrow(y) == 0) {
    z<-y
  } else {
    y<-aggregate(.~acc,y[,c("acc","ratio")],sum)
    names(y)[2]<-x
    z<-y
  }
  
  l1[[i]]<-z
  l.tmp[[i]]<-h
  i=i+1
  
}

#Validation 
path<-"input/wet.input/validation.virome/my.analysis/contigs/"
name<-list.files(path,pattern = "^D")
l2=c()
i=1
for (x in name) {
  y<-read.delim(paste(path,x,"/",x,".hmm.out",sep=""))
  y<-y %>% 
    filter(evalue<1e-5) %>%
    separate(col=acc,into = c("acc","acc.number"),sep = "\\.")
  y<-y[!duplicated(y$query),]
  pfam.ta<-read.delim("input/database/pfamA_tax_depth.txt",header=F)
  names(pfam.ta)<-c("acc","taxa")
  y<-merge(y[,c("acc","query")],pfam.ta[,c("acc","taxa")]) %>% filter(taxa=="Viruses")
  z<-read.delim(paste(path,x,"/",x,".orf.count",sep=""),skip = 1)
  names(z)[1]<-"query"
  names(z)[7]<-"read"
  y<-merge(y,z,all.x=T,by="query")
  z<-read.delim(paste(path,x,"/",x,".orf.count.summary",sep=""),skip = 1,header=F)
  y$total<-z[1,2]
  y<-y%>%mutate(ratio=read/total)
  if (nrow(y) == 0) {
    z<-y
  } else {
    y<-aggregate(.~acc,y[,c("acc","ratio")],sum)
    names(y)[2]<-x
    z<-y
  }
  
  l2[[i]]<-y
  i=i+1
  
}


#Discovery and Validation together 
host <- read.delim("input/meta.data/viral.host.txt") # virus host data
com<-Reduce(function(x, y) merge(x, y, all=TRUE,by=c("acc")), c(l1,l2)) 
l.pfam<-left_join(com,pfam.nc[,c("acc","tax.id")],by="acc")
l.pfam<-cbind(l.pfam,getTaxonomy(l.pfam$tax.id,"input/database/ac2tax/accessionTaxa.sql")[,c(1,5)])
l.pfam<-l.pfam%>%filter(superkingdom=="Viruses")
l.pfam<-merge(l.pfam,host[,1:2],by="family")
l.pfam<-l.pfam%>%filter(host=="Bacteria")
l.pfam<-l.pfam[!duplicated(l.pfam$acc),]

meta<-read.delim(paste0(meta.pa,"meta.data.txt"),colClasses = c(rep("character",12),rep("numeric",4),"character",rep("numeric",2),rep("character",4)),stringsAsFactors = T)
meta$Infant_feeding_type<-gsub(meta$Infant_feeding_type,pattern = "Mixed",replacement ="Breastmilk")
meta$Infant_feeding_type<-gsub(meta$Infant_feeding_type,pattern = "Breastmilk",replacement ="Breastmilk + Mixed")
meta$Infant_delivery_type<-gsub(meta$Infant_delivery_type,pattern = "C-Section with labor",replacement = "C-Section")
meta$Infant_delivery_type<-gsub(meta$Infant_delivery_type,pattern = "C-Section without labor",replacement = "C-Section")
meta$library_type<-factor(meta$library_type,levels=c("DNA","RNA"))
#meta$Study_group<-gsub(meta$Study_group,pattern = "Month 3",replacement = "Month 4")
meta$Study_group<-factor(meta$Study_group,levels=c("Month 0","Month 1","Month 3","Month 4","Year 2~5","Negative Control","Germ-free Mice","Positive Control"))
meta$Infant_feeding_type<-factor(meta$Infant_feeding_type,levels=c("Formula","Breastmilk + Mixed"))

# pfam profile
com<-l.pfam[,-c(1,ncol(l.pfam),ncol(l.pfam)-1,ncol(l.pfam)-2)]
com[is.na(com)]<-0



# Pfam phage Age effect 
subdata<-meta%>%filter(Study_group=="Month 0"|Study_group=="Month 1"|Study_group=="Month 3"|Study_group=="Month 4"|Study_group=="Negative Control",library_type=="DNA",Experiment=="Virome")
subdata$Study_group<-droplevels(subdata$Study_group)
keep<-as.character(subdata$library_id)
tcom<-com[,names(com) %in% keep]
tcom<-tcom[rowSums(tcom[,-1]>0)>=1,]
subdata<-subdata%>%filter(library_id %in%as.character(names(tcom)))
subdata$library_id<-factor(subdata$library_id,levels = as.character(subdata$library_id))
tcom<-tcom[,match(as.character(subdata$library_id),colnames(tcom))]
tmp<-as.matrix(vegdist(t(tcom), method="bray"))


#test
adonis(tmp ~ Study_group, data = subdata)

# To plot
compc<-pcoa(tmp)
#biplot(compc,t(com[,names(com) %in% keep]))
compc_df<- merge(meta,compc$vectors[,1:2],by.x="library_id",by.y="row.names")
compc_pct <- round(compc$values$Relative_eig * 100)
g = ggplot(compc_df, aes(x=Axis.1, y=Axis.2))+
  geom_point(aes(shape = Study_group,color= Study_group),size = 3) +
  #stat_ellipse(aes(color= Study_group))+
  xlab(paste0("Axis 1 (", compc_pct[1], "%)")) +
  ylab(paste0("Axis 2 (", compc_pct[2], "%)")) +
  scale_shape_manual(values=c(15:18,8:12))+
  scale_color_manual(values=c(guan.color[1:2],guan.color[4],guan.color[c(3,5)]))+
  annotate("text", x=0.5, y=0.4,size=8, label= "FDR = 0.007")+
  coord_equal() +
  theme_classic()+
  theme(axis.title = element_text(size=25,colour="black"),
        axis.text = element_text(size=20,colour = "black"),
        legend.position = c(0.8,0.15),
        legend.title = element_blank(),
        legend.text = element_text(size=20,colour = "black"))
pdf("output/extended.fig9b.pdf",width=7,height=7,useDingbats = F)
plot(g)
dev.off()


# Month 3,4 samples without negative control
subdata<-meta%>%filter(Study_group=="Month 4",library_type=="DNA")
subdata$Study_group<-droplevels(subdata$Study_group)
keep<-intersect(as.character(subdata$library_id),names(com))
tcom<-com[,names(com) %in% keep]
tcom<-tcom[rowSums(tcom[,-1]>0)>=1,]
subdata<-subdata%>%filter(library_id %in% keep)
subdata$library_id<-factor(subdata$library_id,levels = as.character(subdata$library_id))
tcom<-tcom[,match(as.character(subdata$library_id),colnames(tcom))]
tmp<-as.matrix(vegdist(t(tcom), method="bray"))

adonis(tmp ~ Mother_pregnancy_induce_HTN.diabetes, data = subdata)
adonis(tmp ~ Mother_Chorioamnionitis, data = subdata)


subdata<-meta%>%filter(Study_group=="Month 4"|Study_group=="Month 3",library_type=="DNA")
subdata$Study_group<-droplevels(subdata$Study_group)
keep<-intersect(as.character(subdata$library_id),names(com))
tcom<-com[,names(com) %in% keep]
tcom<-tcom[rowSums(tcom[,-1]>0)>=1,]
subdata<-subdata%>%filter(library_id %in% keep)
subdata$library_id<-factor(subdata$library_id,levels = as.character(subdata$library_id))
tcom<-tcom[,match(as.character(subdata$library_id),colnames(tcom))]
tmp<-as.matrix(vegdist(t(tcom), method="bray"))


#test
adonis(tmp ~ Infant_feeding_type, data = subdata)
adonis(tmp ~ Infant_delivery_type, data = subdata)
adonis(tmp ~ Infant_gender, data = subdata)
adonis(tmp ~ Mother_body_type, data = subdata)
#adonis(tmp ~ Infant_formula_type, data = subdata)# p= 0.41


# plot
compc<-pcoa(tmp)
compc_df<- merge(meta,compc$vectors[,1:2],by.x="library_id",by.y="row.names")
compc_pct <- round(compc$values$Relative_eig * 100)
g = ggplot(compc_df, aes(x=Axis.1, y=Axis.2))+
  geom_point(aes(shape = Infant_feeding_type,fill= Infant_feeding_type),size = 3,shape=21,color="black") +
  stat_ellipse(aes(color= Infant_feeding_type))+
  xlab(paste0("Axis 1 (", compc_pct[1], "%)")) +
  ylab(paste0("Axis 2 (", compc_pct[2], "%)")) +
  #geom_text(aes(label=library_id),hjust=0, vjust=0)
  #scale_shape_manual(values=c(15:18,8:12))+
  scale_fill_manual(values=guan.color[c(10,11)])+
  scale_color_manual(values=guan.color[c(10,11)])+
  annotate("text", x=0.4, y=0.4,size=8, label="FDR = 0.007")+
  #coord_equal() +
  theme_classic()+
  theme(axis.title = element_text(size=25,colour="black"),
        axis.text = element_text(size=20,colour = "black"),
        legend.position = c(0.25,0.9),
        legend.title = element_blank(),
        legend.text = element_text(size=20,colour = "black"))

pdf("output/extended.fig9c.pdf",width=7,height=7,useDingbats = F)
plot(g)
dev.off()


g = ggplot(compc_df, aes(x=Axis.1, y=Axis.2))+
  geom_point(aes(shape = Infant_gender,fill= Infant_gender),size = 3,shape=21,color="black") +
  stat_ellipse(aes(color= Infant_gender))+
  #geom_text(aes(label=library_id),hjust=0, vjust=0)+
  xlab(paste0("Axis 1 (", compc_pct[1], "%)")) +
  ylab(paste0("Axis 2 (", compc_pct[2], "%)")) +
  #geom_text(aes(label=library_id),hjust=0, vjust=0)
  #scale_shape_manual(values=c(15:18,8:12))+
  scale_fill_manual(values=guan.color[c(10,11)])+
  scale_color_manual(values=guan.color[c(10,11)])+
  annotate("text", x=0.4, y=0.4,size=8, label= "FDR = 0.23")+
  #coord_equal() +
  theme_classic()+
  theme(axis.title = element_text(size=25,colour="black"),
        axis.text = element_text(size=20,colour = "black"),
        legend.position = c(0.4,0.9),
        legend.title = element_blank(),
        legend.text = element_text(size=20,colour = "black"))
pdf("output/extended.fig9e.pdf",width=7,height=7,useDingbats = F)
plot(g)
dev.off()

# g = ggplot(compc_df, aes(x=Axis.1, y=Axis.2))+
#   geom_point(aes(shape = Mother_body_type,fill= Mother_body_type),size = 3,shape=21,color="black") +
#   stat_ellipse(aes(color=Mother_body_type))+
#   #geom_text(aes(label=library_id),hjust=0, vjust=0)+
#   xlab(paste0("Axis 1 (", compc_pct[1], "%)")) +
#   ylab(paste0("Axis 2 (", compc_pct[2], "%)")) +
#   #geom_text(aes(label=library_id),hjust=0, vjust=0)
#   #scale_shape_manual(values=c(15:18,8:12))+
#   scale_fill_manual(values=guan.color[c(10,11)])+
#   scale_color_manual(values=guan.color[c(10,11)])+
#   annotate("text", x=0.4, y=0.4,size=8, label= "FDR = 0.56")+
#   #coord_equal() +
#   theme_classic()+
#   theme(axis.title = element_text(size=25,colour="black"),
#         axis.text = element_text(size=20,colour = "black"),
#         legend.position = c(0.4,0.9),
#         legend.title = element_blank(),
#         legend.text = element_text(size=20,colour = "black"))
# pdf("output/extended.fig9e.pdf",width=7,height=7,useDingbats = F)
# plot(g)
# dev.off()

g = ggplot(compc_df, aes(x=Axis.1, y=Axis.2))+
  geom_point(aes(shape = Infant_delivery_type,fill= Infant_delivery_type),size = 3,shape=21,color="black") +
  stat_ellipse(aes(color= Infant_delivery_type))+
  #geom_text(aes(label=library_id),hjust=0, vjust=0)+
  xlab(paste0("Axis 1 (", compc_pct[1], "%)")) +
  ylab(paste0("Axis 2 (", compc_pct[2], "%)")) +
  #geom_text(aes(label=library_id),hjust=0, vjust=0)
  #scale_shape_manual(values=c(15:18,8:12))+
  scale_fill_manual(values=guan.color[c(10,11)])+
  scale_color_manual(values=guan.color[c(10,11)])+
  annotate("text", x=0.39, y=0.4,size=8, label="FDR = 0.93")+
  #coord_equal() +
  theme_classic()+
  theme(axis.title = element_text(size=25,colour="black"),
        axis.text = element_text(size=20,colour = "black"),
        legend.position = c(0.375,0.9),
        legend.title = element_blank(),
        legend.text = element_text(size=20,colour = "black"))
pdf("output/extended.fig9d.pdf",width=7,height=7,useDingbats = F)
plot(g)
dev.off()

#### End ####

# rm(list=ls())
# options(warn=0)










