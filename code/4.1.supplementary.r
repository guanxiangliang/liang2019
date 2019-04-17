#### Set up working directory ####
data.pa<-"./input/wet.input/"
meta.pa<-"./input/meta.data/"
database.pa<-"./input/database/"
options(warn=0)
#### End ####

#### Libraries ####
library(ggplot2)
library(dplyr)
library(plyr) #rbind.fill function
library(reshape)
library(ggsci)
library(viridis)
library(viridisLite)
library(viridis)
library(tidyr) # separate function
library(gtools) # mixedsort
library(scales) # for trans_form in ggplot
library(grid)
library(RColorBrewer) # brew.col
library(taxonomizr)### convert accession id to tax id and to names ###
# getNamesAndNodes()
# getAccession2taxid()
# getAccession2taxid(types='prot')
# read.names.sql('names.dmp','./input/database/ac2tax/accessionTaxa.sql')
# read.nodes.sql('nodes.dmp','./input/database/ac2tax/accessionTaxa.sql')
# read.accession2taxid(list.files('.','accession2taxid.gz$'),'./input/database/ac2tax/accessionTaxa.sql')
# file.remove(list.files('.','accession2taxid.gz$'))
# file.remove(list.files('.','dmp$'))
# #### End ####

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

#### Fig. S3 Lytic phage prediction ####
source("./code/dark.R")
# 
# # make viral.contig.annotation and write out contig's D3400.contig.ann
names<-list.files(paste(data.pa,"/discovery.virome/my.analysis/contigs/",sep=""), pattern="^D")
for (i in names){
  tryCatch({
    tmp<-dark(i)
    assign(paste0(i,".contig.ann"),tmp) #paste0 can be used for objective with variables
    #     write.table(tmp[,"contig"],paste("rtable/",i,".viral.contig.txt",sep=""),quote = FALSE,row.names = F,  col.names = FALSE )
  }, error=function(e){})
}
# 
# ## move to shell to extract orf fasta and predict lifestyle
# #cat D3400.viral.contig.txt | while read p; do awk "/^>/{x = /${p}_/;}(x)" /media/lorax/users/guanxiang/1_igram/1_igram_virome_hiseq/analysis/1_vlp_analysis/my.analysis/orfs/D3400/D3400.contig.3k.orf.fa > ${p}.fa; perl ~/software/PHACTS/phacts.pl -f ${p}.fa --classes ~/software/PHACTS/classes_lifestyle -r 100 > ${p}.lifestyle; remove ${p}.fa;done


samples<-list.files(paste(data.pa,"discovery.virome/my.analysis/life",sep=""))
for (j in samples) {
  tryCatch({
    names<-list.files(paste(data.pa,"discovery.virome/my.analysis/life/",j,"/",sep=""),pattern=paste("*style"))
    l<-c()
    for (i in names){
      tryCatch({
        contig.id<-paste(unlist(strsplit(i, "\\."))[1],unlist(strsplit(i, "\\."))[2],sep=".")
        tmp<-read.delim(paste(data.pa,"discovery.virome/my.analysis/life/",j,"/",i,sep=""),skip=3,header=F)[,c("V1","V2")]
        tmp<-as.data.frame(t(tmp))
        names(tmp)<-as.vector(t(tmp[1,]))
        
        rownames(tmp)<-c("style",contig.id)  
        tmp<-tmp[-1,]
        l[[contig.id]]<-tmp
      }, error=function(e){})
    }
    tmp1<-do.call("rbind", l)
    tmp1$contig<-rownames(tmp1)
    tmp2<-get(paste(j,".contig.ann",sep=""))
    tmp<-merge(tmp2,tmp1,by="contig")
    assign(paste0(j,".contig.ann.life"),tmp)
  }, error=function(e){})
}


samples<-list.files(paste(data.pa,"discovery.virome/my.analysis/life",sep=""),pattern = "0$") # $ means 0 is the end of the string
l<-c()
for (i in samples){ 
  tryCatch({
  tmp<-get(paste(i,".contig.ann.life",sep="")) %>% filter(source=="Viruses",uniprot.total.orf.number>=10) 
    tmp$Temperate<-as.numeric(as.character(tmp$Temperate))
  tmp$Lytic<-as.numeric(as.character(tmp$Lytic))
  tmp<-mutate(tmp,
              life=Temperate-0.5)
  tmp<-cbind(tmp,data.frame(group=rep("Month 0",nrow(tmp))))
  l[[i]]<-tmp[,c("contig","contig.length","life","group","vlp.reads.number")]
  }, error=function(e){})
}
tmp1<-do.call("rbind", l)


samples<-list.files(paste(data.pa,"discovery.virome/my.analysis/life",sep=""),pattern = "1$")
l<-c()
for (i in samples){ 
  tmp<-get(paste(i,".contig.ann.life",sep="")) %>% filter(source=="Viruses",uniprot.total.orf.number>=10) 
  tmp$Temperate<-as.numeric(as.character(tmp$Temperate))
  tmp$Lytic<-as.numeric(as.character(tmp$Lytic))
  tmp<-mutate(tmp,
              life=Temperate-0.5)
  tmp<-cbind(tmp,data.frame(group=rep("Month 1",nrow(tmp))))
  l[[i]]<-tmp[,c("contig","contig.length","life","group","vlp.reads.number")]
}
tmp2<-do.call("rbind", l)


samples<-list.files(paste(data.pa,"discovery.virome/my.analysis/life",sep=""),pattern = "4$")
l<-c()
for (i in samples){ 
  tmp<-get(paste(i,".contig.ann.life",sep="")) %>% filter(source=="Viruses",uniprot.total.orf.number>=10)  
  tmp$Temperate<-as.numeric(as.character(tmp$Temperate))
  tmp$Lytic<-as.numeric(as.character(tmp$Lytic))
  tmp<-mutate(tmp,
              life=Temperate-0.5)
  tmp<-cbind(tmp,data.frame(group=rep("Month 4",nrow(tmp))))
  l[[i]]<-tmp[,c("contig","contig.length","life","group","vlp.reads.number")]
}
tmp3<-do.call("rbind", l)
tmp<-rbind(tmp1,tmp2,tmp3)
r1<-round(sum(tmp[tmp$group=="Month 0"&tmp$life<=0,"vlp.reads.number"])/sum(tmp[tmp$group=="Month 0","vlp.reads.number"]),digits = 3)
r2<-round(sum(tmp[tmp$group=="Month 1"&tmp$life<=0,"vlp.reads.number"])/sum(tmp[tmp$group=="Month 1","vlp.reads.number"]),digits = 3)
r3<-round(sum(tmp[tmp$group=="Month 4"&tmp$life<=0,"vlp.reads.number"])/sum(tmp[tmp$group=="Month 4","vlp.reads.number"]),digits = 3)




l<-c()
for (i in 1:nrow(tmp))
{
  tmp1=data.frame(reads.count=rep(1,round(tmp[i,5]/min(tmp[,5]))),contig=rep(tmp[i,1],round(tmp[i,5]/min(tmp[,5]))),contig.length=rep(tmp[i,2],round(tmp[i,5]/min(tmp[,5]))),life=rep(tmp[i,3],round(tmp[i,5]/min(tmp[,5]))),group=rep(tmp[i,4],round(tmp[i,5]/min(tmp[,5]))))
  l[[i]]<-tmp1
}
#tmp2<-do.call("rbind", l) # too slow
tmp2<-rbind.fill(l)# faster
tmp2<-separate(tmp2,col = contig,sep = "\\.",into = "contig")
tmp2<-mutate(tmp2,subject=substr(contig,2,4)) # substring a whole column


pdf("output/Fig.S3.pdf",width=10,height=10)
ggplot(tmp2, aes(x=life,fill=group)) +
  geom_density()+

  scale_fill_manual(values = age.col)+
  coord_cartesian(xlim = c(-0.14, 0.14)) +
   labs(x="",y="Reads density")+
  geom_vline(xintercept = 0)+
  scale_x_continuous(breaks=c(-0.14,-0.07,0,0.07,0.14),
                     labels=c("1","0.75","0.5","0.75","1"))+ # scale the probability 
  facet_grid(group ~ .)+
  theme_classic()+
  theme(
    axis.title = element_text(size=20,colour="black"),
    axis.text = element_text(size=15,colour="black"),
    strip.background = element_blank(),
    strip.text=element_blank(),
    legend.position="none"
  )

grid.text(percent(r1), x = unit(0.45, "npc"), y = unit(0.90, "npc"),gp=gpar(fontsize=20, col="black"))
grid.text(percent(1-r1), x = unit(0.62, "npc"), y = unit(0.90, "npc"),gp=gpar(fontsize=20, col="black"))
grid.text(percent(r2), x = unit(0.45, "npc"), y = unit(0.60, "npc"),gp=gpar(fontsize=20, col="black"))
grid.text(percent(1-r2), x = unit(0.62, "npc"), y = unit(0.60, "npc"),gp=gpar(fontsize=20, col="black"))
grid.text(percent(r3), x = unit(0.45, "npc"), y = unit(0.3, "npc"),gp=gpar(fontsize=20, col="black"))
grid.text(percent(1-r3), x = unit(0.62, "npc"), y = unit(0.3, "npc"),gp=gpar(fontsize=20, col="black"))
grid.text("Month 0", x = unit(0.18, "npc"), y = unit(0.95, "npc"),gp=gpar(fontsize=20, col="black"))
grid.text("Month 1", x = unit(0.18, "npc"), y = unit(0.65, "npc"),gp=gpar(fontsize=20, col="black"))
grid.text("Month 4", x = unit(0.18, "npc"), y = unit(0.35, "npc"),gp=gpar(fontsize=20, col="black"))
grid.text("Lytic", x = unit(0.1, "npc"), y = unit(0.02, "npc"),gp=gpar(fontsize=20, col="black"))
grid.text("Temperate", x = unit(0.92, "npc"), y = unit(0.02, "npc"),gp=gpar(fontsize=20, col="black"))
grid.text("Probability", x = unit(0.53, "npc"), y = unit(0.02, "npc"),gp=gpar(fontsize=20, col="black"))


grid.lines(x = unit(c(0.25, 0.45), "npc"),
           y = unit(c(0.02, 0.02), "npc"),
           gp = gpar(fill="black"),
           arrow = arrow(length = unit(0.2, "inches"), angle = 5,
                         ends="first", type="closed"))
grid.lines(x = unit(c(0.6, 0.8), "npc"),
           y = unit(c(0.02, 0.02), "npc"),
           gp = gpar(fill="black"),
           arrow = arrow(length = unit(0.2, "inches"), angle = 5,
                         ends="last", type="closed"))
dev.off()
#### End ####

#### Fig. S2 Herv reads####
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


### Human herv
#tmp<-repeatp[(!repeatp$group2=="Control")&(repeatp$library=="DNA"),c("group","hervp")]
#tmp<-repeatp[repeatp$library=="DNA",c("group","hervp")] # this for include negative control
tmp<-(repeatp %>% filter(study_group=="Month 0"|study_group=="Month 1"|study_group=="Month 4"|study_group=="Negative Control",library_type=="DNA"))[,c("study_group","hervp")]
tmp1<-aggregate(. ~ study_group,tmp, median)
#tmp2<-aggregate(. ~ group,tmp, sd)
tmp2<-aggregate(. ~ study_group,tmp, function(x) {sd(x)/sqrt(length(x))}) # this is for standard erro
tmp<-cbind(tmp1,tmp2[,2])
names(tmp)<-c("Age","mean","se")
tmp$library<-"DNA"

tmp3<-(repeatp %>% filter(study_group=="Month 0"|study_group=="Month 1"|study_group=="Month 4"|study_group=="Negative Control",library_type=="RNA"))[,c("study_group","hervp")]
#tmp3<-repeatp[repeatp$library=="RNA",c("group","hervp")] # this is for include negative control
tmp1<-aggregate(. ~ study_group,tmp3, median)
#tmp2<-aggregate(. ~ group,tmp, sd)
tmp2<-aggregate(. ~ study_group,tmp3, function(x) {sd(x)/sqrt(length(x))}) # this is for standard erro
tmp3<-cbind(tmp1,tmp2[,2])
names(tmp3)<-c("Age","mean","se")
tmp3$library<-"RNA"
tmp<-rbind(tmp,tmp3)
ggplot(tmp,aes(Age,100*mean,width=0.5,fill=library))+
  geom_bar(stat="identity", color="black",size=1,position=position_dodge())+
  geom_errorbar(aes(ymin=mean*100, ymax=(mean+se)*100), width=0.25,size=1,
                position=position_dodge(0.5)) +
  scale_y_continuous(expand = c(0, 0)) +
  #scale_x_discrete(breaks=c("Month 0","Month 1","Month 4"),
  #                 labels=c("0", "1", "4"))+
  scale_fill_manual(values=cohort.col)+
  #geom_hline(yintercept = 6.67E+06,color="bisque3",linetype="dashed",size=2)+ # this is for dash line for limitation
  ylab("HERV reads percentage")+
  xlab("")+
  coord_cartesian(ylim = c(0,10.2)) +
  theme_classic()+
  theme(axis.title = element_text(size=20,colour="black"),
        axis.text = element_text(size=15,colour = "black"),
        legend.direction = "horizontal",
        legend.position = c(0.8,0.8),
        legend.text = element_text(size = 15,colour = "black"),
        legend.title = element_blank()
  )
ggsave("output/Fig.S2.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 7, height = 5, units = c("in"),
       dpi = 300, limitsize = TRUE)
#### End ####

#### Fig. S4 crAssphage ####
meta <- read.delim(paste(meta.pa,"meta.data.txt",sep=""))
veo.meta<-read.delim(paste(meta.pa,"veo.meta.txt",sep=""))
read_bwa_crass <- function(filename1,filename2) {
  p = read.delim(filename1, header=F,
                 stringsAsFactors = FALSE)
  names(p)<-c("family", 
              "read",
              "align.len","genome.len", "cov")
  p<-p[,c(1,2,5)]
    x<- p %>% filter(read>0)
  if (nrow(x) == 0) {
    x<-p[1,c("family","cov")]
    x$cov<-0
  } else {
    x<- aggregate(.~family,x[,c("family","cov")],sum)
 
    y = read.delim(filename2, header=F,
                   stringsAsFactors = FALSE)
    names(y)<-c("family", 
                "genome.len",
                "mapped","unmapped")
    y<-y %>% mutate (total.read=sum(mapped)+sum(unmapped))
    y<-y[-nrow(y),c(1,2,3,5)]
    z<-merge(x,y,all.y=T,by="family")
    z<-z %>% mutate (rpkm=mapped/(genome.len/1000*total.read/1000000))
    z<- z %>% arrange(desc(cov))
    z<-z[,c(1,2,4,6)]
    z[is.na(z$cov),"cov"] <- 0
    z$sample_id<-strsplit(strsplit(filename1,"\\/")[[1]][length(strsplit(filename1,"\\/")[[1]])],"\\.")[[1]][1]
    names(z)<-c("crass.family","cov","reads","rpkm","library_id")
    return(z)
  }
}
# function to make table 
table_bwa_crass<- function(cohort,viralname){
  path1 = paste(data.pa,"coverage/",cohort,"/",viralname,sep="")
  filename.cov <- list.files(path = path1,full.names = TRUE,pattern = ".coverage")
  filename.reads <- list.files(path = path1,full.names = TRUE,pattern = ".read")
  bwa.cov.data <- mapply(read_bwa_crass,filename.cov,filename.reads,SIMPLIFY = FALSE) # mapply do not give list default as lapply, need simplify option
  bwa.cov<-do.call(rbind,bwa.cov.data)
  rownames(bwa.cov)<-NULL
  return(bwa.cov)
}  

dis.crass<-table_bwa_crass("dis","crassphage")
dis.crass<-merge(dis.crass,meta[,c("library_id","study_group")],by="library_id")
dis.crass$Cohort<-"Discovery"
dis.crass$Family<-"crAssphage"

vali.crass<-table_bwa_crass("vali","crassphage")
vali.crass<-merge(vali.crass,meta[,c("library_id","study_group")],by="library_id")
vali.crass$Cohort<-"Validation"
vali.crass$Family<-"crAssphage"

veo.crass<-table_bwa_crass("veo","crassphage")
veo.crass<-merge(veo.crass,veo.meta[,c("library_id","study_group")],by="library_id")
veo.crass$Cohort<-"VEO"
veo.crass$Family<-"crAssphage"

l<-list(dis.crass,vali.crass,veo.crass)
tmp<-do.call(rbind,l)

tmp$inf<-"Uninfected"
tmp[tmp$cov>=(1/3),"inf"]<-"Infected"
tmp$inf<-factor(tmp$inf,levels=c("Infected","Uninfected"))

tmp<-tmp %>% 
  filter(Cohort=="Discovery"|Cohort=="VEO",
         study_group=="Month 0"|study_group=="Month 1"|study_group=="Month 4"|study_group=="healthy"
  )
tmp<-tmp[which(tmp$crass.family %in% unique(tmp[tmp$inf=="Infected",]$crass.family)),]
tmp$group<-gsub(tmp$study_group,pattern = "healthy",replacement = "Year 2~5")
tmp$crass.family<-gsub(tmp$crass.family,
                       pattern = c("CDZN01024782"),
                       replacement = c("crAssphage member 1"))
tmp$crass.family<-gsub(tmp$crass.family,
                       pattern = c("CEAR01029167"),
                       replacement = c("crAssphage member 2"))
tmp$crass.family<-gsub(tmp$crass.family,
                       pattern = c("Chlamydia_CVNZ01000019ext"),
                       replacement = c("crAssphage member 3"))

tmp$crass.family<-factor(tmp$crass.family,levels=c("crAssphage member 3","crAssphage member 2","crAssphage member 1","crAssphage"))
g1<-  ggplot(tmp,aes(x = library_id, y = crass.family)) +
  geom_tile(aes(fill=inf),color="grey", size=0.1)+
  scale_fill_manual(values = c("black","white"))+
  labs(x="",y="")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.ticks =element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=15,colour="black"),
        strip.background = element_rect(linetype = 0),
        strip.text.x  = element_text(size=20,colour="black"),
        strip.text.y  = element_text(size=15,colour="black",angle = 360),
        legend.position = "none"
        
  )+
  facet_grid(.~study_group,scales = "free", space = "free")+
  theme(aspect.ratio = 1)


g <- ggplot_gtable(ggplot_build(g1))
stript <- which(grepl('strip-t', g$layout$name))
fills <- age.col
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

pdf("output/Fig.S4.pdf", width = 20, height = 5)
grid.draw(g)
dev.off()

#### End ####


# rm(list=ls())
# options(warn=0)
