#### Set up working directory ####
data.pa<-"./input/wet.input/"
meta.pa<-"./input/meta.data/"
database.pa<-"./input/database/"
options(warn=-1)
#### End ####

#### Libraries ####
library(ggplot2)
library(dplyr)
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
# 
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

#### Fig.2B Induction phage number ####
induce<-read.delim(paste(data.pa,"induced.phage.number.txt",sep=""))
induce$inducer<-factor(induce$inducer,levels=c("No iducer","Mitomycin C"))
induce<-induce %>% arrange(number)
induce$id<-factor(induce$id,levels=unique(as.character(induce$id)),ordered = T)

g1<-ggplot(induce, aes(x = inducer, y = id, fill = log(number+1,10))) +
  geom_tile(color="grey", size=0.1) +
  labs(x="",y="")+
  facet_grid(Genus~condition,scales = "free",space = "free",switch="y")+
  scale_fill_viridis(name="Particle number per ml media",guide=guide_colorbar(title.position = "top"),
                     breaks = c(log(1,10),log(101,10),log(10001,10),log(1000001,10),log(100000001,10),log(10000000001,10) ),
                     labels=c(0,expression("10"^2),expression("10"^4),expression("10"^6),expression("10"^8),expression("10"^10)))+
  theme(aspect.ratio = 1)+
  theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        axis.ticks =element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.background.x = element_blank(),
        strip.background.y = element_rect(fill="Black",color="white"),
        strip.text.x  = element_text(size=25,colour="black"),
        strip.text.y  = element_text(size=0,colour="white",angle = 180),
        legend.position="bottom",
        legend.text=element_text(size=15),
        legend.title = element_text(size=15,margin = margin(t = 0, unit = "pt")),
        legend.key=element_rect(color='black'),
        legend.key.height = unit(0.1, "inch"),
        legend.key.width=unit(0.5, "inch"),
        plot.margin = margin(l=100)
  )




g<-ggplotGrob(g1) 
stripr <- which(grepl('strip-l', g$layout$name))
fills <- color.tron[c(1:3,5:7)]
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

pdf("output/Fig.2B.pdf",width = 7, height = 7)
grid.draw(g)
grid.text(levels(induce$Genus),
          gp=gpar(fontsize=17), 
          x = c(0.135,0.12,0.115,0.130,0.145,0.115), 
          y = c(0.93,0.813,0.565,0.365,0.265,0.20))
grid.text(c("Inducer","-","+","-","+"),
          gp=gpar(fontsize=17), 
          x = c(0.165,0.35,0.52,0.72,0.9), 
          y = rep(0.16,5))
dev.off()


#### End #### 

#### Fig.S1 Totalshotgun ####
meta <- read.delim(paste(meta.pa,"meta.data.txt",sep=""))
kraken <- read.delim(paste(data.pa,"shotgun.metagenome/sunbeam_output/classify/all_samples.tsv",sep=""), skip = 1)
names(kraken) <- gsub(x = names(kraken),
                      pattern = "\\.taxa",
                      replacement = "\\")
names(kraken) <- gsub(x = names(kraken),
                      pattern = "\\.",
                      replacement = "\\")
names(kraken)[1]<-"tax_id" ## kraken is table with 1 is tax id, 2:157 are libraries, 158 is names
lin <- read.delim(paste(data.pa,"shotgun.metagenome/sunbeam_output/classify/lineages.txt",sep="")) # lineages data by tax to baltimore
host <- read.delim("input/meta.data/viral.host.txt") # virus host data
kraken<-merge(kraken,lin,by="tax_id") # merge kraken and lineage data

shotgun.sam<-meta %>% 
  filter(Cohort=="Discovery",study_group=="Month 0"|study_group=="Month 1"|study_group=="Month 4",library_type=="shotgun") 

kraken<-
  kraken[,c("superkingdom",as.character(shotgun.sam$library_id))]
kraken$superkingdom[1]<-"Bacteria" # the ifrst row did not show anything, its bacteria
tmp<-aggregate(.~superkingdom,kraken,sum)

filename<-list.files(path =paste(data.pa,"shotgun.metagenome/sunbeam_output/qc/log",sep="") ,full.names = TRUE,pattern = "1.txt")
read_log<-function(x){
y<-read.delim(x)
y<-y[,c(1,4)]
y$library_id<-strsplit(strsplit(x,"\\/")[[1]][8],"\\_")[[1]][1]
return(y)
}
l<-lapply(filename,read_log)
data<-do.call(rbind,l)
names(data)<-c("Human", "nonhost" ,"library_id")

tmp1<-as.data.frame(t(tmp[,-1])) # transpose number and string together will have problem
tmp1$V4<-rownames(tmp1)
rownames(tmp1)<-NULL
names(tmp1)<-c(as.character(tmp$superkingdom), "library_id")


tmp<- merge(tmp1,data,by="library_id")
tmp<-tmp%>%mutate(Archaea.p=Archaea/(nonhost+Human),
                  Bacteria.p=Bacteria/(nonhost+Human),
                  Viruses.p=Viruses/(nonhost+Human),
                  Human.p=Human/(nonhost+Human),
                  Unassigned=1-(Archaea+Bacteria+Viruses+Human)/(nonhost+Human)
                  )
tmp<-tmp[,c("library_id","Archaea.p","Bacteria.p","Viruses.p","Unassigned","Human.p")]                 
names(tmp)<-c("library_id","Archaea","Bacteria","Viruses","Unassigned","Human")

data<-melt(tmp)
names(data)<-c("library_id","Source","Percentage")
# merge with prop data come from # total shot gun human reads vs 16s
data<-merge(data,meta[,c("library_id","study_group")],by="library_id")
data<-arrange(data,Source,Percentage) # sort by Source then by Percentage
data$library_id<-factor(data$library_id,levels=unique(data[data$Source=="Bacteria",]$library_id),ordered = T) # change sample id order by bacteria percentage
data$Source<-factor(data$Source,levels=c("Bacteria","Archaea","Viruses","Human","Unassigned"))

ggplot(data,aes(y = Percentage*100, x = interaction(library_id,study_group), fill = Source))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=color[1:6])+
  labs(x ="",y="")+
  facet_wrap(~study_group,scales ="free_x",nrow = 1 )+
  theme(
    axis.text.x=element_blank(),
    axis.text.y = element_text(size=15,colour="black"),
    axis.ticks.x=element_blank(),
    panel.background = element_blank(),
    legend.title=element_blank(),
    legend.text = element_text(size=15,colour="black"),
    strip.background = element_blank(),
    strip.text = element_text(size=20,colour="black") 
  )

ggsave("output/Fig.S1.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 7, height = 3.5, units = c("in"),
       dpi = 300, limitsize = TRUE)

#### End ####

#### Fig.2,C,D,E ####
# funtions needed
uniprottax<-read.delim(paste(database.pa,"uniprot.virome/uniprot.virome.meta",sep=""),header=F)
names(uniprottax)<-c("uniprot.hit.id","taxid","name")
uniprottax$taxid<-as.numeric(as.character(uniprottax$taxid))
wgs.contig.annotation<-function (x){
  pa=paste(data.pa,"wgs/sunbeam_output/assembly/sspace/",sep="")
  # Contig (>3k) reads number (decontam read map contig) ###
  datapath <- file.path(paste(pa,x,"/contigs/",x,".contig.reads.txt",sep=""))
  map <- read.delim(datapath,sep="\t",header=F)
  names(map)<-c("contig","contig.length","vlp.reads.number","unmmaped.reads.number")
  #data<-map[-nrow(map),] this command will remove only the last row, we will remove all rows lower than 3000
  contig <- map %>% filter(contig.length>=3000) %>%
    mutate(fpkm=vlp.reads.number*10e9/contig.length/sum(map$vlp.reads.number))
  
  # Contig (>3k) coverage data (decontam read map contig) ###
  datapath <- file.path(paste(pa,x,"/contigs/",x,".contig.coverage.txt",sep=""))
  cov <- read.delim(datapath,sep="\t",header=F)
  names(cov)<-c("contig","vlp.coverage")
  contig<-merge(contig,cov,all.x=T,by="contig")
  
  # Contig orf number data ###
  datapath <- file.path(paste(pa,x,"/contigs/",x,".orf.number",sep=""))
  orfn<-read.delim(datapath,sep="\t",header=F)
  orfn<-separate(data = orfn, col = V1, into = c("contig"), sep = "\\_")
  tmp<-data.frame(table(orfn$contig))
  names(tmp)<-c("contig","orf.number")
  contig<-merge(contig,tmp,all.x=T,by="contig")
  
  # Uniprot data process ###
  datapath <- file.path(paste(pa,x,"/contigs/",x,".3k.uniprot",sep=""))
  uniprot<-read.delim(datapath,header=F,sep="\t")
  names(uniprot)<-c("contig", "sseqid", "pident", "qlen","slen","alignment.len" ,"mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  uniprot<- uniprot  %>% 
    separate(col = contig, into = c("contig","orfs.id"), sep = "\\_")  %>%  
    separate(col = sseqid, into = c("1","2","uniprot.hit.id"), sep = "\\|")
  uniprot<-merge(uniprot,uniprottax,by="uniprot.hit.id")
  uniprot<-cbind(uniprot,getTaxonomy(uniprot$taxid,"./input/database/ac2tax/accessionTaxa.sql"))
  # Uniprot orf number ###
  tmp<-data.frame(table(uniprot$contig))
  names(tmp)<-c("contig","uniprot.total.orf.number")
  contig<-merge(contig,tmp,all.x=T,by="contig")
  contig$uniprot.total.orf.number[is.na(contig$uniprot.total.orf.number)]<-0
  # Calculate uniprot species ###
  tmp<-uniprot[,c("contig","species")]
  tmp<-table(tmp)
  tmp1<-data.frame(colnames(tmp)[apply(tmp,1,which.max)])
  tmp<-uniprot[,c("contig","species")]
  tmp<-data.frame(table(tmp))
  tmp2<-data.frame(tapply(tmp$Freq, tmp$contig, max))
  tmp3<-cbind(tmp1,tmp2)
  names(tmp3)<-c("uniprot.species","uniprot.species.orf.number")
  tmp3$contig<-rownames(tmp3)
  rownames(tmp3)<-NULL
  contig<-merge(contig,tmp3,all.x=T,by="contig")
  contig$uniprot.species.orf.number[is.na(contig$uniprot.species.orf.number)]<-0
  
  # Prophage data process ###
  datapath <- file.path(paste(pa,x,"/contigs/",x,".3k.prophage",sep=""))
  prophage<-read.delim(datapath,header=F,sep="\t")
  names(prophage)<-c("contig", "sseqid", "pident", "qlen","slen","alignment.len" ,"mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  prophage<- prophage  %>% 
    separate(col = contig, into = c("contig","orfs.id"), sep = "\\_")  %>%  
    separate(col = sseqid, into = c("1","2","3","prophage.hit.id"), sep = "\\|")
  
  
  taxaId<-accessionToTaxa(prophage$prophage.hit.id,"./input/database/ac2tax/accessionTaxa.sql")
  prophage<-cbind(prophage,getTaxonomy(taxaId,"./input/database/ac2tax/accessionTaxa.sql"))
  # Prophage orf number ###
  tmp<-data.frame(table(prophage$contig))
  names(tmp)<-c("contig","prophage.total.orf.number")
  contig<-merge(contig,tmp,all.x=T,by="contig")
  contig$prophage.total.orf.number[is.na(contig$prophage.total.orf.number)]<-0
  # Calculate uniprot species ###
  tmp<-prophage[,c("contig","species")]
  tmp<-table(tmp)
  tmp1<-data.frame(colnames(tmp)[apply(tmp,1,which.max)])
  tmp<-prophage[,c("contig","species")]
  tmp<-data.frame(table(tmp))
  tmp2<-data.frame(tapply(tmp$Freq, tmp$contig, max))
  tmp3<-cbind(tmp1,tmp2)
  names(tmp3)<-c("prophage.species","prophage.species.orf.number")
  tmp3$contig<-rownames(tmp3)
  rownames(tmp3)<-NULL
  contig<-merge(contig,tmp3,all.x=T,by="contig")
  contig$prophage.species.orf.number[is.na(contig$prophage.species.orf.number)]<-0
  
  
  contig <- contig %>% mutate(uniprot.orf.ratio=uniprot.total.orf.number/orf.number,
                              uniprot.species.orf.ratio=uniprot.species.orf.number/uniprot.total.orf.number,
                              prophage.orf.ratio=prophage.total.orf.number/orf.number,
                              prophage.species.orf.ratio=prophage.species.orf.number/prophage.total.orf.number)
  
  # final output ####
  l<-list(contig,uniprot,prophage)
  return(l)
  
}
coverage<-function (sample,genome,path,wgs.contig){
  bg<-read.delim(paste(path,genome,"/map/",sample,".bg",sep=""),header=F) # induction data
  names(bg)<-c("contig","start","end","coverage")
  bg<-merge(bg,wgs.contig[,c("contig","contig.length")])
  bg$contig<-droplevels(bg$contig)
  bg<-bg[with(bg, order(contig, start)), ]
  
  bg1=NULL
  for (i in 1:length(levels(bg$contig))) {
    for (j in 1:nrow(bg[bg$contig==levels(bg$contig)[i],])) {
      poseq<-seq(bg[bg$contig==levels(bg$contig)[i],]$start[j],bg[bg$contig==levels(bg$contig)[i],]$end[j],by=100)
      bg1<-rbind(bg1,data.frame(contig=rep(levels(bg$contig)[i],length(poseq)),
                                pos=poseq,
                                coverage=rep(bg[bg$contig==levels(bg$contig)[i],"coverage"][j],length(poseq))))
      
    }
  }
  bg1<-bg1[!duplicated(bg1[c("pos","contig")]),]
  bg1$pos<-as.numeric(as.character(bg1$pos))
  
  return(bg1)
}
prophage<-function (genome,path){
  orf.pro<-wgs.contig.annotation(genome)[[3]]
  orf.pro<-orf.pro[,c("contig","orfs.id","prophage.hit.id")]
  orf.coords<-read.delim(paste(path,genome,"/contigs/orf.coords",sep=""),header=F)
  orf.coords <- orf.coords %>% separate(col = V1, into = c("contig","orfs.id"), sep = "\\_")
  orf.coords$orfs.id<-gsub(" ", "", orf.coords$orfs.id)
  pro<-left_join(orf.pro,orf.coords,by=c("contig","orfs.id"))
  
  pro<-pro %>% mutate(protein=ifelse(!is.na(prophage.hit.id),1,0))
  pro<-pro[,c("contig","protein","V2","V3")]
  pro<-pro[with(pro, order(contig, V2)), ]
  pro$contig<-as.factor(pro$contig)
  pro$contig<-droplevels(pro$contig)
  names(pro)<-c("contig","protein","start","end")
  pro<-pro %>% drop_na(start)
  
  tmp<-NULL
  tmp1<-NULL
  tmp2<-NULL
  contig<-wgs.contig.annotation(genome)[[1]]
  for (i in 1:length(levels(pro$contig))) {
    tmp1<-seq(0,contig[contig$contig==levels(pro$contig)[i],"contig.length"]+10000,by=10000)
    for (k in 1:(length(tmp1)-1))
    {
      tryCatch({
        tmp2<-aggregate(.~contig,pro[pro$contig==levels(pro$contig)[i]&pro$start<tmp1[k+1]&pro$start>=tmp1[k],],sum)[,1:2]
        tmp2$pos<-tmp1[k+1]
        tmp<-rbind(tmp,tmp2)
      }, error=function(e){})      
    }
  }
  
  pro1=NULL
  for (i in 1:length(levels(pro$contig))) {
    for (j in 1:nrow(pro[pro$contig==levels(pro$contig)[i],])) {
      poseq<-seq((j-1)*10000,j*10000,by=100)
      pro1<-rbind(pro1,data.frame(contig=rep(levels(pro$contig)[i],length(poseq)),
                                  pos=poseq,
                                  protein=rep(tmp[tmp$contig==levels(tmp$contig)[i],"protein"][j],length(poseq))))
    }
  }
  pro1<-pro1[with(pro1, order(contig, pos)), ]
  pro1<-pro1[!duplicated(pro1[c("pos","contig")]),]
  pro1$pos<-as.numeric(as.character(pro1$pos)) 
  names(pro1)<-c("Contig","Position","Protein")
  pro1[is.na(pro1$Protein)|pro1$Protein<=3,"Protein"]<-0 # lower than 3 assaign 0
  return(pro1)
  
}
coverage.plot<-function (com){
  
  com<-com[with(com, order(Contig, Position)), ]
  com[is.na(com)]<-0
  com<-com[! com$Contig %in% levels(tmp$Contig), ]
  com$Contig<-droplevels(com$Contig)
  com$Contig<-factor(com$Contig,levels=mixedsort(levels(com$Contig)))
  mtmp<-melt(com,id=c("Contig","Position"))
  mtmp$Contig<-factor(mtmp$Contig,levels=mixedsort(levels(mtmp$Contig)))
  color=pal_npg(palette = c("nrc"), alpha = 1)(10)
  #### figure  ####
  
  p<-ggplot(mtmp, aes(x=Position, y=value,color=variable)) + 
    geom_line(size=2)+
    scale_x_continuous(breaks=c(0,100000),
                       labels=trans_format("log10", math_format(10^.x)))+
    scale_color_manual(values=color[c(1:3,5)])+
    facet_grid(variable~Contig,scales = "free",space="free_x",switch = "y")+
    #ggtitle(unlist(strsplit(as.character(com$Contig[1]),"\\."))[1]) +
    theme_bw()+
    theme(
      panel.border = element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.line.y = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title = element_blank(),
      legend.position='none',
      strip.background = element_rect(color="white"),
      #strip.background =element_blank(),
      strip.text.x = element_blank(),
      strip.text.y = element_text(size=30,angle = 180)
    )
  return(p) 
  
}
genome.plot<-function (wgs){
  genome.color<-brewer.pal(n=11, name = 'Paired')
  names(genome.color)<-c("Attachment site",
                         "Head protein",
                         "Hypothetical protein",
                         "Integrase",
                         "Lysis protein",
                         "Phage-like protein",
                         "Portal protein",
                         "Tail protein",
                         "Terminase",
                         "tRNA","Fiber protein"
  )
  tmp<-read.delim(paste(data.pa,"wgs/sunbeam_output/assembly/sspace/",wgs,".prophage.genome.map.txt",sep=""))
  tmp$protein_name<-factor(tmp$protein_name,levels=c("Attachment site",
                                                     "Head protein",
                                                     "Hypothetical protein",
                                                     "Integrase",
                                                     "Lysis protein",
                                                     "Phage-like protein",
                                                     "Portal protein",
                                                     "Tail protein",
                                                     "Terminase",
                                                     "tRNA","Fiber protein"
  ))
  p<-ggplot(tmp) +
    geom_rect(aes(xmin = start, xmax = end, fill = protein_name, color = protein_name,ymin=ymin, ymax=ymax)) +
    scale_fill_manual(values=genome.color,drop=FALSE)+
    scale_color_manual(values=genome.color,drop=FALSE)+
    geom_segment(aes(x=min(tmp$start),xend=max(tmp$end),y=0,yend=0),color="grey50")+
    annotate("text",x=c(min(tmp$start)-500,max(tmp$end)+1000),y = c(0,0),label=c("0",paste(round((max(tmp$end)-min(tmp$start))/1000),"K",sep="")))+
    annotate("segment", x = c(quantile(tmp$start, .35),quantile(tmp$start, .65)),
             xend = c(quantile(tmp$start, .65),quantile(tmp$start, .35)), 
             y =c (1.2,-1.2), yend = c(1.2,-1.2), colour = "black", size=0.5,
             alpha=0.8, arrow=arrow(ends = "last",angle = 5,type = "closed"))+
    guides(fill= guide_legend(nrow = 2))+
    theme_classic()+
    theme(axis.line=element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.title = element_blank(),
          legend.position = "bottom"
          
    )
  return(p)
}

### all files saved in path
path=paste(data.pa,"wgs/sunbeam_output/assembly/sspace/",sep="")


# wgs1 for genome coverage and phage genome annotation ####
wgs1.contig<-wgs.contig.annotation("wgs1")[[1]]
BM3451.cov<-coverage("BM3451","wgs1",path,wgs1.contig)
names(BM3451.cov)<-c("Contig","Position","Mitomycin C")
wgs1.D3451.cov<-coverage("D3451","wgs1",path,wgs1.contig)
names(wgs1.D3451.cov)<-c("Contig","Position","Purified VLP")
BN3451.cov<-coverage("BN3451","wgs1",path,wgs1.contig)
names(BN3451.cov)<-c("Contig","Position","No Inducer")
wgs1.pro<-prophage("wgs1",path)
names(wgs1.pro)<-c("Contig","Position","Prophage protein number")
mito<-BM3451.cov
nomito<-BN3451.cov
vlp<-wgs1.D3451.cov
pro<-wgs1.pro

com<-NULL
com<-Reduce(function(x, y) merge(x, y, all.x=TRUE,by=c("Contig","Position")), list(mito,nomito,vlp,pro)) # merge multiple dataframe
wgs1.p<-coverage.plot(com)


pdf("output/Fig.2C.1.pdf",width = 30, height = 8)
color=pal_npg(palette = c("nrc"), alpha = 1)(10)
wgs1.p
g <- ggplot_gtable(ggplot_build(wgs1.p))
stripl <- which(grepl('strip-l', g$layout$name))
fills <- color[c(1:3,5)]
k <- 1
for (i in stripl) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)
grid.rect(x = 0.593, y = 0.5,
          width = 0.03, height = 0.98,
          just = "centre", hjust = NULL, vjust = NULL,
          default.units = "npc", name = NULL,
          gp=gpar(lwd = 4), draw = TRUE, vp = NULL)
dev.off()

g<-genome.plot("wgs1")
pdf("output/Fig.2C.2.pdf",width = 10, height = 2)
g
dev.off()
# end ####
# wgs2 for genome coverage and phage genome annotation ####
wgs2.contig<-wgs.contig.annotation("wgs2")[[1]]
BC3751.cov<-coverage("BC3751","wgs2",path,wgs2.contig)
names(BC3751.cov)<-c("Contig","Position","Mitomycin C")
BN3751.cov<-coverage("BN3751","wgs2",path,wgs2.contig)
names(BN3751.cov)<-c("Contig","Position","No Inducer")
D3751.cov<-coverage("D3751","wgs2",path,wgs2.contig)
names(D3751.cov)<-c("Contig","Position","Purified VLP")
wgs2.pro<-prophage("wgs2",path)
names(wgs2.pro)<-c("Contig","Position","Prophage protein number")
mito<-BC3751.cov
nomito<-BN3751.cov
vlp<-D3751.cov
pro<-wgs2.pro


com<-NULL
com<-Reduce(function(x, y) merge(x, y, all.x=TRUE,by=c("Contig","Position")), list(mito,nomito,vlp,pro)) # merge multiple dataframe
wgs2.p<-coverage.plot(com)


pdf("output/Fig.2D.1.pdf",width = 30, height = 8)
color=pal_npg(palette = c("nrc"), alpha = 1)(10)
wgs2.p
g <- ggplot_gtable(ggplot_build(wgs2.p))
stripl <- which(grepl('strip-l', g$layout$name))
fills <- color[c(1:3,5)]
k <- 1
for (i in stripl) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)
grid.rect(x = 0.55, y = 0.5,
          width = 0.03, height = 0.98,
          just = "centre", hjust = NULL, vjust = NULL,
          default.units = "npc", name = NULL,
          gp=gpar(lwd = 4), draw = TRUE, vp = NULL)
dev.off()

g<-genome.plot("wgs2")
pdf("output/Fig.2D.2.pdf",width = 10, height = 2)
g
dev.off()
# end ####
# wgs3 for genome coverage and phage genome annotation ####
wgs3.contig<-wgs.contig.annotation("wgs3")[[1]]
BM3754.cov<-coverage("BM3754","wgs3",path,wgs3.contig)
names(BM3754.cov)<-c("Contig","Position","Mitomycin C")
mito<-BM3754.cov
BN3754.cov<-coverage("BN3754","wgs3",path,wgs3.contig)
names(BN3754.cov)<-c("Contig","Position","No Inducer")
nomito<-BN3754.cov
wgs3.D3754.cov<-coverage("D3754","wgs3",path,wgs3.contig)
names(wgs3.D3754.cov)<-c("Contig","Position","Purified VLP")
vlp<-wgs3.D3754.cov
wgs3.pro<-prophage("wgs3",path)
names(wgs3.pro)<-c("Contig","Position","Prophage protein number")
pro<-wgs3.pro
com<-NULL
com<-Reduce(function(x, y) merge(x, y, all.x=TRUE,by=c("Contig","Position")), list(mito,nomito,vlp,pro)) # merge multiple dataframe
wgs3.p<-coverage.plot(com)


pdf("output/Fig.2E.1.pdf",width = 30, height = 8)
color=pal_npg(palette = c("nrc"), alpha = 1)(10)
wgs3.p
g <- ggplot_gtable(ggplot_build(wgs3.p))
stripl <- which(grepl('strip-l', g$layout$name))
fills <- color[c(1:3,5)]
k <- 1
for (i in stripl) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)
grid.rect(x = 0.285, y = 0.5,
          width = 0.03, height = 0.98,
          just = "centre", hjust = NULL, vjust = NULL,
          default.units = "npc", name = NULL,
          gp=gpar(lwd = 4), draw = TRUE, vp = NULL)
dev.off()

g<-genome.plot("wgs3")
pdf("output/Fig.2E.2.pdf",width = 10, height = 2)
g
dev.off()
# end ####
#### End ####

#### Fig.2F and 2G, S6 ####
pa<-paste(data.pa,"induced.virome/my.analysis/induced.vlp.to.stool.vlp",sep="")
names<-grep(list.files(path=pa), pattern='.txt', inv=T, value=T) # exclude .txt files
read_idx<-function(x){
  y<-read.delim(x, header=F)
  names(y)<-c("contig", 
              "contig.len",
              "mapped","unmapped")  
  y<-y %>% mutate (total.read=sum(mapped)+sum(unmapped))
  y<-y[-nrow(y),c("contig","mapped","total.read")]
  y$vlp<-strsplit(strsplit(x,"\\/")[[1]][length(strsplit(x,"\\/")[[1]])],"\\.")[[1]][1]
  return(y)
}
read_induced_vlp<-function(x){
  pa<-paste(data.pa,"induced.virome/my.analysis/induced.vlp.to.stool.vlp",sep="")
  pa1<-paste(pa,"/",x,"/",sep="")
  z=list.files(pa1,patter=".txt",full.names = T)
  y<-lapply(z,read_idx)
  y<-do.call(rbind,y)
  return(y)
}  
data<-lapply(names,read_induced_vlp)
data<-do.call(rbind,data)
data<-data %>% separate(col=contig,into = c("induced","contig"))
tmp1<-data[which(substring(data$vlp,2)==substring(data$induced,3)),]
tmp2<-data[which(substring(data$vlp,2)!=substring(data$induced,3)),]
tmp1$group<-"Autologous"
tmp2$group<-"Heterologous"
data<-rbind(tmp1,tmp2)
data<-data %>% mutate(percentage=mapped/total.read)
tmp<-data[,c("group","percentage")]
tmp1<-aggregate(. ~ group,tmp, mean)
tmp2<-aggregate(. ~ group,tmp, function(x) {sd(x)/sqrt(length(x))}) # this is for standard erro
mtmp<-cbind(tmp1,tmp2[,2])
names(mtmp)<-c("Group","mean","se")
g1<-ggplot(mtmp,aes(Group,mean*100,width=0.5))+
  geom_bar(stat="identity", color="black",fill="grey",size=1,position=position_dodge())+
  geom_errorbar(aes(ymin=100*mean, ymax=100*mean+100*se), width=0.2,size=1,
                position=position_dodge(.9)) +
  scale_y_continuous(breaks=c(0,1,2),limits = c(0,2))+
  #scale_fill_manual(values=two.col)+
  ylab("Percentage (%)")+
  xlab("")+
  annotate("text", x=c(1.5), y=c(2),size=5, label= c(expression(paste(italic(P)," < 0.0001"))))+
  annotate("segment", x = c(1.05), xend = c(1.95), y = c(1.9), yend = c(1.9))+
  theme_classic()+
  theme(
    axis.title.y  = element_text(size=35,colour="black", margin =margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text = element_text(size=30,colour = "black"),
    legend.position = "none"
  )
ggsave("output/Fig.2F.pdf",width = 7,height =5 )
wilcox.test(data[data$group=="Autologous","percentage"],tmp[tmp$group=="Heterologous","percentage"])

# Fig.2G
tmp<-read.delim(paste(data.pa,"combined.txt",sep=""))
tmp$Inducer<-factor(tmp$Inducer,levels=c("Mitomycin C","Carbadox","Cipro","No inducer"))
tmp<- tmp %>% filter(uniprot.total.orf.number>0)
tmp<- tmp[!is.na(tmp$wgs.contig),]
tmp<-mutate(tmp,
            unip.p=uniprot.total.orf.number/orf.number,
            len1.p=alignment.len/ind.len,
            len2.p=alignment.len/wgs.len,
            ind.p=ind.phage.total/ind.total,
            vlp.p=vlp.reads.number/vlp.total,
            shot.p=shot.maptowgs.total/shotgun.total)
tmp<- tmp %>% filter(unip.p>=0.5,
                     len1.p>=0.5|len2.p>=0.5,
                     ind.p>=0.05)

ggplot(tmp[tmp$Inducer=="Mitomycin C",],aes(x = shot.p+10e-10,y = vlp.p+10e-10))+
  geom_point(colour="black",size=3,alpha=0.8)+
  geom_smooth(method='lm',formula=y~x,colour="black",linetype = "dashed")+
  scale_x_continuous(trans='log2',breaks=c(0.001,0.008,0.06,0.5),labels=c(0.1,0.8,6,50))+
  scale_y_continuous(trans='log2',breaks=c(10e-10,10e-8,10e-5,0.1),labels=c(0,expression("10"^-6),expression("10"^-3),expression("10"^1)))+
  annotate("text", x=0.001, y=0.2,size=10, label= expression(paste(italic(P)," = 0.0008")))+
  annotate("text", x=0.001, y=0.05,size=10, label= expression(paste(italic(R)," = 0.84")),hjust = 0.63)+
  labs(y="Induced phage abundance\nin purified VLP (%)",x="Isolated bacterial abundance (%)")+
  theme_classic()+
  theme( 
    axis.title = element_text(size=35,colour="black"),
    axis.text = element_text(size=30,colour = "black"),
    axis.title.y=element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
  )

ggsave("output/Fig.2G.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 9, height = 7, units = c("in"),
       dpi = 300, limitsize = TRUE)


tmp<-read.delim(paste(data.pa,"combined.txt",sep=""))
tmp$Inducer<-factor(tmp$Inducer,levels=c("Mitomycin C","Carbadox","Cipro","No inducer"))
tmp<- tmp %>% filter(uniprot.total.orf.number>0)
tmp<- tmp[!is.na(tmp$wgs.contig),]
tmp<-mutate(tmp,
            unip.p=uniprot.total.orf.number/orf.number,
            len1.p=alignment.len/ind.len,
            len2.p=alignment.len/wgs.len,
            ind.p=ind.phage.total/ind.total,
            vlp.p=vlp.reads.number/vlp.total,
            shot.p=shot.maptowgs.total/shotgun.total)
tmp<- tmp %>% filter(unip.p>=0.5,
                     len1.p>=0.5|len2.p>=0.5,
                     ind.p>=0.05)

ggplot(tmp[tmp$Inducer=="No inducer",],aes(x = shot.p+10e-10,y = vlp.p+10e-10))+
  geom_point(colour="black",size=3,alpha=0.8)+
  geom_smooth(method='lm',formula=y~x,colour="black",linetype = "dashed")+
  scale_x_continuous(trans='log2',breaks=c(0.001,0.008,0.06,0.5),labels=c(0.1,0.8,6,50))+
  scale_y_continuous(trans='log2',breaks=c(10e-10,10e-8,10e-5,0.1),labels=c(0,expression("10"^-6),expression("10"^-3),expression("10"^1)))+
  annotate("text", x=0.001, y=0.2,size=10, label= expression(paste(italic(P)," = 0.0004")))+
  annotate("text", x=0.001, y=0.05,size=10, label= expression(paste(italic(R)," = 0.83")),hjust = 0.63)+
  labs(y="Induced phage abundance\nin purified VLP (%)",x="Isolated bacterial abundance (%)")+
  theme_classic()+
  theme( 
    axis.title = element_text(size=35,colour="black"),
    axis.text = element_text(size=30,colour = "black"),
    axis.title.y=element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
  )

ggsave("output/Fig.S6.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 9, height = 7, units = c("in"),
       dpi = 300, limitsize = TRUE)

cor.test(tmp[tmp$Inducer=="No inducer",]$vlp.p,tmp[tmp$Inducer=="No inducer",]$shot.p,method = "spearman")


#### End ####


# rm(list=ls())
# options(warn=0)
