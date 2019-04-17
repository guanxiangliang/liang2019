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
library(RColorBrewer)
library(grid)
library(gtable) # need for gtable filter
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

#### Fig.1D viral richness ####
meta <- read.delim(paste(meta.pa,"meta.data.txt",sep=""))
kraken <- read.delim(paste(data.pa,"dis.vlp.kraken.tsv",sep=""), skip = 1)
names(kraken) <- gsub(x = names(kraken),
                      pattern = "\\.taxa",
                      replacement = "\\")
names(kraken) <- gsub(x = names(kraken),
                      pattern = "\\.",
                      replacement = "\\")
names(kraken)[1]<-"tax_id" ## kraken is table with 1 is tax id, 2:157 are libraries, 158 is names
lin <- read.delim("input/wet.input/dis.vlp.lin.txt") # lineages data by tax to baltimore
host <- read.delim("input/meta.data/viral.host.txt") # virus host data
kraken<-merge(kraken,lin,by="tax_id") # merge kraken and lineage data

# Prepare DNA and RNA virus table
virus <- kraken[kraken$superkingdom=="Viruses",] # only pick viruses 
virus <- virus[rowSums(is.na(virus))!=ncol(virus), ] # remove all NA in a row
virus<-merge(virus,host,by="family",all.x=T)

# Remove contamnation
virus<-virus[-which(virus$name %in% "Klebsiella virus 0507KN21"),]
virus<-virus[-which(virus$name %in% "Burkholderia phage ST79"),]
virus<-virus[-which(virus$name %in% "Pseudomonas phage PPpW-4"),]
virus<-virus[-which(virus$name %in% "Burkholderia phage KS9"),]
virus<-virus[-which(virus$name %in% "Bacillus virus phi29"),]


# DNA Virus species 
dna.sam<-meta %>% 
  filter(Cohort=="Discovery",study_group=="Month 0"|study_group=="Month 1"|study_group=="Month 4",library_type=="DNA") 
dna<-virus[,c(as.character(dna.sam$library_id))]
dna<-dna[rowSums(dna>=100)>=1,]
# RNA virus species
rna.sam<-meta %>% 
  filter(Cohort=="Discovery",study_group=="Month 0"|study_group=="Month 1"|study_group=="Month 4",library_type=="RNA") 
rna<-virus[,c(as.character(rna.sam$library_id))]
rna<-rna[rowSums(rna>=100)>=1,]

# To combine
names(dna)<-substring(names(dna), 2)
names(rna)<-substring(names(rna), 2)
com<-rbind(dna,rna)
names(com)<-paste0("D", names(com))
rich<-data.frame(librayr_id=names(com),colSums(com>=100))
names(rich)<-c("library_id","richness")

# merge with meta
df<-merge(rich,meta[,c("library_id","study_group")])
df$study_group<-factor(df$study_group, levels = c("Month 0", "Month 1", "Month 4"))
mtmp<-melt(df[,c("richness","study_group")])
g1<-ggplot(mtmp,aes(x=study_group,y=value+1)) +
  geom_violin(aes(fill=study_group))+
  geom_beeswarm(data=. %>% filter(value==0),cex=1,alpha=0.8,size=3,shape=21,colour="white",fill="black",groupOnX = T)+
  geom_beeswarm(data=. %>% filter(value>0),cex=4,alpha=0.8,size=3,shape=21,colour="white",fill="black",groupOnX = T)+
  scale_y_continuous(trans='log2', breaks=c(1,5,15,40),labels = c("0",expression(5),expression(15),expression(40)), expand = c(0.1,0))+
  scale_x_discrete(breaks=c("Month 0","Month 1","Month 4"),
                   labels=c("0", "1", "4"))+
  scale_fill_manual(values = age.col)+
  expand_limits(y = c(0, 50))+
  ylab("Viral richness")+
  xlab("Month")+
  annotate("text", x=c(1.5,2,2.5), y=c(40,58,40),size=5, label= c(expression(paste(italic(P)," = 0.0006")),expression(paste(italic(P)," = 0.0001")),"NS"))+
  annotate("segment", x = c(1.05,1.05,2.05), xend = c(1.95,2.95,2.95), y = c(35,50,35), yend = c(35,50,35))+
  theme_classic()+
  theme(axis.title = element_text(size=35,colour="black"),
        axis.text = element_text(size=30,colour = "black"),
        legend.position = "none"
  )
ggsave("output/Fig.1D.pdf",  width = 5, height = 5)

#### End ####

#### Fig.1E Kraken heatmap ####
meta <- read.delim(paste(meta.pa,"meta.data.txt",sep=""))
kraken <- read.delim(paste(data.pa,"dis.vlp.kraken.tsv",sep=""), skip = 1)
names(kraken) <- gsub(x = names(kraken),
                      pattern = "\\.taxa",
                      replacement = "\\")
names(kraken) <- gsub(x = names(kraken),
                      pattern = "\\.",
                      replacement = "\\")
names(kraken)[1]<-"tax_id" ## kraken is table with 1 is tax id, 2:157 are libraries, 158 is names
lin <- read.delim("input/wet.input/dis.vlp.lin.txt") # lineages data by tax to baltimore
host <- read.delim("input/meta.data/viral.host.txt") # virus host data
kraken<-merge(kraken,lin,by="tax_id") # merge kraken and lineage data

# Prepare DNA and RNA virus table
virus <- kraken[kraken$superkingdom=="Viruses",] # only pick viruses 
virus <- virus[rowSums(is.na(virus))!=ncol(virus), ] # remove all NA in a row
virus<-merge(virus,host,by="family",all.x=T)

# Remove contamnation
virus<-virus[-which(virus$name %in% "Klebsiella virus 0507KN21"),]
virus<-virus[-which(virus$name %in% "Burkholderia phage ST79"),]
virus<-virus[-which(virus$name %in% "Pseudomonas phage PPpW-4"),]
virus<-virus[-which(virus$name %in% "Burkholderia phage KS9"),]
virus<-virus[-which(virus$name %in% "Bacillus virus phi29"),]

# Family level
com.sam<- meta %>%
  filter(Cohort=="Discovery",study_group=="Month 0"|study_group=="Month 1"|study_group=="Month 4",library_type=="DNA"|library_type=="RNA")

family <- aggregate(.~family,virus[,c("family",as.character(com.sam$library_id))],sum)

#family<-family[,c(1,match(rownames(vlp.meta[vlp.meta$group2!="Control",]),colnames(family)))]


# DNA Virus species 
dna.sam<-meta %>% 
  filter(Cohort=="Discovery",study_group=="Month 0"|study_group=="Month 1"|study_group=="Month 4",library_type=="DNA") 
dna<-family[,c("family",as.character(dna.sam$library_id))]
dna<-dna[rowSums(dna[,-1]>=100)>=1,]
# RNA virus species
rna.sam<-meta %>% 
  filter(Cohort=="Discovery",study_group=="Month 0"|study_group=="Month 1"|study_group=="Month 4",library_type=="RNA") 
rna<-family[,c("family",as.character(rna.sam$library_id))]
rna<-rna[rowSums(rna[,-1]>=100)>=1,]
rna<-rna[-which(rna$family %in% c("Podoviridae","Siphoviridae")),] # remove contamination

names(dna)<-substring(names(dna), 2)
names(rna)<-substring(names(rna), 2)

com<-rbind(dna,rna)
names(com)<-paste0("D", names(com))

mcom<-melt(com)
names(mcom)<-c("family","library_id","reads")
mcom<-merge(mcom,host[,c("family","host.name")],by="family",all.x=T)
mcom<-merge(mcom,meta,by="library_id",all.x=T)
mcom$host.name<-factor(mcom$host.name,levels=c("Animal DNA Viruses","Animal RNA Viruses","Plant RNA Viruses","Bacteriophages"))

g1<- ggplot(mcom, aes(x = library_id, y = family, fill = log(reads+1,2))) +
  geom_tile(color="grey", size=0.1) +
  labs(x="",y="")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
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
  facet_grid(host.name~study_group,scales = "free", space = "free")+
  scale_fill_continuous(name="Read Number",na.value="white",guide=guide_colorbar(title.position = "top"),
                        limits = c(log(101,2),log(max(mcom$reads)+1,2)),
                        low=brewer.pal(9, name = 'OrRd')[2],high=brewer.pal(9, name = 'OrRd')[8],
                        breaks = c(log(101,2),log(1001,2),log(10001,2),log(100001,2),log(1000001,2) ),
                        labels=c(expression("10"^2),expression("10"^3),expression("10"^4),expression("10"^5),expression("10"^6)))+
  
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

pdf("output/Fig.1E.pdf", width = 15, height = 10)
grid.draw(g)
dev.off()
#### End ####

#### Fig.3A read coverage data and plot coverage heatmap####
# function to read files
read_bwa_coverage <- function(filename1,filename2) {
  p = read.delim(filename1, header=F,
                 stringsAsFactors = FALSE)
  names(p)<-c("family", 
              "read",
              "align.len","genome.len", "cov")
  p<-p[,c(1,2,5)]
  x<- p %>% filter(read>0)
  #n<-data.frame(virus_id=NA,cov=0,reads=0,total.read=0,sample_id=strsplit(strsplit(filename1,"\\/")[[1]][10],"\\.")[[1]][1])
  if (nrow(x) == 0) {
    x<-p[1,c("family","cov")]
    x$cov<-0
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
  z<- z %>% arrange(desc(cov))
  z<-z[1,c(1,2,4,5)]
  z$sample_id<-strsplit(strsplit(filename1,"\\/")[[1]][length(strsplit(filename1,"\\/")[[1]])],"\\.")[[1]][1]
  names(z)<-c("virus_id","cov","reads","total.read","library_id")
  return(z)
  
}
# function to make table 
table_bwa_coverage<- function(cohort,viralname){
  path1 = paste(data.pa,"coverage/",cohort,"/",viralname,sep="")
  filename.cov <- list.files(path = path1,full.names = TRUE,pattern = ".coverage")
  filename.reads <- list.files(path = path1,full.names = TRUE,pattern = ".read")
  bwa.cov.data <- mapply(read_bwa_coverage,filename.cov,filename.reads,SIMPLIFY = FALSE) # mapply do not give list default as lapply, need simplify option
  bwa.cov<-do.call(rbind,bwa.cov.data)
  rownames(bwa.cov)<-NULL
  return(bwa.cov)
}

# load meta
meta <- read.delim(paste(meta.pa,"meta.data.txt",sep=""))

# merge coverage and reads data
dis.adeno<-table_bwa_coverage("dis","adeno")
dis.adeno<-merge(dis.adeno,meta[,c("library_id","study_group","Infant_feeding_type","subject_id")],by="library_id")
dis.adeno$Cohort<-"Discovery"
dis.adeno$Family<-"Adenoviridae"
dis.anello<-table_bwa_coverage("dis","anello")
dis.anello<-merge(dis.anello,meta[,c("library_id","study_group","Infant_feeding_type","subject_id")],by="library_id")
dis.anello$Cohort<-"Discovery"
dis.anello$Family<-"Anelloviridae"
dis.calici<-table_bwa_coverage("dis","calici")
dis.calici<-merge(dis.calici,meta[,c("library_id","study_group","Infant_feeding_type","subject_id")],by="library_id")
dis.calici$Cohort<-"Discovery"
dis.calici$Family<-"Caliciviridae"
dis.parvo<-table_bwa_coverage("dis","parvo")
dis.parvo<-merge(dis.parvo,meta[,c("library_id","study_group","Infant_feeding_type","subject_id")],by="library_id")
dis.parvo$Cohort<-"Discovery"
dis.parvo$Family<-"Parvoviridae"
dis.picona<-table_bwa_coverage("dis","picona")
dis.picona<-merge(dis.picona,meta[,c("library_id","study_group","Infant_feeding_type","subject_id")],by="library_id")
dis.picona$Cohort<-"Discovery"
dis.picona$Family<-"Picornaviridae"
dis.polymavi<-table_bwa_coverage("dis","polymavi")
dis.polymavi<-merge(dis.polymavi,meta[,c("library_id","study_group","Infant_feeding_type","subject_id")],by="library_id")
dis.polymavi$Cohort<-"Discovery"
dis.polymavi$Family<-"Polyomaviridae"
dis.astro<-table_bwa_coverage("dis","astro")
dis.astro<-merge(dis.astro,meta[,c("library_id","study_group","Infant_feeding_type","subject_id")],by="library_id")
dis.astro$Cohort<-"Discovery"
dis.astro$Family<-"Astroviridae"

vali.adeno<-table_bwa_coverage("vali","adeno")
vali.adeno<-merge(vali.adeno,meta[,c("library_id","study_group","Infant_feeding_type","subject_id")],by="library_id")
vali.adeno$Cohort<-"Validation"
vali.adeno$Family<-"Adenoviridae"
vali.anello<-table_bwa_coverage("vali","anello")
vali.anello<-merge(vali.anello,meta[,c("library_id","study_group","Infant_feeding_type","subject_id")],by="library_id")
vali.anello$Cohort<-"Validation"
vali.anello$Family<-"Anelloviridae"
vali.calici<-table_bwa_coverage("vali","calici")
vali.calici<-merge(vali.calici,meta[,c("library_id","study_group","Infant_feeding_type","subject_id")],by="library_id")
vali.calici$Cohort<-"Validation"
vali.calici$Family<-"Caliciviridae"
vali.parvo<-table_bwa_coverage("vali","parvo")
vali.parvo<-merge(vali.parvo,meta[,c("library_id","study_group","Infant_feeding_type","subject_id")],by="library_id")
vali.parvo$Cohort<-"Validation"
vali.parvo$Family<-"Parvoviridae"
vali.picona<-table_bwa_coverage("vali","picona")
vali.picona<-merge(vali.picona,meta[,c("library_id","study_group","Infant_feeding_type","subject_id")],by="library_id")
vali.picona$Cohort<-"Validation"
vali.picona$Family<-"Picornaviridae"
vali.polymavi<-table_bwa_coverage("vali","polymavi")
vali.polymavi<-merge(vali.polymavi,meta[,c("library_id","study_group","Infant_feeding_type","subject_id")],by="library_id")
vali.polymavi$Cohort<-"Validation"
vali.polymavi$Family<-"Polyomaviridae"
vali.astro<-table_bwa_coverage("vali","astro")
vali.astro<-merge(vali.astro,meta[,c("library_id","study_group","Infant_feeding_type","subject_id")],by="library_id")
vali.astro$Cohort<-"Validation"
vali.astro$Family<-"Astroviridae"

# put all data together
l<-list(vali.adeno,vali.anello,vali.calici,vali.parvo,vali.picona,vali.polymavi,vali.astro,
        dis.adeno,dis.anello,dis.calici,dis.parvo,dis.picona,dis.polymavi,dis.astro)
tmp<-do.call(rbind,l)
tmp$inf<-"Uninfected"
tmp[tmp$cov>=(1/3),"inf"]<-"Infected"

# To make discovery cohort heatmap for Fig.3A


tmp$Infant_feeding_type<-gsub(tmp$Infant_feeding_type,pattern ="Mixed" ,replacement = "Breastmilk")
tmp$Infant_feeding_type<-gsub(tmp$Infant_feeding_type,pattern ="Breastmilk" ,replacement = "Breastmilk + Mixed")
tmp$Infant_feeding_type<-factor(tmp$Infant_feeding_type,levels=c("Formula","Breastmilk + Mixed"))
tmp$study_group<-factor(tmp$study_group, levels = c("Month 0", "Month 1", "Month 4","Negative control","Positive control"))
tmp<- tmp %>% filter(Cohort=="Discovery",Family!="Polyomaviridae",Family!="Astroviridae",Family!="Anelloviridae",
                     study_group=="Month 4"|study_group=="Month 1")

  
  
  g1<-ggplot(tmp,aes(x = subject_id, y = Family)) +
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
  facet_grid(.~study_group+Infant_feeding_type,scales = "free", space = "free")+
  theme(aspect.ratio = 1)


  z <- ggplotGrob(g1) 
  locations <- grep("strip-t", z$layout$name)
  strip <- gtable_filter(z, "strip-t", trim = FALSE)
  top <- strip$layout$t[1]
  l   <- strip$layout$l[c(1, 3)] # adjust left starting point
  r   <- strip$layout$r[c(2, 4)] # adjust right ending point
  mat   <- matrix(vector("list", length = 6), nrow = 2)
  mat[] <- list(zeroGrob())
  
  # The separator for the facets has zero width
  res <- gtable_matrix("toprow", mat, unit(c(1, 0, 1), "null"), unit(c(1, 1), "null"))
  zz <- res %>%
    gtable_add_grob(z$grobs[[locations[1]]]$grobs[[1]], 1, 1, 1, 3) %>%
    gtable_add_grob(z, ., t = top,  l = l[1],  b = top,  r = r[1], name = c("add-strip"))
  
  
  # Adding the second layer (note the indices) locations[4] indicate the label you need 
  pp <- gtable_add_grob(res, z$grobs[[locations[4]]]$grobs[[1]], 1, 1, 1, 3) %>%
    gtable_add_grob(zz, ., t = top,  l = l[2],  b = top,  r = r[2], name = c("add-strip"))
  
  # add color to bigger grip
  stript <- which(grepl('add-strip', pp$layout$name))
  fills <- age.col[2:3] # give the age color 
  k <- 1
  for (i in stript) {
    j <- which(grepl('rect', pp$grobs[[i]]$grobs[[7]]$childrenOrder)) # this is to find the background
    pp$grobs[[i]]$grobs[[7]]$children[[j]]$gp$fill <- fills[k] # this is for the top grip
    pp$grobs[[i]]$grobs[[7]]$children[[j]]$gp$lwd <-0
    k <- k+1
  }
  
  # add color to small grip
  stript <- which(grepl('strip-t', pp$layout$name))
  fills <- c(feed.col[1:2],feed.col[1:2]) # give the age color 
  k <- 1
  for (i in stript) {
    j <- which(grepl('rect', pp$grobs[[i]]$grobs[[2]]$childrenOrder))
    pp$grobs[[i]]$grobs[[2]]$children[[j]]$gp$fill <- fills[k] # this is for the second grip
    pp$grobs[[i]]$grobs[[2]]$children[[j]]$gp$lwd <- 0 # remove outline
    j<-which(grepl('title', pp$grobs[[i]]$grobs[[2]]$childrenOrder)) # this is to find the tile
    pp$grobs[[i]]$grobs[[2]]$children[[j]]$children[[1]]$gp$fontsize<-10 # change the strip text size as 15
    k <- k+1
  }

pdf("output/Fig.3A.pdf", width = 12, height = 10)
grid.draw(pp)
dev.off()






#### End ####

#### Fig.3B validation plot coverage heatmap####
# put all data together
l<-list(vali.adeno,vali.anello,vali.calici,vali.parvo,vali.picona,vali.polymavi,vali.astro,
        dis.adeno,dis.anello,dis.calici,dis.parvo,dis.picona,dis.polymavi,dis.astro)
tmp<-do.call(rbind,l)
tmp$inf<-"Uninfected"
tmp[tmp$cov>=(1/3),"inf"]<-"Infected"

# To make discovery cohort heatmap for Fig.3A


tmp$Infant_feeding_type<-gsub(tmp$Infant_feeding_type,pattern ="Mixed" ,replacement = "Breastmilk")
tmp$Infant_feeding_type<-gsub(tmp$Infant_feeding_type,pattern ="Breastmilk" ,replacement = "Breastmilk + Mixed")
tmp$Infant_feeding_type<-factor(tmp$Infant_feeding_type,levels=c("Formula","Breastmilk + Mixed"))
tmp$study_group<-factor(tmp$study_group, levels = c("Month 0", "Month 1", "Month 4","Negative control","Positive control"))
tmp<- tmp %>% filter(Cohort=="Validation",Family!="Polyomaviridae",
                     study_group=="Month 4"|study_group=="Month 1")



g1<-ggplot(tmp,aes(x = subject_id, y = Family)) +
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
  facet_grid(.~study_group+Infant_feeding_type,scales = "free", space = "free")+
  theme(aspect.ratio = 1)

z <- ggplotGrob(g1) 
locations <- grep("strip-t", z$layout$name)
strip <- gtable_filter(z, "strip-t", trim = FALSE)
top <- strip$layout$t[1]
l   <- strip$layout$l[c(1)] # adjust left starting point
r   <- strip$layout$r[c(2)] # adjust right ending point
mat   <- matrix(vector("list", length = 6), nrow = 2)
mat[] <- list(zeroGrob())

# The separator for the facets has zero width
res <- gtable_matrix("toprow", mat, unit(c(1, 0, 1), "null"), unit(c(1, 1), "null"))
zz <- res %>%
  gtable_add_grob(z$grobs[[locations[1]]]$grobs[[1]], 1, 1, 1, 3) %>%
  gtable_add_grob(z, ., t = top,  l = l[1],  b = top,  r = r[1], name = c("add-strip")) # add the bigger layer on the top

pp<-zz
# Adding the second layer (note the indices) locations[4] indicate the label you need 
# pp <- gtable_add_grob(res, z$grobs[[locations[2]]]$grobs[[1]], 1, 1, 1, 3) %>%
#   gtable_add_grob(zz, ., t = top,  l = l[1],  b = top,  r = r[1], name = c("add-strip"))

# add color to bigger grip
stript <- which(grepl('add-strip', pp$layout$name))
fills <- age.col[3] # give the age color 
k <- 1
for (i in stript) {
  j <- which(grepl('rect', pp$grobs[[i]]$grobs[[7]]$childrenOrder)) # this is to find the background
  pp$grobs[[i]]$grobs[[7]]$children[[j]]$gp$fill <- fills[k] # this is for the top grip
  pp$grobs[[i]]$grobs[[7]]$children[[j]]$gp$lwd <-0
  k <- k+1
}

# add color to small grip
stript <- which(grepl('strip-t', pp$layout$name))
fills <- feed.col[1:2] # give the age color 
k <- 1
for (i in stript) {
  j <- which(grepl('rect', pp$grobs[[i]]$grobs[[2]]$childrenOrder))
  pp$grobs[[i]]$grobs[[2]]$children[[j]]$gp$fill <- fills[k] # this is for the second grip
  pp$grobs[[i]]$grobs[[2]]$children[[j]]$gp$lwd <- 0 # remove outline
  j<-which(grepl('title', pp$grobs[[i]]$grobs[[2]]$childrenOrder)) # this is to find the tile
  pp$grobs[[i]]$grobs[[2]]$children[[j]]$children[[1]]$gp$fontsize<-15 # change the strip text size as 15
  k <- k+1
}

pdf("output/Fig.3B.pdf", width = 12, height = 10)
grid.draw(pp)
dev.off()






#### End ####

#### Fig.3C Infection percentage ####
# put all data together
l<-list(vali.adeno,vali.anello,vali.calici,vali.parvo,vali.picona,vali.polymavi,vali.astro,
        dis.adeno,dis.anello,dis.calici,dis.parvo,dis.picona,dis.polymavi,dis.astro)
tmp<-do.call(rbind,l)
tmp$inf<-"Uninfected"
tmp[tmp$cov>=(1/3),"inf"]<-"Infected"
tmp$Infant_feeding_type<-gsub(tmp$Infant_feeding_type,pattern ="Mixed" ,replacement = "Breastmilk")
tmp$Infant_feeding_type<-gsub(tmp$Infant_feeding_type,pattern ="Breastmilk" ,replacement = "Breastmilk")

tmp1<-tmp %>% filter(study_group=="Month 4")
tmp1<-tmp1[,c("subject_id","inf","Infant_feeding_type","Cohort")]
tmp2<-tmp1[tmp1$inf=="Infected",] %>% distinct(subject_id,.keep_all= TRUE)
tmp1<-tmp1 %>% distinct(subject_id,.keep_all= TRUE)
tmp1<-tmp1[-which(tmp1$subject_id %in% unique(tmp2$subject_id)),]
tmp<-rbind(tmp1,tmp2)

tmp<-as.data.frame(tmp %>% count(inf,Infant_feeding_type,Cohort))
tmp1<-tmp%>%filter(Infant_feeding_type=="Formula")
tmp2<-tmp%>%filter(Infant_feeding_type=="Breastmilk")
tmp2<-rbind(tmp2,data.frame(inf="Infected",Infant_feeding_type="Breastmilk",Cohort="Discovery",n=0))

tmp1<- tmp1%>% filter(Cohort=="Discovery")%>% mutate(freq=n/sum(n))
tmp2<- tmp2%>% filter(Cohort=="Discovery")%>%mutate(freq=n/sum(n))
tmp<-rbind(tmp1,tmp2)
ci.f<-binom.test(tmp[tmp$Infant_feeding_type=="Formula","n"][1],sum(tmp[tmp$Infant_feeding_type=="Formula","n"]))$conf[2] # calvculate 95% confidence interval proportion
ci.b<-binom.test(tmp[tmp$Infant_feeding_type=="Breastmilk","n"][1],sum(tmp[tmp$Infant_feeding_type=="Breastmilk","n"]))$conf[2] # calvculate 95% confidence interval proportion
cil.f<-binom.test(tmp[tmp$Infant_feeding_type=="Formula","n"][1],sum(tmp[tmp$Infant_feeding_type=="Formula","n"]))$conf[1] # calvculate 95% confidence interval proportion
cil.b<-binom.test(tmp[tmp$Infant_feeding_type=="Breastmilk","n"][1],sum(tmp[tmp$Infant_feeding_type=="Breastmilk","n"]))$conf[1] # calvculate 95% confidence interval proportion
data1<-tmp %>% filter(inf=="Infected")
data1$ci<-NA
data1[data1$Infant_feeding_type=="Formula",]$ci<-ci.f
data1[data1$Infant_feeding_type=="Breastmilk",]$ci<-ci.b
data1$cil<-NA
data1[data1$Infant_feeding_type=="Formula",]$cil<-cil.f
data1[data1$Infant_feeding_type=="Breastmilk",]$cil<-cil.b
data1$Infant_feeding_type<-factor(data1$Infant_feeding_type,levels=c("Formula","Breastmilk"))




l<-list(vali.adeno,vali.anello,vali.calici,vali.parvo,vali.picona,vali.polymavi,vali.astro,
        dis.adeno,dis.anello,dis.calici,dis.parvo,dis.picona,dis.polymavi,dis.astro)
tmp<-do.call(rbind,l)
tmp$inf<-"Uninfected"
tmp[tmp$cov>=(1/3),"inf"]<-"Infected"
tmp$Infant_feeding_type<-gsub(tmp$Infant_feeding_type,pattern ="Mixed" ,replacement = "Breastmilk")
tmp$Infant_feeding_type<-gsub(tmp$Infant_feeding_type,pattern ="Breastmilk" ,replacement = "Breastmilk")

tmp1<-tmp %>% filter(study_group=="Month 4")
tmp1<-tmp1[,c("subject_id","inf","Infant_feeding_type","Cohort")]
tmp2<-tmp1[tmp1$inf=="Infected",] %>% distinct(subject_id,.keep_all= TRUE)
tmp1<-tmp1 %>% distinct(subject_id,.keep_all= TRUE)
tmp1<-tmp1[-which(tmp1$subject_id %in% unique(tmp2$subject_id)),]
tmp<-rbind(tmp1,tmp2)

tmp<-as.data.frame(tmp %>% count(inf,Infant_feeding_type,Cohort))
tmp1<-tmp%>%filter(Infant_feeding_type=="Formula")
tmp2<-tmp%>%filter(Infant_feeding_type=="Breastmilk")

tmp1<- tmp1%>% filter(Cohort=="Validation")%>% mutate(freq=n/sum(n))
tmp2<- tmp2%>% filter(Cohort=="Validation")%>%mutate(freq=n/sum(n))
tmp<-rbind(tmp1,tmp2)
ci.f<-binom.test(tmp[tmp$Infant_feeding_type=="Formula","n"][1],sum(tmp[tmp$Infant_feeding_type=="Formula","n"]))$conf[2] # calvculate 95% confidence interval proportion
ci.b<-binom.test(tmp[tmp$Infant_feeding_type=="Breastmilk","n"][1],sum(tmp[tmp$Infant_feeding_type=="Breastmilk","n"]))$conf[2] # calvculate 95% confidence interval proportion
cil.f<-binom.test(tmp[tmp$Infant_feeding_type=="Formula","n"][1],sum(tmp[tmp$Infant_feeding_type=="Formula","n"]))$conf[1] # calvculate 95% confidence interval proportion
cil.b<-binom.test(tmp[tmp$Infant_feeding_type=="Breastmilk","n"][1],sum(tmp[tmp$Infant_feeding_type=="Breastmilk","n"]))$conf[1] # calvculate 95% confidence interval proportion
data2<-tmp %>% filter(inf=="Infected")
data2$ci<-NA
data2[data2$Infant_feeding_type=="Formula",]$ci<-ci.f
data2[data2$Infant_feeding_type=="Breastmilk",]$ci<-ci.b
data2$cil<-NA
data2[data2$Infant_feeding_type=="Formula",]$cil<-cil.f
data2[data2$Infant_feeding_type=="Breastmilk",]$cil<-cil.b
data2$Infant_feeding_type<-factor(data2$Infant_feeding_type,levels=c("Formula","Breastmilk"))

data<-rbind(data1,data2)
data[2,6:7]<-0

g1<-ggplot(data,aes(y = freq*100, x = Infant_feeding_type,fill=Infant_feeding_type)) +
  
  geom_bar( color="black",stat="identity",width = 0.5)+
  geom_errorbar(aes(ymin=cil*100, ymax=ci*100), width=0.2,size=0.5,
                position=position_dodge(.9)) +
  scale_fill_manual(values=feed.col[1:2])+
 #expand_limits(y = c(0, 37))+
  labs(title="",x ="", y = "Percentage of subjects (%)\nwith human-cell viruses")+
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
g <- ggplot_gtable(ggplot_build(g1))
stript <- which(grepl('strip-t', g$layout$name))
fills <- cohort.col
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

pdf("output/Fig.3C.pdf", width = 8, height = 6)
grid.draw(g)
dev.off()




#### End ####






# rm(list=ls())
# options(warn=0)










