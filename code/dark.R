##usage: dark(sampleID)
#library(tidyr) 
#library(taxonomizr)### convert accession id to tax id and to names ###
#define the commonly used path
#### Set up working directory ####
data.pa<-"./input/wet.input/"
meta.pa<-"./input/meta.data/"
database.pa<-"./input/database/"
options(warn=-1)
#### End ####
uniprottax<-read.delim(paste(database.pa,"uniprot.virome/uniprot.virome.meta",sep=""),header=F)
names(uniprottax)<-c("uniprot.hit.id","taxid","name")
uniprottax$taxid<-as.numeric(as.character(uniprottax$taxid))
hostd<-read.delim(paste(meta.pa,"viral.host.txt",sep=""))
names(hostd)<-c("uniprot.family","host","host.name")
dark.pa<-paste(data.pa,"discovery.virome/my.analysis/",sep="")

dark<-function (x){
  #load contig reads number (decontam read map contig)
  datapath <- file.path(paste(dark.pa,"contigs/",x,"/",x,".contig.reads.txt",sep=""))
  map <- read.delim(datapath,sep="\t",header=F)
  names(map)<-c("contig","contig.length","vlp.reads.number","unmmaped.reads.number")
  data<-map[map$contig.length>=3000,]
  data<-mutate(data, fpkm=vlp.reads.number*10e9/contig.length/sum(data$vlp.reads.number))
  
  #load contig coverage data (decontam read map contig)
  datapath <- file.path(paste(dark.pa,"contigs/",x,"/",x,".contig.coverage.txt",sep=""))
  cov <- read.delim(datapath,sep="\t",header=F)
  names(cov)<-c("contig","vlp.coverage")
  data<-merge(data,cov,all.x=T,by="contig")
  
  #load orf number data
  datapath <- file.path(paste(dark.pa,"orfs/",x,"/",x,".orf.number",sep=""))
  orfn<-read.delim(datapath,sep="\t",header=F)
  orfn<-separate(data = orfn, col = V1, into = c("contig"), sep = "\\_")
  tmp<-data.frame(table(orfn$contig))
  names(tmp)<-c("contig","orf.number")
  data<-merge(data,tmp,all.x=T,by="contig")
  
  #load nt blast result
  datapath <- file.path(paste(dark.pa,"blast.data/",x,"/",x,".3k.nt",sep=""))
  nt <- read.delim(datapath,header=F)
  names(nt)<-c("contig", "sseqid", "pident", "qlen","slen","alignment.len" ,"mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  nt<-nt[nt$evalue<1e-10,] 
  nt<-nt[!duplicated(nt$contig),][,c("contig","sseqid")]
  names(nt)<-c("contig","nt.top.hit")
  nt$nt.top.hit<-separate(data = nt, col = nt.top.hit, into = c("1", "2","3","ntname"), 
                   sep = "\\|")$ntname
  taxaId<-accessionToTaxa(nt$nt.top.hit,"./input/database/ac2tax/accessionTaxa.sql")
  nt<-mutate(nt, nt.taxa.kingdom=getTaxonomy(taxaId,"./input/database/ac2tax/accessionTaxa.sql")[,1],nt.taxa.speices=getTaxonomy(taxaId,"./input/database/ac2tax/accessionTaxa.sql")[,7])
  data<-merge(data,nt,all.x=T,by="contig")
  
  #load viral genome blast 
  datapath <- file.path(paste(dark.pa,"blast.data/",x,"/",x,".3k.viral.genome",sep=""))
  viral <- read.delim(datapath,header=F)
  names(viral)<-c("contig", "sseqid", "pident", "qlen","slen","alignment.len" ,"mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  viral<-viral[viral$evalue<0.00001,] 
  viral<-viral[!duplicated(viral$contig),][,c("contig","sseqid")]
  names(viral)<-c("contig","viral.top.hit")
  viral$viral.top.hit<-separate(data = viral, col = viral.top.hit, into = c("1", "viral.name"), 
                          sep = "\\|")$viral.name
  taxaId<-accessionToTaxa(viral$viral.top.hit,"./input/database/ac2tax/accessionTaxa.sql")
  viral<-mutate(viral, viral.taxa.family=getTaxonomy(taxaId,"./input/database/ac2tax/accessionTaxa.sql")[,5],viral.taxa.species=getTaxonomy(taxaId,"./input/database/ac2tax/accessionTaxa.sql")[,7])
  data<-merge(data,viral,all.x=T,by="contig")
  
  # #load contig mapping results
  # datapath <- file.path(paste(dark.pa,"orfs/",x,"/",x,".shot.contig.blast",sep=""))
  # contigshot <- read.delim(datapath,header=F)
  # names(contigshot)<-c("contig", "shotcontig", "pident", "qlen","slen","alignment.len" ,"mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  # contigshot<-contigshot[(contigshot$alignment.len/contigshot$qlen>0.9|contigshot$alignment.len/contigshot$slen>0.9)&(contigshot$pident>0.9),]
  # contigshot<-contigshot[!duplicated(contigshot$contig),]
  # tmp<-contigshot[,c("contig","shotcontig")]
  # data<-merge(data,tmp,all.x=T,by="contig")
  
  
  #load bacterial genome blast result
  datapath <- file.path(paste(dark.pa,"orfs/",x,"/",x,".contig.3k.bacteria",sep=""))
  bac <- read.delim(datapath,header=F)
  names(bac)<-c("contig", "sseqid", "pident", "qlen","slen","alignment.len" ,"mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  bac<-bac[(bac$alignment.len/bac$qlen>0.9)&(bac$pident>0.9),] 
  bac<-bac[!duplicated(bac$contig),]
  tmp<-data.frame(table(bac$contig))
  names(tmp)<-c("contig","bacterial.genome")
  data<-merge(data,tmp,all.x=T,by="contig")
  data$bacterial.genome[is.na(data$bacterial.genome)]<-0
  
  #load integrase data
  datapath <- file.path(paste(dark.pa,"orfs/",x,"/",x,".3k.integrase",sep=""))
  inte<-read.delim(datapath,header=T)
  names(inte)<-c("contig", "sseqid", "pident", "qlen","slen","alignment.len" ,"mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  inte<-separate(data = inte, col = contig, into = c("contig"), sep = "\\_")
  tmp<-data.frame(table(inte$contig))
  names(tmp)<-c("contig","integrase")
  data<-merge(data,tmp,all.x=T,by="contig")
  data$integrase[is.na(data$integrase)]<-0
  
  #load aclme data
  datapath <- file.path(paste(dark.pa,"orfs/",x,"/",x,".3k.aclame",sep=""))
  acl<-read.delim(datapath,header=T)
  names(acl)<-c("contig", "sseqid", "pident", "qlen","slen","alignment.len" ,"mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  acl<-separate(data = acl, col = contig, into = c("contig"), sep = "\\_")
  tmp<-data.frame(table(acl$contig))
  names(tmp)<-c("contig","prophage.protein")
  data<-merge(data,tmp,all.x=T,by="contig")
  data$prophage.protein[is.na(data$prophage.protein)]<-0
  
  
  #load negative control reads data
  datapath <- file.path(paste(dark.pa,"nc.coverage/",x,"/",x,".nc.reads.txt",sep=""))
  tmp<-read.delim(datapath,sep="\t",header=F)[,c("V1","V3")]
  names(tmp)<-c("contig","nc.reads.number")
  data<-merge(data,tmp,all.x=T,by="contig")
  
  #load negative control coverage data
  datapath <- file.path(paste(dark.pa,"nc.coverage/",x,"/",x,".nc.coverage.txt",sep=""))
  tmp<-read.delim(datapath,sep="\t",header=F)
  names(tmp)<-c("contig","nc.coverage")
  data<-merge(data,tmp,all.x=T,by="contig")
  data$nc.coverage[is.na(data$nc.coverage)]<-0
  
  #load uniprot data
  datapath <- file.path(paste(dark.pa,"orfs/",x,"/",x,".3k.uniprot",sep=""))
  uniprot<-read.delim(datapath,header=F,sep="\t")
  names(uniprot)<-c("contig", "sseqid", "pident", "qlen","slen","alignment.len" ,"mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  uniprot<-separate(data = uniprot, col = contig, into = c("contig","orfs.id"), sep = "\\_")
  uniprot<-separate(data = uniprot, col = sseqid, into = c("1","2","uniprot.hit.id"), sep = "\\|")
  uniprot<-merge(uniprot,uniprottax,by="uniprot.hit.id")
  uniprot<-cbind(uniprot,getTaxonomy(uniprot$taxid,"./input/database/ac2tax/accessionTaxa.sql"))
  tmp<-data.frame(table(uniprot$contig))
  names(tmp)<-c("contig","uniprot.total.orf.number")
  data<-merge(data,tmp,all.x=T,by="contig")
  data$uniprot.total.orf.number[is.na(data$uniprot.total.orf.number)]<-0
  
  #calculate uniprot family 
  tmp<-uniprot[,c("contig","family")]
  tmp<-table(tmp)
  tmp1<-data.frame(colnames(tmp)[apply(tmp,1,which.max)])
  tmp<-uniprot[,c("contig","family")]
  tmp<-data.frame(table(tmp))
  tmp2<-data.frame(tapply(tmp$Freq, tmp$contig, max))
  tmp3<-cbind(tmp1,tmp2)
  names(tmp3)<-c("uniprot.family","uniprot.family.orf.number")
  tmp3$contig<-rownames(tmp3)
  rownames(tmp3)<-NULL
  data<-merge(data,tmp3,all.x=T,by="contig")
  data$uniprot.family.orf.number[is.na(data$uniprot.family.orf.number)]<-0
  #calculate uniprot species
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
  data<-merge(data,tmp3,all.x=T,by="contig")
  data$uniprot.species.orf.number[is.na(data$uniprot.species.orf.number)]<-0
  

  #catergrize virues
  
  data = within(data, {
    source = ifelse((nc.coverage >= 0.5)|(nc.reads.number>=100), "Contamination", ifelse((uniprot.total.orf.number/orf.number>=0.5)&(uniprot.total.orf.number/contig.length*10000>=1),
                                                                      "Viruses",ifelse(!is.na(viral.top.hit),"Viruses",nt.taxa.kingdom
                                                                      )))})
  data$source[is.na(data$source)]<-"Unassigned"
  
  # Predict name
  data$uniprot.species<-as.character(data$uniprot.species)
  data = within(data, {
    name = ifelse(source!="Viruses", "NoVir", ifelse(!is.na(viral.taxa.species),viral.taxa.species,
                                                     ifelse(uniprot.species.orf.number>=2,uniprot.species,"Others")
    ))})
  
  
  #  data = within(data, {
  #    source = ifelse(nc.vlp.cov.5reads >= 0.7, "Contamination", ifelse((uniprot.total.orf.number/orf.number>=0.25)&(uniprot.total.orf.number/contig.length*10000>=minimal.orf),
  #                                                                      "Viruses",nt.taxa.kingdom
  #                                                                      ))})
  #
  
  data<-merge(data,hostd[,c("uniprot.family","host")],all.x=T,by="uniprot.family")
  data= within(data, {
    viral.type=ifelse((prophage.protein >=2|integrase>0|bacterial.genome>0),"Temperate Phages",
                      ifelse((uniprot.total.orf.number<2)|(uniprot.family.orf.number<2),"Uncharacterized Viruses",
                             ifelse(as.character(data$host)=="Bacteriophages","Other phages",as.character(data$host))))
  })
  
  
  
  # data= within(data, {
  #   viral.type=ifelse((prophage.protein >=2|integrase>0|bacterial.genome>0),"Temperate Phages",
  #   ifelse(as.character(data$host)=="Bacteriophages","Other phages",as.character(data$host)))
  # })
  #data= within(data, {
  #  induction=ifelse(ind.vlp.number>=5,"Induction","Others")
  #})
  
  ##load head protein data
  #datapath <- file.path(paste(pa,"head/",x,"/","head_best_hits.csv",sep=""))
  #head<-read.delim(datapath,sep=",",header=T)
  #head<-separate(data = head, col = query.def, into = c("contig"), sep = "\\ ")
  #head<-separate(data = head, col = contig, into = c("contig"), sep = "\\_")
  #tmp<-data.frame(table(head$contig))
  #names(tmp)<-c("contig","phage.head")
  #data<-merge(data,tmp,all.x=T,by="contig")
  #data$phage.head[is.na(data$phage.head)]<-0
  ##load tail protein data
  #datapath <- file.path(paste(pa,"tail/",x,"/","tail_best_hits.csv",sep=""))
  #tail<-read.delim(datapath,sep=",",header=T)
  #tail<-separate(data = tail, col = query.def, into = c("contig"), sep = "\\ ")
  #tail<-separate(data = tail, col = contig, into = c("contig"), sep = "\\_")
  #tmp<-data.frame(table(tail$contig))
  #names(tmp)<-c("contig","phage.tail")
  #data<-merge(data,tmp,all.x=T,by="contig")
  #data$phage.tail[is.na(data$phage.tail)]<-0
  
  ## combine coverage data
  #datapath <- file.path(paste(pa,"total.vlp/",x,"/",x,".read.coverage.txt",sep=""))
  #tmp1<-read.delim(datapath,sep="\t",header=F)
  #names(tmp1)<-c("contig","vlp.cov.5reads")
  #datapath <- file.path(paste(pa,"total.vlp/",x,"/nc.",x,".read.number.txt",sep=""))
  #tmp2<-read.delim(datapath,sep="\t",header=F)[,c(1,3)]
  #names(tmp2)<-c("contig","nc.vlp.number")
  #datapath <- file.path(paste(pa,"total.vlp/",x,"/nc.",x,".read.coverage.txt",sep=""))
  #tmp3<-read.delim(datapath,sep="\t",header=F)
  #names(tmp3)<-c("contig","nc.vlp.cov.5reads")
  #datapath <- file.path(paste(pa,"total.vlp/",x,"/ind.",x,".read.number.txt",sep=""))
  #tmp4<-read.delim(datapath,sep="\t",header=F)[,c(1,3)]
  #names(tmp4)<-c("contig","ind.vlp.number")
  #datapath <- file.path(paste(pa,"total.vlp/",x,"/ind.",x,".read.coverage.txt",sep=""))
  #tmp5<-read.delim(datapath,sep="\t",header=F)
  #names(tmp5)<-c("contig","ind.vlp.cov.5reads")
  ## Calculate EVENESS
  #datapath <- file.path(paste(pa,"total.vlp/",x,"/",x,".full.coverage",sep=""))
  #tmp6<-read.delim(datapath,sep="\t",header=F)
  #even<-matrix(rep("NA",2*length(levels(tmp6$V1))), nrow=length(levels(tmp6$V1)), ncol=2)
  #j=1
  #for (i in levels(tmp6$V1)){
  #  D=tmp6[tmp6$V1==i,"V3"]
  #  C=round(mean(D))
  #  D2=D[D<=C]
  #  E=1-(length(D2)-sum(D2)/C)/length(D)
  #  even[j,1]=i
  #  even[j,2]=E
  #  j=j+1
  #}
  #even<-print.data.frame(data.frame(even), 
  #                       quote=FALSE)
  #names(even)<-c("contig","eveness")
  ##com<-Reduce(function(x, y) merge(x, y, all=TRUE,by="contig"), list(tmp1,tmp2,tmp3, tmp4,tmp5,even)) # merge multiple dataframe
  #com<-Reduce(function(x, y) merge(x, y, all=TRUE,by="contig"), list(tmp1,tmp2,tmp3, tmp4,tmp5)) # merge multiple dataframe
  #if("nc.vlp.cov.5reads" !%in% colnames(com))
  #{
  #  com$nc.vlp.cov.5reads<-NA
  #}
  #com[is.na(com)]<-0
  #data<-merge(data,com,all.x=T,by="contig")
  
  return(data)
}

