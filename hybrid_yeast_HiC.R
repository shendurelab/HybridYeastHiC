# hybrid_yeast_HiC.R
# plotting of HiC data
# Seungsoo Kim

# load libraries ----
library(ggplot2)
library(dplyr)
library(reshape)
library(RColorBrewer)
library(data.table)
library(grid)
library(gplots)
library(scales)
library(permute)
library(gridExtra)
library(Cairo)

# output directory ----
dir.create("nobackup/figures")
outdir <- "nobackup/figures"
date <- "2016-12-05"

# color palette
brewercols <- brewer.pal(4,"Set1")
heatmapcols <- c("white","orange","red","black")
heatmaprange <- c(0,3)

# plot themes ----

# font size 7
paper.fig <- theme(panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank(),
                   text=element_text(size=7),
                   axis.text=element_text(size=7))
# with rotated y-axis labels
paper.fig.rot <- theme(panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  text=element_text(size=7),
                  axis.text=element_text(size=7),
                  axis.text.y=element_text(angle=90,hjust=0.5))
# heat maps
paperhm <- theme(axis.line = element_blank(),
                 axis.title = element_blank(), 
                 legend.key.size = unit(0.5,"cm"), 
                 legend.title = element_text(size=7), 
                 legend.text = element_text(size=7))
# heat maps 2
hm <- theme(axis.line = element_blank(), 
            axis.title = element_blank(), 
            legend.key.size = unit(1,"cm"), 
            legend.title = element_text(size=12), 
            legend.text = element_text(size=12))

# chromosome labels ----
chrlabs.all <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI")
chrlabs.odd <- c("I","","III","","V","","VII","","IX","","XI","","XIII","","XV","")

# mapping rates ----

# laod data
mapping_ScSu <- data.frame(hybrid="S. cerevisiae x S. uvarum",
                           maprate=read.table("nobackup/mappability/ScSu_w10_r80_q30.map_rate",header=FALSE))
mapping_ScSp <- data.frame(hybrid="S. cerevisiae x S. paradoxus",
                           maprate=read.table("nobackup/mappability/ScSp_w10_r80_q30.map_rate",header=FALSE))
mapping_ScSc <- data.frame(hybrid="S. cerevisiae x S. cerevisiae",
                           maprate=read.table("nobackup/mappability/Y12xDBVPG6044_w10_r150_q30.map_rate",header=FALSE))
mapping <- do.call("rbind", list(mapping_ScSu, mapping_ScSp, mapping_ScSc))
colnames(mapping) <- c("hybrid","maprate")

# labels
hybrids <- c("S. cerevisiae x S. uvarum","S. cerevisiae x S. paradoxus", "S. cerevisiae x S. cerevisiae")
hyblabs <- c(expression(paste(italic("S. cerevisiae"), " x ", italic("S. uvarum"))),
             expression(paste(italic("S. cerevisiae"), " x ", italic("S. paradoxus"))),
             expression(paste(italic("S. cerevisiae"), " x ", italic("S. cerevisiae"))))

# save figure with labels
pdf(paste(outdir,"/",date,"_map_rates.pdf",sep=""),5,3)
ggplot(mapping) + geom_bar(aes(x=factor(hybrid,levels=hybrids),y=maprate),stat="identity") + 
  theme_classic() + scale_y_continuous(limits=c(0,1),expand=c(0,0)) + 
  xlab("") + ylab(expression("Proportion of simulated reads mapping with MAPQ" >= 30)) + 
  theme(text=element_text(size=7)) + scale_x_discrete(labels=hyblabs)
dev.off()

# save figure without labels
pdf(paste(outdir,"/",date,"_map_rates_nolab.pdf",sep=""),3.5,3)
ggplot(mapping) + geom_bar(aes(x=factor(hybrid,levels=hybrids),y=maprate),stat="identity") + theme_classic() + scale_y_continuous(limits=c(0,1),expand=c(0,0)) + xlab("") + ylab(expression("Proportion of simulated reads mapping with MAPQ" >= 30)) + theme(text=element_text(size=7),axis.text.x=element_blank())
dev.off()

# homologous pairing schematic ----

plot_refs <- read.table("plotting/plotting_refs.txt",stringsAsFactors=FALSE,sep = "\t")
colnames(plot_refs) <- c("ref","g1","g2","g1p","g2p","rDNA1","rDNA2","bin1","bin2","c13b1","c13b2")

samp <- "ILY456_galactose_Sau3AI"
ref <- "ScSu"
bsize <- 32000
mask <- "Surep-rDNA"

# get reference information
g1 <- plot_refs[plot_refs$ref==ref,]$g1
g2 <- plot_refs[plot_refs$ref==ref,]$g2
g1p <- plot_refs[plot_refs$ref==ref,]$g1p
g2p <- plot_refs[plot_refs$ref==ref,]$g2p
nbin1 <- plot_refs[plot_refs$ref==ref,]$bin1
nbin2 <- plot_refs[plot_refs$ref==ref,]$bin2
c13b1 <- plot_refs[plot_refs$ref==ref,]$c13b1
c13b2 <- plot_refs[plot_refs$ref==ref,]$c13b2
rDNA1 <- plot_refs[plot_refs$ref==ref,]$rDNA1
rDNA2 <- plot_refs[plot_refs$ref==ref,]$rDNA2
rDNA <- c(rDNA1,rDNA2)

# load ref annotations
chrannot <- read.table(paste("nobackup/annotations/",ref,".",bsize,".chr_annotations",sep=""))
colnames(chrannot) <- c("chr","start","end","cen","chrbin")
binannot <- read.table(paste("nobackup/annotations/",ref,".",bsize,".bin_annotations",sep=""))

# load data
x <- 0:(nbin2-1)
empty <- data.frame(cbind(rep(x,each=nbin2),rep(x,nbin2)))
colnames(empty) <- c("bin1","bin2")

#add bin annotations
annotated <- merge(empty,binannot,by.x="bin1",by.y="V1")
colnames(annotated) <- c("bin1","bin2","chr1","cen1","arm1","chrn1")
annotated <- merge(annotated,binannot,by.x="bin2",by.y="V1")
colnames(annotated) <- c("bin2","bin1","chr1","cen1","arm1","chrn1","chr2","cen2","arm2","chrn2")
annotated <- annotated[order(annotated$bin1,annotated$bin2),]
annotated$bins <- paste(annotated$bin1,annotated$bin2,sep="-")
annotated$tel1 <- annotated$arm1-annotated$cen1
annotated$tel2 <- annotated$arm2-annotated$cen2

a <- read.table("nobackup/annotations/ScSu.32000.homology_noisolated")
b <- read.table("nobackup/annotations/ScSu.32000.homology_neighbors")

png(paste(outdir,"/",date,"_hom_prox_analysis_all.png",sep=""),2.4,2.4,"in",res=600)
ggplot(a) + geom_tile(aes(x=V1,y=V2)) + 
  scale_x_continuous(limits=c(0,nbin1),breaks=chrannot[1:16,]$cen,labels=chrlabs.odd ,expand=c(0,1)) + 
  scale_y_continuous(limits=c(nbin1,nbin2),breaks=chrannot[17:32,]$cen,labels=chrlabs.odd ,expand=c(0,1)) + paper.fig.rot + 
  geom_rect(aes(xmin=50,xmax=150,ymin=430,ymax=530),size=0.25,color="black",fill=NA)+
  geom_tile(data=subset(annotated,(bin1==100 | bin2==480)),aes(x=bin1,y=bin2,fill=(bin1==100 & bin2==480)),color=NA) + coord_fixed() + xlab(expression(italic("S. cerevisiae"))) + ylab(expression(italic("S. uvarum"))) + scale_fill_manual(values=c("grey","red")) + theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"cm"))
dev.off()

png(paste(outdir,"/",date,"_hom_prox_analysis_dcen.png",sep=""),2.4,2.4,"in",res=600)
ggplot(a) + geom_tile(aes(x=V1,y=V2)) + 
  scale_x_continuous(limits=c(0,nbin1),breaks=chrannot[1:16,]$cen,labels=chrlabs.odd ,expand=c(0,1)) + 
  scale_y_continuous(limits=c(nbin1,nbin2),breaks=chrannot[17:32,]$cen,labels=chrlabs.odd ,expand=c(0,1)) + paper.fig.rot + 
  geom_rect(aes(xmin=50,xmax=150,ymin=430,ymax=530),size=0.25,color="black",fill=NA)+
  geom_tile(data=subset(annotated,(bin1==100 | bin2==480) & cen1==4 & cen2==4),aes(x=bin1,y=bin2,fill=(bin1==100 & bin2==480),height=1,width=1),color=NA) + coord_fixed() + xlab(expression(italic("S. cerevisiae"))) + ylab(expression(italic("S. uvarum"))) + scale_fill_manual(values=c("grey","red")) + theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"cm"))
dev.off()

png(paste(outdir,"/",date,"_hom_prox_analysis_armlen.png",sep=""),2.4,2.4,"in",res=600)
ggplot(a) + geom_tile(aes(x=V1,y=V2)) + 
  scale_x_continuous(limits=c(0,nbin1),breaks=chrannot[1:16,]$cen,labels=chrlabs.odd,expand=c(0,1)) + 
  scale_y_continuous(limits=c(nbin1,nbin2),breaks=chrannot[17:32,]$cen,labels=chrlabs.odd,expand=c(0,1)) + paper.fig.rot + 
  geom_rect(aes(xmin=50,xmax=150,ymin=430,ymax=530),size=0.25,color="black",fill=NA)+
  geom_tile(data=subset(annotated,(bin1==100 | bin2==480) & cen1==4 & cen2==4 & (abs(arm1-14)/14 <= .25 & abs(arm2-13)/13 <= .25)),aes(x=bin1,y=bin2,fill=(bin1==100 & bin2==480),height=1,width=1),color=NA) + coord_fixed() + xlab(expression(italic("S. cerevisiae"))) + ylab(expression(italic("S. uvarum"))) + scale_fill_manual(values=c("grey","red")) + theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"cm"))
dev.off()

png(paste(outdir,"/",date,"_hom_prox_analysis_all_zoom.png",sep=""),2.4,2.4,"in",res=600)
ggplot(a) + geom_tile(aes(x=V1,y=V2)) + 
  scale_x_continuous(limits=c(50,150),breaks=chrannot[1:16,]$cen,labels=chrlabs.all,expand=c(0,1)) + 
  scale_y_continuous(limits=c(430,530),breaks=chrannot[17:32,]$cen,labels=chrlabs.all,expand=c(0,1)) + paper.fig.rot + 
  geom_tile(data=subset(annotated,(bin1==100 | bin2==480)),aes(x=bin1,y=bin2,fill=(bin1==100 & bin2==480)),color=NA) + coord_fixed() + xlab(expression(italic("S. cerevisiae"))) + ylab(expression(italic("S. uvarum"))) + scale_fill_manual(values=c("grey","red")) + theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"cm"))
dev.off()

png(paste(outdir,"/",date,"_hom_prox_analysis_dcen_zoom.png",sep=""),2.4,2.4,"in",res=600)
ggplot(a) + geom_tile(aes(x=V1,y=V2)) + 
  scale_x_continuous(limits=c(50,150),breaks=chrannot[1:16,]$cen,labels=chrlabs.all,expand=c(0,1)) + 
  scale_y_continuous(limits=c(430,530),breaks=chrannot[17:32,]$cen,labels=chrlabs.all,expand=c(0,1)) + paper.fig.rot + 
  geom_tile(data=subset(annotated,(bin1==100 | bin2==480) & cen1==4 & cen2==4),aes(x=bin1,y=bin2,fill=(bin1==100 & bin2==480),height=1,width=1),color=NA) + coord_fixed() + xlab(expression(italic("S. cerevisiae"))) + ylab(expression(italic("S. uvarum"))) + scale_fill_manual(values=c("grey","red")) + theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"cm"))
dev.off()

png(paste(outdir,"/",date,"_hom_prox_analysis_armlen_zoom.png",sep=""),2.4,2.4,"in",res=600)
ggplot(a) + geom_tile(aes(x=V1,y=V2)) + 
  scale_x_continuous(limits=c(50,150),breaks=chrannot[1:16,]$cen,labels=chrlabs.all,expand=c(0,1)) + 
  scale_y_continuous(limits=c(430,530),breaks=chrannot[17:32,]$cen,labels=chrlabs.all,expand=c(0,1)) + paper.fig.rot + 
  geom_tile(data=subset(annotated,(bin1==100 | bin2==480) & cen1==4 & cen2==4 & (abs(arm1-14)/14 <= .25 & abs(arm2-13)/13 <= .25)),aes(x=bin1,y=bin2,fill=(bin1==100 & bin2==480),height=1,width=1),color=NA) + coord_fixed() + xlab(expression(italic("S. cerevisiae"))) + ylab(expression(italic("S. uvarum"))) + scale_fill_manual(values=c("grey","red")) + theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"cm"))
dev.off()

# synteny dot plots ----

Sc <- read.table("nobackup/annotations/Scer.chrom_lengths")
Sc$cum <- cumsum(c(0,Sc$V2))[-17]
Sc$cum2 <- cumsum(Sc$V2)

Sc.cen <- read.table("references/Scer_cen.gff")
colnames(Sc.cen) <- c("chr","start","end")
Sc.cen <- merge(Sc.cen, Sc, by.x="chr", by.y="V1")
Sc.cen$value <- Sc.cen$cum + (Sc.cen$start + Sc.cen$end)/2
Sc.cen <- Sc.cen[order(Sc.cen$value),]

genomes <- data.frame(cbind(c("Spar_original","Spar","Suva_original","Suva"),
                            c("S. paradoxus","S. paradoxus","S. uvarum","S. uvarum")))

refnames <- c("Spar_original","Spar","Suva_original","Suva")
refroot <- c("Spar","Spar_revised","Suva","Suva_revised")

for (i in 1:4) {
  g2 <- refnames[i]
  g2r <- refroot[i]
  Sp <- read.table(paste("nobackup/annotations/",g2,".chrom_lengths",sep=""))
  Sp$cum <- cumsum(c(0,Sp$V2))[-17]
  Sp$cum2 <- cumsum(Sp$V2)
  
  ScSp <- read.table(paste("nobackup/annotations/Scer_",g2,"_gene_homology.txt",sep=""))
  colnames(ScSp) <- c("SGD","Scer_chr","Scer_start","Scer_end","Scer_dir","Spar_chr","Spar_start","Spar_end","Spar_dir")
  ScSp <- merge(ScSp, Sc, by.x="Scer_chr", by.y="V1")
  ScSp <- merge(ScSp, Sp, by.x="Spar_chr", by.y="V1")
  
  Sp.cen <- read.table(paste("references/",g2r,"_cen.gff",sep=""))
  colnames(Sp.cen) <- c("chr","start","end")
  Sp.cen <- merge(Sp.cen, Sp, by.x="chr", by.y="V1")
  Sp.cen$value <- Sp.cen$cum + (Sp.cen$start + Sp.cen$end)/2
  Sp.cen <- Sp.cen[order(Sp.cen$value),]
  
  ScSp$Scer_start <- ScSp$cum.x + ScSp$Scer_start
  ScSp$Scer_end <- ScSp$cum.x + ScSp$Scer_end
  ScSp$Spar_start <- ScSp$cum.y + ScSp$Spar_start
  ScSp$Spar_end <- ScSp$cum.y + ScSp$Spar_end
  
  png(paste(outdir,"/",date,"_Scer-",g2,"_homology_dotplot.png",sep=""),3,3,"in",res=1200)
  print(ggplot(ScSp) + geom_point(aes(x=(Scer_start+Scer_end)/2,y=(Spar_start+Spar_end)/2),size=0.1) + geom_vline(data=Sc,aes(xintercept=cum2),size=0.25) + geom_hline(data=Sp,aes(yintercept=cum2),size=0.25) + theme_classic() + theme(text=element_text(size=7),axis.text=element_text(size=7),axis.title=element_text(face="italic")) + scale_x_continuous(expand=c(0,0),breaks=Sc.cen$value,labels=chrlabs.odd,limits=c(0,max(Sc$cum2))) + scale_y_continuous(expand=c(0,0),breaks=Sp.cen$value,labels=chrlabs.odd,limits=c(0,max(Sp$cum2))) + coord_fixed() + xlab("S. cerevisiae") + ylab(genomes[genomes$X1==g2,]$X2))
  dev.off()
}

# chr13 interactions ----

chr13data1 <- read.table("nobackup/assigned/ScSu/ILY456_saturated_Sau3AI.chr13.assigned")
chr13data2 <- read.table("nobackup/assigned/ScSu/ILY456_exponential_Sau3AI.chr13.assigned")

chr13data1$cer <- 0
chr13data1$uva <- 0
chr13data1[chr13data1$V1=="Scer_13",]$cer <- chr13data1[chr13data1$V1=="Scer_13",]$V5
chr13data1[chr13data1$V1=="Scer_13",]$uva <- chr13data1[chr13data1$V1=="Scer_13",]$V11
chr13data1[chr13data1$V1=="Suva_13",]$uva <- chr13data1[chr13data1$V1=="Suva_13",]$V5
chr13data1[chr13data1$V1=="Suva_13",]$cer <- chr13data1[chr13data1$V1=="Suva_13",]$V11
chr13data1$cond <- "sat"
chr13data2$cer <- 0
chr13data2$uva <- 0
chr13data2[chr13data2$V1=="Scer_13",]$cer <- chr13data2[chr13data2$V1=="Scer_13",]$V5
chr13data2[chr13data2$V1=="Scer_13",]$uva <- chr13data2[chr13data2$V1=="Scer_13",]$V11
chr13data2[chr13data2$V1=="Suva_13",]$uva <- chr13data2[chr13data2$V1=="Suva_13",]$V5
chr13data2[chr13data2$V1=="Suva_13",]$cer <- chr13data2[chr13data2$V1=="Suva_13",]$V11
chr13data2$cond <- "exp"
chr13data <- rbind(chr13data1,chr13data2)
chr13data[chr13data$cond=="exp",]$cond <- "Exponential growth"
chr13data[chr13data$cond=="sat",]$cond <- "Saturated"

pdf(paste(outdir,"/",date,"_chr13_highres_heatmaps.pdf",sep=""),3.2,1.8)
ggplot(chr13data) + stat_bin_2d(aes(x=cer/1000,y=uva/1000),binwidth = 20) + scale_fill_gradientn(colors = heatmapcols) + coord_fixed() + facet_wrap(~cond) + theme_classic() + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + xlab("") + ylab("") + paper.fig.rot + theme(strip.background=element_blank(),plot.margin=unit(c(0,0,0,0),"cm"))
dev.off()

# browser shot ----

region <- read.table("plotting/840k-860k.txt",header=FALSE,stringsAsFactors=FALSE)
region$end <- as.numeric(region$V7=="+")*region$V5 + as.numeric(region$V7=="-")*region$V4
region$end2 <- as.numeric(region$V7=="+")*(region$V5+400) + as.numeric(region$V7=="-")*(region$V4-400)

dels <- read.table("plotting/has1dels.txt",header=FALSE,stringsAsFactors=FALSE)

pdf(paste(outdir,"/",date,"_chr13_browser.pdf",sep=""),3.6,1.5)
ggplot(region) + geom_rect(aes(xmin=V4,xmax=V5,ymin=V10,ymax=V10+.3),fill="royalblue",color=NA) + 
  geom_text(aes(x=V12,y=V10+.65,label=V9),fontface="italic",size=2.469444) + 
  xlab("") + paper.fig.rot + 
  theme(axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) + 
  ylab("") + geom_segment(aes(x=end,xend=end2,y=V10+.15,yend=V10+.15),color="black",arrow = arrow(length = unit(0.02, "npc"))) + ylim(0,4.7) + scale_x_continuous(labels=comma) + geom_segment(data=dels,aes(x=V2,xend=V3,y=V4,yend=V4)) + geom_segment(data=dels,aes(x=V2,xend=V2,y=V4-.1,yend=V4+.1)) + geom_segment(data=dels,aes(x=V3,xend=V3,y=V4-.1,yend=V4+.1)) + theme(plot.margin=unit(c(0,0,0,0),"cm"))
dev.off()


# heatmap color bar ----

pdf(paste(outdir,"/",date,"_colorbar.pdf",sep=""),2,2)
ggplot(region) + geom_tile(aes(x=1,y=1,fill=1)) + scale_fill_gradientn(colors=heatmapcols,limits=heatmaprange,name="") + paper.fig.rot + theme(legend.key.height=unit(0.3,"in"),legend.key.width=unit(0.2,"in"))
dev.off()

# raw heatmaps for genome revisions ----

matrices <- read.table("plotting/matrices_revisions.txt",stringsAsFactors = FALSE)

# S. paradoxus chromosome 4
for (i in 1:2) {
  data <- read.table(paste("nobackup/matrix/",matrices[i,]$V2,"/",matrices[i,]$V3,"/",matrices[i,]$V1,".",matrices[i,]$V4,".matrix",sep=""))
  chrannot <- read.table(paste("nobackup/annotations/",matrices[i,]$V2,".",matrices[i,]$V3,".chr_annotations",sep=""))
  png(paste(outdir,"/",date,"_",matrices[i,]$V2,"_chr4_heatmap.png",sep=""),3,3,"in",res=1200)
  print(ggplot(data) + 
          geom_raster(aes(x=V1,y=V2,fill=log10(V3+1))) + 
          paper.fig.rot + hm + coord_fixed() + 
          scale_fill_gradientn(colors=heatmapcols, name = "") +
          scale_x_continuous(limits=c(chrannot[4,]$V2+.5,
                                      chrannot[4,]$V3-.5),
                             breaks=seq(chrannot[4,]$V2,
                                        chrannot[4,]$V3,20)+.5,
                             labels=seq(0,chrannot[4,]$V3-chrannot[4,]$V2,20)*16,
                             expand=c(0,0)) + 
          scale_y_continuous(limits=c(chrannot[4,]$V2+.5,
                                      chrannot[4,]$V3-.5),
                             breaks=seq(chrannot[4,]$V2,
                                        chrannot[4,]$V3,20)+.5,
                             labels=seq(0,chrannot[4,]$V3-chrannot[4,]$V2,20)*16,
                             expand=c(0,0)) +
          theme(legend.text=element_text(size=7),
                legend.key.height=unit(0.75,"cm"),
                legend.key.width=unit(0.25,"cm"))
  )
  dev.off()
}

# S. uvarum chromosome 3
for (i in 3:4) {
  data <- read.table(paste("nobackup/matrix/",matrices[i,]$V2,"/",matrices[i,]$V3,"/",matrices[i,]$V1,".",matrices[i,]$V4,".matrix",sep=""))
  chrannot <- read.table(paste("nobackup/annotations/",matrices[i,]$V2,".",matrices[i,]$V3,".chr_annotations",sep=""))
  png(paste(outdir,"/",date,"_",matrices[i,]$V2,"_chr3_heatmap.png",sep=""),3,3,"in",res=1200)
  print(ggplot(data) + 
          geom_raster(aes(x=V1,y=V2,fill=log10(V3+1))) + 
          paper.fig.rot + hm + coord_fixed() + 
          scale_fill_gradientn(colors=heatmapcols, name = "") +
          scale_x_continuous(limits=c(chrannot[3,]$V2+.5,
                                      chrannot[3,]$V3-.5),
                             breaks=seq(chrannot[3,]$V2,
                                        chrannot[3,]$V3,5)+.5,
                             labels=seq(0,chrannot[3,]$V3-chrannot[3,]$V2,5)*16,
                             expand=c(0,0)) + 
          scale_y_continuous(limits=c(chrannot[3,]$V2+.5,
                                      chrannot[3,]$V3-.5),
                             breaks=seq(chrannot[3,]$V2,
                                        chrannot[3,]$V3,5)+.5,
                             labels=seq(0,chrannot[3,]$V3-chrannot[3,]$V2,5)*16,
                             expand=c(0,0)) +
          theme(legend.text=element_text(size=7),
                legend.key.height=unit(0.75,"cm"),
                legend.key.width=unit(0.25,"cm"))
  )
  dev.off()
}

# generalized analysis/plotting pipeline ----

# load reference data
plot_refs <- read.table("plotting/plotting_refs.txt",stringsAsFactors=FALSE,sep = "\t")
colnames(plot_refs) <- c("ref","g1","g2","g1p","g2p","rDNA1","rDNA2","bin1","bin2","c13b1","c13b2")

# experimental samples
plotting_expt <- read.table("matrices.txt",stringsAsFactors = FALSE)
colnames(plotting_expt) <- c("sample","ref","bsize","mask","filt")
plotting_expt$matr <- paste("nobackup/matrix/",plotting_expt$ref,"/",plotting_expt$bsize,"/",plotting_expt$sample,".",plotting_expt$mask,".matrix.txt",sep="")
plotting_expt$rsum <- paste("nobackup/matrix/",plotting_expt$ref,"/",plotting_expt$bsize,"/",plotting_expt$sample,".",plotting_expt$mask,".rowsums.txt",sep="")

# polymer model samples
plotting_model <- read.table("model/matrices_model.txt",stringsAsFactors = FALSE)
colnames(plotting_model) <- c("sample","ref","bsize","mask","filt","root")
plotting_model$matr <- paste(plotting_model$root,".matrix.txt",sep="")
plotting_model$rsum <- paste(plotting_model$root,".rowsums.txt",sep="")
plotting_model$root <- NULL

# combine all samples
plotting <- rbind(plotting, plotting_model)

# options whether to save each type of plot or run each type of analysis
saverowsums <- TRUE
savenormheatmap <- TRUE
savenormheatmaphom <- TRUE
savenorm13heatmap <- TRUE
savenorm13comp <- TRUE
savenorm13hist <- TRUE
saverDNAheatmap <- TRUE
saverDNA5heatmap <- TRUE
savechr14knockin <- TRUE

# options of which samples to process
jmin <- 1
jmax <- nrow(plotting)
js <- jmin:jmax

# homolog proximity analysis options
hom_analysis <- TRUE
exclusions <- c("none","rDNA-cen")
cen_thrs <- c(0,35)
arm_thrs <- c(0.25,20)
min_n <- 2

# loop through samples
for (j in js) {
  sample <- plotting[j,]$sample
  print(sample)
  mask <- plotting[j,]$mask
  ref <- plotting[j,]$ref
  bsize <- plotting[j,]$bsize
  filt <- plotting[j,]$filt
  matr <- plotting[j,]$matr
  rsum <- plotting[j,]$rsum
  
  # get reference information for appropriate reference
  g1 <- plot_refs[plot_refs$ref==ref,]$g1
  g2 <- plot_refs[plot_refs$ref==ref,]$g2
  g1p <- plot_refs[plot_refs$ref==ref,]$g1p
  g2p <- plot_refs[plot_refs$ref==ref,]$g2p
  nbin1 <- plot_refs[plot_refs$ref==ref,]$bin1
  nbin2 <- plot_refs[plot_refs$ref==ref,]$bin2
  c13b1 <- plot_refs[plot_refs$ref==ref,]$c13b1
  c13b2 <- plot_refs[plot_refs$ref==ref,]$c13b2
  rDNA1 <- plot_refs[plot_refs$ref==ref,]$rDNA1
  rDNA2 <- plot_refs[plot_refs$ref==ref,]$rDNA2
  rDNA <- c(rDNA1,rDNA2)
  
  # homologous order for S. uvarum
  hom_order <- data.frame(cbind(seq(0,nbin2-1),seq(0,nbin2-1)))
  if (g2 == "Suva") {
    hom_order <- data.frame(cbind(c(seq(0,nbin1-1),
                                    seq(nbin1,nbin1+5),                                 # Scer_1 = Suva_1
                                    seq(nbin1+6,nbin1+14), seq(nbin1+71,nbin1+87),      # Scer_2 = Suva_2, Suva_4
                                    seq(nbin1+47,nbin1+56),                             # Scer_3 = Suva_3
                                    seq(nbin1+57,nbin1+70), seq(nbin1+15,nbin1+46),     # Scer_4 = Suva_4, Suva_2
                                    seq(nbin1+88,nbin1+104),                            # Scer_5 = Suva_5
                                    seq(nbin1+105,nbin1+110), seq(nbin1+195,nbin1+195), # Scer_6 = Suva_6, Suva_10
                                    seq(nbin1+122,nbin1+154),                           # Scer_7 = Suva_7
                                    seq(nbin1+155,nbin1+158), seq(nbin1+324,nbin1+337), # Scer_8 = Suva_8, Suva_15
                                    seq(nbin1+181,nbin1+194),                           # Scer_9 = Suva_9
                                    seq(nbin1+121,nbin1+111), seq(nbin1+196,nbin1+208), # Scer_10 = Suva_6, Suva_10
                                    seq(nbin1+209,nbin1+313),                           # chr11, chr12, chr13, chr14
                                    seq(nbin1+314,nbin1+323), seq(nbin1+159,nbin1+180), # Scer_15 = Suva_15, Suva_8
                                    seq(nbin1+338,nbin1+366)),                          # Scer_16 = Suva_16
                                  seq(0,nbin2-1))) 
  }
  hom_order <- hom_order[order(hom_order$X1),]
  
  # set path to save output
  outdir <- "nobackup/figures"
  dir.create(paste(outdir,"/",ref,sep=""))
  dir.create(paste(outdir,"/",ref,"/",bsize,sep=""))
  outdir <- paste(outdir,ref,bsize,sep="/")
  
  # load ref annotations
  chrannot <- read.table(paste("nobackup/annotations/",ref,".",bsize,".chr_annotations",sep=""))
  colnames(chrannot) <- c("chr","start","end","cen","chrbin")
  binannot <- read.table(paste("nobackup/annotations/",ref,".",bsize,".bin_annotations",sep=""))
  
  # make directories
  dir.create(paste(outdir,"/rowsums",sep=""))
  dir.create(paste(outdir,"/normed_heatmap",sep=""))
  dir.create(paste(outdir,"/normed_heatmap_homorder",sep=""))
  dir.create(paste(outdir,"/normed_heatmap_chr13",sep=""))
  dir.create(paste(outdir,"/normed_heatmap_13-14",sep=""))
  dir.create(paste(outdir,"/normed_heatmap_rDNA",sep=""))
  dir.create(paste(outdir,"/normed_heatmap_rDNA5",sep=""))
  dir.create(paste(outdir,"/chr13_interaction_hist",sep=""))
  dir.create(paste(outdir,"/chr13_knockin_hist",sep=""))
  for (exclusion in exclusions) {
    dir.create(paste(outdir,"/",exclusion,sep=""))
    for (cen_thr in cen_thrs) {
      dir.create(paste(outdir,"/",exclusion,"/cen",cen_thr,sep=""))
      for (arm_thr in arm_thrs) {
        dir.create(paste(outdir,"/",exclusion,"/cen",cen_thr,"/arm",arm_thr,sep=""))
        dir.create(paste(outdir,"/",exclusion,"/cen",cen_thr,"/arm",arm_thr,"/permtest",sep=""))
        dir.create(paste(outdir,"/",exclusion,"/cen",cen_thr,"/arm",arm_thr,"/genome",sep=""))
      }
    }
  }
  
  # load data
  data <- read.table(matr)
  colnames(data) <- c("bin1","bin2","raw","norm")
  # label positions
  yfrac = 0.08
  xfrac = 0.1
  
  #rowsums
  raw.rowsum <- read.table(rsum)
  colnames(raw.rowsum) <- c("bin","sum")

  if (saverowsums) {
    pdf(paste(outdir,"/rowsums/",date,"_",sample,"_",mask,".pdf",sep=""),3.5,1.7)
    p <- ggplot(raw.rowsum) + geom_line(aes(x=bin,y=sum/1000)) + paper.fig.rot + 
      scale_x_continuous(breaks=c(chrannot$cen),labels=c(chrlabs.odd,chrlabs.odd),expand=c(0,0)) + 
      scale_y_continuous(expand=c(0,0)) + xlab("") + ylab("Reads per bin (thousands)") +
      geom_point(data=subset(raw.rowsum,sum<max(chrannot$end)*filt),aes(x=bin,y=sum/1000),color="red") 
    print(p)
    dev.off()
  }
  
  # normed heatmap linear scale
  if (savenormheatmap) {
    png(paste(outdir,"/normed_heatmap/",date,"_",sample,"_",mask,".png",sep=""),3.5,3.5,"in",res=600)
    p <- ggplot(data) + geom_raster(aes(x=bin1,y=bin2,fill=norm)) + 
      paper.fig.rot + 
      hm + 
      coord_fixed() + 
      scale_fill_gradientn(colors=heatmapcols, limits=heatmaprange, oob=squish, na.value="grey50", name = "O/E") + 
      scale_x_continuous(breaks=chrannot$cen,labels=c(chrlabs.odd,chrlabs.odd),expand=c(0,0)) + 
      scale_y_continuous(breaks=chrannot$cen,labels=c(chrlabs.odd,chrlabs.odd),expand=c(0,0)) +
      geom_rect(data=chrannot,aes(xmax=end-.5,xmin=start-.5,ymax=end-.5,ymin=start-.5),size=0.25,color="black",fill=NA) + 
      theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),legend.position="none",axis.text = element_blank()) +
      geom_vline(aes(xintercept=nbin1),color="black",size=0.25) +
      geom_hline(aes(yintercept=nbin1),color="black",size=0.25)
    print(p)
    dev.off()
    
    # zoom into ScXI and ScXII
    png(paste(outdir,"/normed_heatmap/",date,"_",sample,"_",mask,"_ScXI-XII.png",sep=""),2.2,2.2,"in",res=600)
    p <- ggplot(subset(data,bin1 >= 211 & bin1 < 266 & bin2 >= 211 & bin2 < 266)) + geom_raster(aes(x=bin1,y=bin2,fill=norm)) + 
      paper.fig.rot + 
      hm + 
      coord_fixed() + 
      scale_fill_gradientn(colors=heatmapcols, limits=heatmaprange, oob=squish, na.value="grey50", name = "O/E") + 
      scale_x_continuous(limits=c(210.5,265.5),breaks=chrannot$cen,labels=c(chrlabs.all,chrlabs.all),expand=c(0,0)) + 
      scale_y_continuous(limits=c(210.5,265.5),breaks=chrannot$cen,labels=c(chrlabs.all,chrlabs.all),expand=c(0,0)) +
      theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),legend.position="none") +
      geom_vline(aes(xintercept=nbin1),color="black",size=0.25) +
      geom_hline(aes(yintercept=nbin1),color="black",size=0.25)
    print(p)
    dev.off()
  }
  
  #normed heatmap linear scale, homology order
  if (savenormheatmaphom) {
    png(paste(outdir,"/normed_heatmap_homorder/",date,"_",sample,"_",mask,".png",sep=""),3.5,3.5,"in",res=600)
    p <- ggplot(data) + geom_raster(aes(x=hom_order[bin1+1,]$X2,y=hom_order[bin2+1,]$X2,fill=norm)) + 
      paper.fig.rot + 
      hm + 
      coord_fixed() + 
      scale_fill_gradientn(colors=heatmapcols, limits=heatmaprange, oob=squish, na.value="grey50", name = "O/E") + 
      scale_x_continuous(breaks=hom_order[chrannot$cen+1,]$X2,labels=chrannot$V6,expand=c(0,0)) + 
      scale_y_continuous(breaks=hom_order[chrannot$cen+1,]$X2,labels=chrannot$V6,expand=c(0,0)) +
      theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),legend.position="none",axis.text = element_blank()) +
      geom_vline(aes(xintercept=nbin1),color="black",size=0.25) +
      geom_hline(aes(yintercept=nbin1),color="black",size=0.25)
    print(p)  
    dev.off()
  }
  
  #add bin annotations
  annotated <- merge(data,binannot,by.x="bin1",by.y="V1")
  colnames(annotated) <- c("bin1","bin2","raw","norm","chr1","cen1","arm1","chrn1")
  annotated <- merge(annotated,binannot,by.x="bin2",by.y="V1")
  colnames(annotated) <- c("bin2","bin1","raw","norm","chr1","cen1","arm1","chrn1","chr2","cen2","arm2","chrn2")
  annotated <- annotated[order(annotated$bin1,annotated$bin2),]
  annotated$bins <- paste(annotated$bin1,annotated$bin2,sep="-")
  annotated$tel1 <- annotated$arm1-annotated$cen1
  annotated$tel2 <- annotated$arm2-annotated$cen2
  
  #normed chr13 heatmap
  subbins1 <- seq(subset(chrannot,chr==paste(g1,"_13",sep=""))$start,subset(chrannot,chr==paste(g1,"_14",sep=""))$start-1)
  subbins2 <- seq(subset(chrannot,chr==paste(g2,"_13",sep=""))$start,subset(chrannot,chr==paste(g2,"_14",sep=""))$start-1)
  annotated.sub <- subset(annotated,bin1 %in% subbins1 & bin2 %in% subbins2)
  
  if (savenorm13heatmap) {
    pdf(paste(outdir,"/normed_heatmap_chr13/",date,"_",sample,"_",mask,".pdf",sep=""),1.2,1.2)
    p <- ggplot(annotated.sub) + geom_tile(aes(x=bin1,y=bin2,fill=norm)) + 
      paper.fig.rot + 
      paperhm + 
      coord_fixed() + 
      scale_fill_gradientn(colors=heatmapcols, limits=heatmaprange, oob=squish, na.value="grey50", name = "O/E") + 
      theme(legend.position="none",text=element_text(size=6)) +
      scale_x_continuous(breaks=seq(subbins1[8],subbins1[length(subbins1)],8)+.5,labels=seq(bsize*8/1000,bsize/1000*(chrannot[chrannot$chr==paste(g1,"_13",sep=""),]$end-chrannot[chrannot$chr==paste(g1,"_13",sep=""),]$start),bsize*8/1000),expand=c(0,0)) +
      scale_y_continuous(breaks=seq(subbins2[8],subbins2[length(subbins2)],8)+.5,labels=seq(bsize*8/1000,bsize/1000*(chrannot[chrannot$chr==paste(g2,"_13",sep=""),]$end-chrannot[chrannot$chr==paste(g2,"_13",sep=""),]$start),bsize*8/1000),expand=c(0,0))
    print(p)
    dev.off()
  }
  
  # write HAS1 comparisons to file for aggregated figure
  if (savenorm13comp) {
    chr13comp <- subset(annotated,
                        cen1 >= 15 & cen2 >= 15 & 
                          tel1 > 1 & tel2 > 1 & 
                          bin1 < chrannot[chrannot$chr==paste(g2,"_1",sep=""),2] & 
                          bin2 >= chrannot[chrannot$chr==paste(g2,"_1",sep=""),2] & 
                          (chr1 != rDNA1 | chr2 != rDNA2))
    chr13pair <- data.frame(subset(annotated,bin1==c13b1 & bin2==c13b2)$norm)
    chr13comp.df <- as.data.frame(chr13comp$norm)
    
    write.table(chr13comp.df,paste(outdir,"/chr13_interaction_hist/",date,"_",sample,"_",mask,"_comp.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep=" ")
    write.table(chr13pair,paste(outdir,"/chr13_interaction_hist/",date,"_",sample,"_",mask,"_pair.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep=" ")
  }
  
  # HAS1 comparison histogram
  if (savenorm13hist) {
    pdf(paste(outdir,"/chr13_interaction_hist/",date,"_",sample,"_",mask,".pdf",sep=""),5,4)
    print(ggplot(subset(annotated,
                        cen1 >= 15 & cen2 >= 15 & 
                          tel1 > 1 & tel2 > 1 & 
                          bin1 < chrannot[chrannot$chr==paste(g2,"_1",sep=""),2] & 
                          bin2 >= chrannot[chrannot$chr==paste(g2,"_1",sep=""),2] & 
                          (chr1 != rDNA1 | chr2 != rDNA2))) + 
            geom_histogram(aes(x=norm),width=1,binwidth=.025,fill="grey") + 
            geom_vline(aes(xintercept=subset(annotated,bin1==c13b1 & bin2==c13b2)$norm),color="red",size=2) + 
            geom_vline(aes(xintercept=mean(subset(annotated,cen1 >= 15 & cen2 >= 15 & 
                                                    tel1 > 1 & tel2 > 1 & 
                                                    bin1 < chrannot[chrannot$chr==paste(g2,"_1",sep=""),2] & 
                                                    bin2 >= chrannot[chrannot$chr==paste(g2,"_1",sep=""),2] & 
                                                    (chr1 != rDNA1 | chr2 != rDNA2))$norm,na.rm=TRUE)),color="royalblue",size=2) + 
            bl2 + xlab("normalized interaction frequency") + 
            theme(axis.text=element_text(size=20),text=element_text(size=20)) + 
            scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)))
    dev.off()
  }
  
  marg <- 0 # margin (in bins) around rDNA-carrying chromosomes
  #marg <- 15
  subbins1 <- seq(subset(chrannot,chr==rDNA1)$start-marg,subset(chrannot,chr==rDNA1)$end-1+marg)
  subbins2 <- seq(subset(chrannot,chr==rDNA2)$start-marg,subset(chrannot,chr==rDNA2)$end-1+marg)
  annotated.sub <- subset(annotated,bin1 %in% subbins1 & bin2 %in% subbins2)
  
  # heatmap zoomed to rDNA-carrying chromosomes
  if(saverDNAheatmap) {
    pdf(paste(outdir,"/normed_heatmap_rDNA/",date,"_",sample,"_",mask,".pdf",sep=""),1,1)
    p <- ggplot(annotated.sub) + geom_tile(aes(x=bin1,y=bin2,fill=norm)) + 
      paper.fig.rot + 
      paperhm + 
      coord_fixed() + 
      scale_fill_gradientn(colors=heatmapcols, limits=heatmaprange, oob=squish, na.value="grey50", name = "O/E") + 
      theme(plot.margin=unit(c(0,0,0,0),"cm"),legend.position="none",axis.ticks = element_blank(),axis.text=element_blank()) +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0))
    print(p)    
    dev.off()
  }

  # heatmap of interactions between ScV and SuXII  
  if(saverDNA5heatmap) {
    pdf(paste(outdir,"/normed_heatmap_rDNA5/",date,"_",sample,"_",mask,".pdf",sep=""),0.5,1)
    p <- ggplot(annotated) + geom_tile(aes(x=bin1,y=bin2,fill=norm)) + 
      paper.fig.rot + 
      paperhm + 
      coord_fixed() + 
      scale_fill_gradientn(colors=heatmapcols, limits=heatmaprange, oob=squish, na.value="grey50", name = "O/E") + 
      theme(plot.margin=unit(c(0,0,0,0),"cm"),legend.position="none",axis.ticks = element_blank(),axis.text=element_blank()) +
      scale_x_continuous(expand=c(0,0),limits=c(subset(chrannot,chr==paste(g1,"_5",sep=""))$start-.5,subset(chrannot,chr==paste(g1,"_5",sep=""))$end-.5)) +
      scale_y_continuous(expand=c(0,0),limits=c(subset(chrannot,chr==rDNA2)$start-.5,subset(chrannot,chr==rDNA2)$end-.5))
    print(p)
    dev.off()
  }
  
  #normed chr13-chr14 heatmap
  if (savechr14knockin) {
    subbins1 <- seq(subset(chrannot,chr==paste(g1,"_13",sep=""))$start,subset(chrannot,chr==paste(g1,"_15",sep=""))$start-1)
    subbins2 <- seq(subset(chrannot,chr==paste(g2,"_13",sep=""))$start,subset(chrannot,chr==paste(g2,"_14",sep=""))$start-1)
    annotated.sub <- subset(annotated,bin1 %in% subbins1 & bin2 %in% subbins2)
  
    pdf(paste(outdir,"/normed_heatmap_13-14/",date,"_",sample,"_",mask,".pdf",sep=""),2.4,1.6)
    p <- ggplot(annotated.sub) + geom_tile(aes(x=bin1,y=bin2,fill=norm)) + 
      paper.fig.rot + 
      paperhm + 
      coord_fixed() + 
      geom_vline(aes(xintercept=subset(chrannot,chr==paste(g1,"_14",sep=""))$start-.5)) +
      scale_fill_gradientn(colors=heatmapcols, limits=heatmaprange, oob=squish, na.value="grey50", name = "O/E") + 
      theme(legend.position="none",text=element_text(size=6)) +
      scale_x_continuous(breaks=c(273,281,289,302,310,318)+.5,labels=c(256,512,768,256,512,768),expand=c(0,0)) +
      scale_y_continuous(breaks=c(653,661,669)+.5,labels=c(256,512,768),expand=c(0,0))
    print(p)
    dev.off()
  
    pdf(paste(outdir,"/chr13_knockin_hist/",date,"_",sample,"_",mask,".pdf",sep=""),1.6,1.6)
    print(ggplot(subset(annotated,
                        cen1 >= 15 & cen2 >= 15 & 
                          tel1 > 1 & tel2 > 1 & 
                          bin1 < chrannot[chrannot$chr==paste(g2,"_1",sep=""),2] & 
                          bin2 >= chrannot[chrannot$chr==paste(g2,"_1",sep=""),2] & 
                          (chr1 != rDNA1 | chr2 != rDNA2))) + 
            geom_histogram(aes(x=norm),width=1,binwidth=.025,fill="grey") + 
            geom_vline(aes(xintercept=subset(annotated,bin1==299 & bin2==c13b2)$norm),color="red",size=0.5) + 
            geom_vline(aes(xintercept=subset(annotated,bin1==c13b1 & bin2==c13b2)$norm),color="royalblue",size=0.5) + 
            paper.fig.rot + xlab("Normalized interaction frequency") + ylab("Count") +
            scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)))
    dev.off()
  }
  
  # global homolog proximity analysis
  if (hom_analysis) {
    # load list of homologous bins
    hom <- read.table(paste("nobackup/annotations/",ref,".",bsize,".homology_noisolated",sep=""))
    hom$bins <- paste(hom$V1,hom$V2,sep="-")
    hom_neighbors <- read.table(paste("nobackup/annotations/",ref,".",bsize,".homology_neighbors",sep=""))
    hom_neighbors$bins <- paste(hom_neighbors$V1,hom_neighbors$V2,sep="-")
    
    # add homology annotations to raw matrix
    annotated$hom <- annotated$bins %in% hom$bins
    annotated$hom_neighb <- annotated$bins %in% hom_neighbors$bins
    
    for (exclusion in exclusions) {
      if (exclusion=="none") {
        nonhomologs <- subset(annotated, hom==FALSE & 
                                bin1 < nbin1 & bin2 >= nbin1 & hom_neighb==FALSE)
        homologs <- subset(annotated, hom==TRUE)
      }
      else if (exclusion=="rDNA-cen") {
        nonhomologs <- subset(annotated, hom==FALSE & 
                                chr1 != rDNA[1] & chr2 != rDNA[2] & 
                                cen1 != 0 & cen2 != 0 & 
                                bin1 < nbin1 & bin2 >= nbin1 & 
                                hom_neighb==FALSE)
        homologs <- subset(annotated, hom==TRUE & 
                             chr1 != rDNA[1] & chr2 != rDNA[2] & 
                             cen1 != 0 & cen2 != 0)
      }
      
      for (cen_thr in cen_thrs) {
        for (arm_thr in arm_thrs) {
          
          nonhomologous_set <- data.frame()
          homologous_set <- data.frame()
          pairing_data <- data.frame()
          nperm <- 10000
          perm_sums <- rep(0,nperm)
          
          for (i in seq(1,nrow(homologs))) {
            # same bin vs all similar
            new_entries <- subset(nonhomologs, 
                                  (bin1==homologs[i,]$bin1 & 
                                     abs(cen2-homologs[i,]$cen2) <= cen_thr & 
                                     abs(arm2-homologs[i,]$arm2)/homologs[i,]$cen2 <= arm_thr & 
                                     !is.na(norm)) |
                                    (bin2==homologs[i,]$bin2 & 
                                       abs(cen1-homologs[i,]$cen1) <= cen_thr & 
                                       abs(arm1-homologs[i,]$arm1)/homologs[i,]$arm1 <= arm_thr & 
                                       !is.na(norm))
            )

            if (!is.na(homologs[i,]$norm) && nrow(new_entries) >= min_n) {
              new_entries$hombin <- homologs[i,]$bins
              new_entries$harm1 <- homologs[i,]$arm1
              new_entries$harm2 <- homologs[i,]$arm2
              new_entries$hcen1 <- homologs[i,]$cen1
              new_entries$hcen2 <- homologs[i,]$cen2
              
              nonhomologous_set <- rbind(nonhomologous_set,new_entries)
              homologous_set <- rbind(homologous_set,homologs[i,])
              new_data <- data.frame(arm1 = homologs[i,]$arm1, arm2 = homologs[i,]$arm2, mean_diff = homologs[i,]$norm - mean(new_entries$norm),median_diff = homologs[i,]$norm - median(new_entries$norm),percentile = sum(homologs[i,]$norm > new_entries$norm)/nrow(new_entries),num_gt = sum(homologs[i,]$norm >= new_entries$norm),num_lt = sum(homologs[i,]$norm < new_entries$norm))
              pairing_data <- rbind(pairing_data,new_data)
              
              perm_sums <- perm_sums + sample(new_entries$norm,nperm,replace=TRUE)
            }
          }
          homologous_set$hombin <- homologous_set$bins
          homologous_set$harm1 <- homologous_set$arm1
          homologous_set$harm2 <- homologous_set$arm2
          homologous_set$hcen1 <- homologous_set$cen1
          homologous_set$hcen2 <- homologous_set$cen2

          perm_sums.df <- as.data.frame(perm_sums)
          colnames(perm_sums.df) <- "perm_sums"
          
          perm_sums.df$perm_sums <- perm_sums.df$perm_sums/nrow(homologous_set)
          write.table(perm_sums.df,paste(outdir,"/",exclusion,"/cen",cen_thr,"/arm",arm_thr,"/permtest/",date,"_",sample,"_",mask,"_min",min_n,".txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep=" ")
          
          hom_sum <- data.frame(mean(homologous_set$norm))
          write.table(hom_sum,paste(outdir,"/",exclusion,"/cen",cen_thr,"/arm",arm_thr,"/permtest/",date,"_",sample,"_",mask,"_min",min_n,"_hom.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep=" ")
          
          plot_homology <- homologous_set
          plot_homology$median_diff <- pairing_data$median_diff
          
          write.table(plot_homology,paste(outdir,"/",exclusion,"/cen",cen_thr,"/arm",arm_thr,"/genome/",date,"_",sample,"_",mask,"_min",min_n,".txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep=" ")
        }
      }
    }
  }
}

# combined permtest violin plot #1 ----
all_permtest <- data.frame()
outdir <- "nobackup/figures"
subdir <- "/ScSu/32000/none/cen35/arm20/permtest/"
filelist <- c(paste(outdir,subdir,"2016-12-05","_ILY456_saturated_Sau3AI_Surep-rDNA_min2",sep=""),
              paste(outdir,subdir,"2016-12-05","_ILY456_exponential_Sau3AI_Surep-rDNA_min2",sep=""),
              paste(outdir,subdir,"2016-12-05","_ScSu_v1_none_min2",sep="")
)
samplist <- c("saturated none",
              "exponential none",
              "model v1 none"
)
cols=c("Saturated","Exponential growth","Polymer model")
for (i in 1:length(filelist)) {
  toadd <- read.table(paste(filelist[i],".txt",sep=""))
  toadd_h <- read.table(paste(filelist[i],"_hom.txt",sep=""))
  toadd$diff <- (toadd_h[1,]-toadd$V1)/toadd$V1
  toadd$V2 <- samplist[i]
  toadd$cols <- cols[i]
  all_permtest <- rbind(all_permtest,toadd)  
}
customcols <- c(brewercols[1],brewercols[2],"grey50")

pdf(paste(outdir,"/",date,"_hom_prox_violin_1_L.pdf",sep=""),1.4,1.5)
ggplot(all_permtest) + 
  geom_violin(aes(x=factor(V2,levels=samplist),y=diff+1,fill=factor(cols,levels=cols)),alpha=1) + 
  paper.fig.rot + 
  scale_x_discrete(breaks=factor(c("exponential none","exponential cendist","exponential armlen"),levels=samplist)) +
  scale_y_continuous(limits=c(0.5,4.7),expand=c(0,0))+
  xlab("") + ylab("") + theme(axis.text.x=element_blank(),text=element_text(size=7),axis.text.y=element_text(size=7),legend.position="none") + 
  geom_hline(aes(yintercept=1),size=0.25) +
  scale_fill_manual(values = customcols,name="")
dev.off()

all_permtest <- data.frame()
subdir <- "/ScSu/32000/none/cen0/"
filelist <- c(paste(outdir,subdir,"arm20/permtest/2016-12-05_ILY456_saturated_Sau3AI_Surep-rDNA_min2",sep=""),
              paste(outdir,subdir,"arm20/permtest/2016-12-05_ILY456_exponential_Sau3AI_Surep-rDNA_min2",sep=""),
              paste(outdir,subdir,"arm20/permtest/2016-12-05_ScSu_v1_none_min2",sep=""),
              paste(outdir,subdir,"arm0.25/permtest/2016-12-05_ILY456_saturated_Sau3AI_Surep-rDNA_min2",sep=""),
              paste(outdir,subdir,"arm0.25/permtest/2016-12-05_ILY456_exponential_Sau3AI_Surep-rDNA_min2",sep=""),
              paste(outdir,subdir,"arm0.25/permtest/2016-12-05_ScSu_v1_none_min2",sep="")
)
samplist <- c("saturated cendist",
              "exponential cendist",
              "model v1 cendist",
              "saturated armlen",
              "exponential armlen",
              "model v1 armlen"
)
cols=c("Saturated",
       "Exponential growth",
       "Polymer model",
       "Saturated",
       "Exponential growth",
       "Polymer model")
for (i in 1:length(filelist)) {
  toadd <- read.table(paste(filelist[i],".txt",sep=""))
  toadd_h <- read.table(paste(filelist[i],"_hom.txt",sep=""))
  toadd$diff <- (toadd_h[1,]-toadd$V1)/toadd$V1
  toadd$V2 <- samplist[i]
  toadd$cols <- cols[i]
  all_permtest <- rbind(all_permtest,toadd)  
}
customcols <- c(brewercols[1],brewercols[2],"grey50")

pdf(paste(outdir,"/",date,"_hom_prox_violin_1_R.pdf",sep=""),2.3,1.5)
ggplot(all_permtest) + 
  geom_violin(aes(x=factor(V2,levels=samplist),y=diff+1,fill=factor(cols,levels=cols)),alpha=1) + 
  paper.fig.rot + 
  scale_x_discrete(breaks=factor(c("exponential none","exponential cendist","exponential armlen"),levels=samplist)) +
  scale_y_continuous(limits=c(0.75,2.8),expand=c(0,0))+
  xlab("") + ylab("") + theme(axis.text.x=element_blank(),text=element_text(size=7),axis.text.y=element_text(size=7),legend.position=c(0.8,0.9),legend.key.size=unit(0.3,"cm")) + 
  geom_hline(aes(yintercept=1),size=0.25) +
  scale_fill_manual(values = customcols,name="")
dev.off()

# combined genomic homolog proximity plot ----
all_genomic_hom <- data.frame()
subdir <- "/ScSu/32000/none/cen0/arm0.25/genome/"
filelist <- c(paste(outdir,subdir,"2016-12-05","_ILY456_saturated_Sau3AI_Surep-rDNA_min2.txt",sep=""),
              paste(outdir,subdir,"2016-12-05","_ILY456_exponential_Sau3AI_Surep-rDNA_min2.txt",sep=""),
              paste(outdir,subdir,"2016-12-05","_ScSu_v1_none_min2.txt",sep="")
)
samplist <- c("Saturated",
              "Exponential growth",
              "Polymer model")
for (i in 1:length(filelist)) {
  toadd <- read.table(filelist[i],header=FALSE)
  toadd$controls <- samplist[i]
  all_genomic_hom <- rbind(all_genomic_hom,toadd)
}
all_genomic_hom$median <- all_genomic_hom$V4-all_genomic_hom$V23
all_genomic_hom$prop <- all_genomic_hom$V4/all_genomic_hom$median
#all_genomic_hom[all_genomic_hom$prop==Inf,] <- NULL

customcols <- c(brewercols[1],brewercols[2],"grey10")

chrannot <- read.table("nobackup/annotations/ScSu.32000.chr_annotations")
pdf(paste(outdir,"/",date,"_comparative_genomic_hom_prox.pdf",sep=""),3.47,1.3)
ggplot(all_genomic_hom) + paper.fig.rot + geom_bar(aes(x=V2,y=prop,fill=factor(controls,levels=samplist)),width=1,stat="identity",position="identity",alpha=0.5) + 
  scale_x_continuous(breaks=chrannot$start,labels = chrannot$chrbin,expand=c(0,0)) + 
  geom_hline(aes(yintercept=1),size=0.25) + scale_fill_manual(values = customcols,name="") + xlab("") + ylab("") + 
  theme(legend.position=c(0.878,0.95),text=element_text(size=7),legend.key.size=unit(0.3,"cm"),axis.text.x=element_blank()) +
  coord_cartesian(ylim=c(0,17)) + scale_y_continuous(expand=c(0,0))
dev.off()

# combined permtest violin plot #2 ----
all_permtest <- data.frame()
outdir <- "nobackup/figures"
subdir <- "/cen0/arm0.25/permtest/"
excl <- "rDNA-cen"
filelist <- c(paste(outdir,"/ScSu/32000/",excl,subdir,"2016-12-05_ILY456_saturated_Sau3AI_Surep-rDNA_min2",sep=""),
              paste(outdir,"/ScSu/32000/",excl,subdir,"2016-12-05_ILY456_exponential_Sau3AI_Surep-rDNA_min2",sep=""),
              paste(outdir,"/ScSu/32000/",excl,subdir,"2016-12-05_ILY456_nocodazole_Sau3AI_Surep-rDNA_min2",sep=""),
              paste(outdir,"/ScSu/32000/",excl,subdir,"2016-12-05_ScSu_v1_none_min2",sep=""),
              paste(outdir,"/ScSu/32000/",excl,subdir,"2016-12-05_ScSu_v2_none_min2",sep=""),
              paste(outdir,"/ScSu/32000/",excl,subdir,"2016-12-05_ScSu_v3_none_min2",sep=""),
              paste(outdir,"/ScSp/32000/",excl,subdir,"2016-12-05_YMD3264_saturated_Sau3AI_Spgap-rDNA_min2",sep=""),
              paste(outdir,"/ScSp/32000/",excl,subdir,"2016-12-05_YMD3264_exponential_Sau3AI_Spgap-rDNA_min2",sep=""),
              paste(outdir,"/Y12xDBVPG6044/32000/",excl,subdir,"2016-12-05_YMD3271_saturated_Sau3AI_mismapping_min2",sep=""),
              paste(outdir,"/Y12xDBVPG6044/32000/",excl,subdir,"2016-12-05_YMD3271_exponential_Sau3AI_mismapping_min2",sep=""))

samplist <- c("saturated ScSu",
              "exponential ScSu",
              "nocodazole ScSu",
              "model v1",
              "model v2",
              "model v3",
              "saturated ScSp",
              "exponential ScSp",
              "saturated ScSc",
              "exponential ScSc")

cols=c("Saturated",
       "Exponential growth",
       "Nocodazole-arrested",
       "Polymer model",
       "Polymer model",
       "Polymer model",
       "Saturated",
       "Exponential growth",
       "Saturated",
       "Exponential growth")

customcols2 <- c(brewercols[1],brewercols[2],brewercols[3],"grey50")

for (i in 1:length(filelist)) {
  toadd <- read.table(paste(filelist[i],".txt",sep=""))
  toadd_h <- read.table(paste(filelist[i],"_hom.txt",sep=""))
  toadd$diff <- (toadd_h[1,]-toadd$V1)/toadd$V1
  toadd$V2 <- samplist[i]
  toadd$cols <- cols[i]
  all_permtest <- rbind(all_permtest,toadd)  
}

pdf(paste(outdir,"/",date,"_hom_prox_violin_2.pdf",sep=""),3.6,1.5)
ggplot(all_permtest) + 
  geom_violin(aes(x=factor(V2,levels=samplist),y=diff+1,fill=factor(cols,levels=cols)),alpha=1) + 
  paper.fig.rot + 
  xlab("") + ylab("") + theme(axis.text.x=element_blank(),text=element_text(size=7),axis.text.y=element_text(size=7),legend.position=c(0.18,0.9),legend.key.size=unit(0.3,"cm")) + 
  geom_hline(aes(yintercept=1),size=0.25) + 
  scale_fill_manual(values=customcols2,name="") + scale_y_continuous(breaks=c(1,1.5,2))
dev.off()

# combined chr13 violin plot, HAS1 locus deletions ----
all_chr13 <- data.frame()
subdir <- "/ScSu/32000/chr13_interaction_hist/"
filelist <- c(paste(outdir,subdir,"2016-12-05","_ILY456_exponential_Sau3AI_Surep-rDNA",sep=""),
              paste(outdir,subdir,"2016-12-05","_ILY456_saturated_Sau3AI_Surep-rDNA",sep=""),
              paste(outdir,subdir,"2016-12-05","_YMD3266_saturated_Sau3AI_Surep-rDNA",sep=""),
              paste(outdir,subdir,"2016-12-05","_YMD3267_saturated_Sau3AI_Surep-rDNA",sep=""),
              paste(outdir,subdir,"2016-12-05","_YMD3268_saturated_Sau3AI_Surep-rDNA",sep=""),
              paste(outdir,subdir,"2016-12-05","_YMD3269_saturated_Sau3AI_Surep-rDNA",sep="")
)
samplist <- c("exponential",
              "saturated",
              "YMR285-295",
              "HAS1",
              "coding",
              "HAS1-TDA1"
)
samplist2 <- c("exponential",
               "saturated",
               expression(paste(Delta,'NGL2-YMR295C')),
               expression(paste(Delta,'HAS1')),
               expression(paste(Delta,'coding')),
               expression(paste(Delta,'HAS1pr-TDA1pr'))
)

for (i in 1:length(filelist)) {
  toadd <- read.table(paste(filelist[i],"_comp.txt",sep=""))
  toadd$V2 <- samplist[i]
  toadd$V3 <- "comp"
  
  toadd_h <- read.table(paste(filelist[i],"_pair.txt",sep=""))
  toadd_h$V2 <- samplist[i]
  toadd_h$V3 <- "pair"
  all_chr13 <- rbind(all_chr13,toadd)
  all_chr13 <- rbind(all_chr13,toadd_h)
}

pdf(paste(outdir,"/",date,"_chr13_violin_HAS1_deletions.pdf",sep=""),4,1.5)
ggplot(subset(all_chr13,V3=="comp")) + 
  geom_violin(aes(x=factor(V2,levels=samplist),y=V1),fill="darkgrey") + 
  paper.fig.rot + 
  geom_point(data=subset(all_chr13,V3=="pair"),aes(x=factor(V2,levels=samplist),y=V1),color="red",shape="-",size=10)+
  xlab("") + ylab("") + theme(text=element_text(size=7),axis.text.y=element_text(size=7),axis.text.x=element_text(angle=30,hjust=1,size=7)) +
  scale_x_discrete(labels=NULL)
dev.off()

# nup2D violin plot ----

all_chr13 <- data.frame()
subdir <- "/ScSu/32000/chr13_interaction_hist/"
filelist <- c(paste(outdir,subdir,"2016-12-05","_ILY456_exponential_Sau3AI_Surep-rDNA",sep=""),
              paste(outdir,subdir,"2016-12-05","_YMD3377_exponential_Sau3AI_Surep-rDNA",sep=""),
              paste(outdir,subdir,"2016-12-05","_ILY456_galactose_Sau3AI_Surep-rDNA",sep=""),
              paste(outdir,subdir,"2016-12-05","_YMD3377_galactose_Sau3AI_Surep-rDNA",sep=""),
              paste(outdir,subdir,"2016-12-05","_ILY456_saturated_Sau3AI_Surep-rDNA",sep=""),
              paste(outdir,subdir,"2016-12-05","_YMD3377_saturated_Sau3AISurep-rDNA",sep="")
)
samplist <- c("exponential WT",
              "exponential nup2",
              "galactose WT",
              "galactose nup2",
              "stationary WT",
              "stationary nup2"
)

for (i in 1:length(filelist)) {
  toadd <- read.table(paste(filelist[i],"_comp.txt",sep=""))
  toadd$V2 <- samplist[i]
  toadd$V3 <- "comp"
  
  toadd_h <- read.table(paste(filelist[i],"_pair.txt",sep=""))
  toadd_h$V2 <- samplist[i]
  toadd_h$V3 <- "pair"
  all_chr13 <- rbind(all_chr13,toadd)
  all_chr13 <- rbind(all_chr13,toadd_h)
}

pdf(paste(outdir,"/",date,"_nup2D_chr13_violin.pdf",sep=""),3.4,1.5)
ggplot(subset(all_chr13,V3=="comp")) + 
  geom_violin(aes(x=factor(V2,levels=samplist),y=V1),fill="darkgrey") + 
  paperfig + 
  geom_point(data=subset(all_chr13,V3=="pair"),aes(x=factor(V2,levels=samplist),y=V1),color="red",shape="-",size=10)+
  xlab("") + ylab("") + theme(text=element_text(size=7),axis.text.y=element_text(size=7),axis.text.x=element_text(angle=30,hjust=1,size=7)) +
  scale_x_discrete(labels=NULL) + ylim(0,3.5)
dev.off()

# galactose vs. glucose ----
plot_refs <- read.table("plotting/plotting_refs.txt",stringsAsFactors=FALSE,sep = "\t")
colnames(plot_refs) <- c("ref","g1","g2","g1p","g2p","rDNA1","rDNA2","bin1","bin2","c13b1","c13b2")

samp1 <- "ILY456_galactose_Sau3AI"
samp2 <- "ILY456_exponential_Sau3AI"
cond1 <- "gal"
cond2 <- "exp"
ref <- "ScSu"
bsize <- 32000
mask <- "Surep-rDNA"

bins4C <- data.frame(bin=c(16,399),
                     name=c("ScGAL1","SuGAL1"))

# get reference information
g1 <- plot_refs[plot_refs$ref==ref,]$g1
g2 <- plot_refs[plot_refs$ref==ref,]$g2
g1p <- plot_refs[plot_refs$ref==ref,]$g1p
g2p <- plot_refs[plot_refs$ref==ref,]$g2p
nbin1 <- plot_refs[plot_refs$ref==ref,]$bin1
nbin2 <- plot_refs[plot_refs$ref==ref,]$bin2
c13b1 <- plot_refs[plot_refs$ref==ref,]$c13b1
c13b2 <- plot_refs[plot_refs$ref==ref,]$c13b2
rDNA1 <- plot_refs[plot_refs$ref==ref,]$rDNA1
rDNA2 <- plot_refs[plot_refs$ref==ref,]$rDNA2
rDNA <- c(rDNA1,rDNA2)

outdir <- "nobackup/figures"
dir.create(paste(outdir,"/",ref,sep=""))
dir.create(paste(outdir,"/",ref,"/",bsize,sep=""))
outdir <- paste(outdir,ref,bsize,sep="/")

# load ref annotations
chrannot <- read.table(paste("nobackup/annotations/",ref,".",bsize,".chr_annotations",sep=""))
binannot <- read.table(paste("nobackup/annotations/",ref,".",bsize,".bin_annotations",sep=""))

# load data
samp1 <- read.table(paste("nobackup/matrix/",ref,"/",bsize,"/",samp1,".",mask,".matrix",sep=""))
samp2 <- read.table(paste("nobackup/matrix/",ref,"/",bsize,"/",samp2,".",mask,".matrix",sep=""))
samp1[is.na(samp1$V4),]$V4 <- NA
samp2[is.na(samp2$V4),]$V4 <- NA

# combine into single data frame
combined <- merge(samp1,samp2,by=c("V1","V2"))
colnames(combined) <- c("V1","V2","s1raw","s1norm","s2raw","s2norm")

# add annotations
combined.cen <- merge(combined,binannot,by.x="V1",by.y="V1")
colnames(combined.cen) <- c("bin1","bin2","s1raw","s1norm","s2raw","s2norm","chr1","cen1","arm1","chrn1")
combined.cen <- merge(combined.cen,binannot,by.x="bin2",by.y="V1")
colnames(combined.cen) <- c("bin2","bin1","s1raw","s1norm","s2raw","s2norm","chr1","cen1","arm1","chrn1","chr2","cen2","arm2","chrn2")
combined.cen$tel1 <- combined.cen$arm1-combined.cen$cen1
combined.cen$tel2 <- combined.cen$arm2-combined.cen$cen2

# difference heatmap
png(paste(outdir,"/",date,"_",cond1,"-",cond2,"_heatmap.png",sep=""),4.8,2,"in",res=300)
print(ggplot(subset(combined.cen, bin2 >= 385 & bin2 < 452 & bin1 >= 391 & bin1 < 540)) + 
        geom_tile(aes(x=bin1,y=bin2,fill=s1norm-s2norm)) + theme_classic() + 
        paper.fig.rot + hm +
        theme(axis.line = element_blank(), 
              legend.key.size = unit(0.5,"cm"),
              legend.title = element_text(size=7,angle=-90),
              legend.text = element_text(size=7)) + 
        coord_fixed() + 
        geom_hline(data=chrannot,aes(yintercept=start-.5),size=0.25) +
        geom_vline(data=chrannot,aes(xintercept=start-.5),size=0.25) +
        scale_x_continuous(breaks=chrannot$cen,labels=c(chrlabs.all,chrlabs.all),limits=c(391,539),expand=c(0,0)) + 
        scale_y_continuous(breaks=chrannot$cen,labels=c(chrlabs.all,chrlabs.all),limits=c(385,452),expand=c(0,0)) + 
        scale_fill_gradient2(low="royalblue",high="red",limits=c(-1,1),oob=squish,name="\u0394Normalized interaction\nfrequency in galactose",breaks=c(-1,0,1),labels=c("\u2264 -1",0,"\u2265 1"),guide=guide_colorbar(label.hjust=0,title.position="right",title.hjust=0.5)) +
        xlab(expression(paste(italic("S. uvarum"),"centromere")))
)
dev.off()

# box plots of interaction frequencies with GAL1 alleles, aggregated by centromere distance (in bins)

pdf(paste(outdir,"/",date,"_",cond1,"-",cond2,"_Suva_GAL1_byCEN.pdf",sep=""),2.4,2)
print(ggplot(subset(combined.cen,bin1==399 & chr1 != chr2)) + 
        geom_boxplot(aes(x=factor(pmin(cen2,8)),y=s1norm-s2norm),outlier.size=0.5) + 
        geom_hline(aes(yintercept=0),size=0.5)+
        scale_x_discrete(expand=c(0.05,0),labels=c(0,1,2,3,4,5,6,7,expression(phantom(x) >=8))) + 
        paper.fig + coord_cartesian(ylim=c(-2,2)) + 
        theme(legend.position="none") + xlab("") + ylab(""))
dev.off()

pdf(paste(outdir,"/",date,"_",cond1,"-",cond2,"_Scer_GAL1_byCEN.pdf",sep=""),2.4,2)
print(ggplot(subset(combined.cen,bin1==16 & chr1 != chr2)) + 
        geom_boxplot(aes(x=factor(pmin(cen2,8)),y=s1norm-s2norm),outlier.size=0.5) + 
        geom_hline(aes(yintercept=0),size=0.5)+
        scale_x_discrete(expand=c(0.05,0),labels=c(0,1,2,3,4,5,6,7,expression(phantom(x) >=8))) + 
        paper.fig + coord_cartesian(ylim=c(-2,2)) + 
        theme(legend.position="none") + xlab("") + ylab(""))
dev.off()

# Wilcoxon rank test for box plots
for (bin in c(16,399)){
  for (w in 0:7) {
    sub <- subset(combined.cen,bin1==bin & chr1 != chr2 & cen2>=8)
    h <- wilcox.test(sub$s1norm-sub$s2norm)
    print(h$p.value)
  }
}

# scatter plot of interaction frequencies in galactose vs. glucose
pdf(paste(outdir,"/",date,"_",cond1,"_v_",cond2,"_GAL1_scatter.pdf",sep=""),3.6,3)
print(ggplot(subset(combined.cen[order(combined.cen$bin1==16 & combined.cen$bin2==399),],cen1==1 & cen2==1 & chr1 != chr2)) + geom_point(aes(x=s1norm,y=s2norm,color=!(bin1==16 & bin2==399)),size=1) + scale_color_manual(values=c("red","grey"),labels=c(expression(italic("GAL1-GAL1")),expression("Similar interactions")),name=c("")) + paper.fig + theme(legend.position=c(0.18,0.8),legend.text.align=0,legend.key=element_blank(),plot.margin=unit(c(0,0,0,0),"cm")) + coord_fixed() + xlab(paste("Normalized interaction frequency in ","overnight galactose",sep="")) + ylab(paste("Normalized interaction frequency in ","glucose",sep="")))
dev.off()

# "4C" like plots of interaction frequencies with each GAL1 allele, zoomed into the homologous chromosome
pdf(paste(outdir,"/",date,"_",cond1,"_v_",cond2,"_","ScGAL1","_4C_zoom.pdf",sep=""),3.5,2)
print(ggplot(subset(combined,V1==16)) + 
        geom_line(aes(x=V2,y=s1norm),color="red") +
        geom_line(aes(x=V2,y=s2norm),color="blue") + 
        scale_x_continuous(limits=c(391,432),breaks=seq(391,432,5),labels=seq(0,41*32,5*32),expand=c(0,0)) + 
        scale_y_continuous(expand=c(0,0)) +
        coord_cartesian(ylim=c(0,max(c(subset(combined.cen,bin1==16 & chr1 != chr2 & (cen2 == 0 | tel2 > 0))$s1norm,subset(combined.cen,bin1==16 & chr1 != chr2 & (cen2 == 0 | tel2 > 0))$s2norm),na.rm = T)+1)) + 
        theme_classic() +
        geom_vline(data=chrannot,aes(xintercept=cen),linetype=3) + 
        paper.fig.rot + xlab(expression(paste(italic("S. uvarum")," chromosome IV (kb)"))) + ylab("Normalized interaction frequency")
)
dev.off()

pdf(paste(outdir,"/",date,"_",cond1,"_v_",cond2,"_","SuGAL1","_4C_zoom.pdf",sep=""),3.5,2)
print(ggplot(subset(combined,V1==399)) + 
        geom_line(aes(x=V2,y=s1norm),color="red") +
        geom_line(aes(x=V2,y=s2norm),color="blue") + 
        scale_x_continuous(limits=c(8,34),breaks=seq(8,34,5),labels=seq(0,26*32,5*32),expand=c(0,0)) + 
        scale_y_continuous(expand=c(0,0)) +
        coord_cartesian(ylim=c(0,max(c(subset(combined.cen,bin1==399 & chr1 != chr2 & (cen2 == 0 | tel2 > 0))$s1norm,subset(combined.cen,bin1==399 & chr1 != chr2 & (cen2 == 0 | tel2 > 0))$s2norm),na.rm = T)+1)) + 
        theme_classic() +
        geom_vline(data=chrannot,aes(xintercept=V4),linetype=3) + 
        paper.fig.rot + xlab(expression(paste(italic("S. cerevisiae")," chromosome II (kb)"))) + ylab("Normalized interaction frequency")
)
dev.off()