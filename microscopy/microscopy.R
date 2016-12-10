# load libraries ----
library(ggplot2)
library(scales)
library(RColorBrewer)
library(dplyr)
library(data.table)

# color palette and condition names ----
brewercols <- brewer.pal(7,"Set1")
cols=brewercols[c(2,4,1)]
conds=c("Glucose","Galactose","Saturated")

# HAS1 allele clustering ----

# load data
data <- read.csv("clustering_data.csv",header=FALSE)

# frequency polygon plot
pdf("HAS1_clustering_freqpoly.pdf",2.2,1.7)
ggplot(data) + geom_freqpoly(aes(x=V1,color=factor(V2,levels=conds)),binwidth=0.15,boundary=-0.001,size=.5) + 
  scale_color_manual(values = cols,name="",labels=conds) + 
  theme_classic() + 
  scale_x_continuous(expand=c(0,0)) + 
  theme(legend.position=c(0.55,0.98),text=element_text(size=7),axis.text=element_text(size=7),legend.key.size=unit(0.25,"cm")) + 
  xlab(expression(paste(italic("HAS1"),"-",italic("HAS1")," distance (",mu,"m)",sep=""))) + ylab("Percent of cells") +
  theme(plot.margin=unit(c(0.2,0.2,0,0),"cm"))
dev.off()

# calculate frequency d < 0.55 µm
dt <- data.table(data)
dt.out <- dt[,list(prop=sum(V1 < .55)/sum(V1 > 0)), by=V2]

# clustering bar plot
pdf("HAS1_clustering_bar.pdf",1.3,2)
ggplot(data.table(dt.out)) + geom_bar(aes(x=factor(V2,levels=conds),y=100*prop,fill=factor(V2,levels=conds)),stat="identity",color="black") + scale_fill_manual(values = cols,name="",labels=conds) +
  theme_classic() + xlab("") + theme(axis.ticks.x=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),text=element_text(size=7),axis.text=element_text(size=7),legend.position="none",legend.key.size=unit(0.25,"cm"),plot.margin=unit(c(0.1,0.1,0,0),"cm")) + 
  scale_y_continuous(expand=c(0,0),limits=c(0,32)) + ylab("\n") +
  scale_x_discrete(labels=conds) + 
  geom_segment(aes(x=1,xend=2,y=27.5,yend=27.5),size=0.25) + geom_segment(aes(x=1,xend=3,y=31.5,yend=31.5),size=0.25)
dev.off()

# HAS1 peripheral localization ----

# load data, calculate summary statistics
data2 <- read.csv("peripheral_data.csv",header=TRUE)
data2mat <- as.matrix(data2[,-1:-2])
data2$ave <- apply(data2mat,1,mean)
data2$sd <- apply(data2mat,1,sd)
data2$sem <- data2$sd/sqrt(3)

# plotting parameters
dodge <- position_dodge(width=0.8)
wid=0.8
spac=0.4
data2$x <- c(3*wid+spac,0,5*wid+2*spac,7*wid+3*spac,9*wid+4*spac,
             4*wid+spac,wid,6*wid+2*spac,8*wid+3*spac,10*wid+4*spac,
             2*wid)
data2$comb <- paste(data2$strain,data2$cond,sep="-")

# strain names and labels
strains <- c("WTalpha","WTa","mlp2","nup2","nup100")
strainlabs <- c(expression(italic(MAT * alpha)),
                expression(italic("MATa")),
                expression(italic(mlp2 * Delta)),
                expression(italic(nup2 * Delta)),
                expression(italic(nup100) * Delta))

# statistical significance
data2$mark <- c("","","","","","*","*","*","","*","*")

# peripheral localization bar plot
pdf("HAS1_peripheral.pdf",3.7,1.7)
ggplot(data2) + 
  geom_errorbar(aes(x=x,ymin=ave-sem,ymax=ave+sem),width=0.5) +
  geom_bar(aes(x=x,fill=factor(cond,levels=conds),y=ave),width=wid,stat="identity",color="black") + 
  scale_fill_manual(values=cols,name="",labels=conds) +
  theme_classic() + 
  scale_y_continuous(expand=c(0,0),limits=c(0,89)) + 
  ylab("Peripheral localization\n(% of cells)") + xlab("") +
  scale_x_continuous(breaks=c(wid,3.5*wid+spac,5.5*wid+2*spac,7.5*wid+3*spac,9.5*wid+4*spac),labels=strainlabs) + 
  theme(text=element_text(size=7),axis.ticks.x=element_blank(),axis.text=element_text(size=7),legend.key.size=unit(0.25,"cm"),
        plot.margin=unit(c(0.2,0,0,0),"cm"),legend.position=c(0.742,0.9)) +
  geom_text(aes(label=mark,x=x,y=ave+sem+3),size=6)
dev.off()
