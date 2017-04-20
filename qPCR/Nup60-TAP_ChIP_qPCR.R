# Nup60-TAP_ChIP_qPCR.R
# plotting script for ChIP qPCR data
# Seungsoo Kim

# load libraries
library(ggplot2)
library(dplyr)
library(reshape)
library(RColorBrewer)
library(grid)
library(scales)
library(stringr)

# colors
brewercols <- brewer.pal(4, "Set1")
cols = brewercols[c(2, 4, 1)]

# load data
data1 <- read.table("2017-03-25_ViiA7.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
data2 <- read.table("2017-03-31_ViiA7.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
data <- rbind(data1, data2)

# group technical replicates together and average quantity
data$sample <- paste(data$Sample.Name, data$Target.Name, sep = " ")
data <- select(data, sample, Quantity)
data_grouped <- group_by(data, sample)
averaged <- summarize(data_grouped, Quantity.Mean = sum(Quantity))

# add columns for sample characteristics
samp.data <- data.frame(str_split_fixed(averaged$sample, ' ', 4))
averaged$cond <- samp.data$X1
averaged$rep <- samp.data$X2
averaged$IP <- samp.data$X3
averaged$target <- samp.data$X4
averaged <- select(averaged, target, cond, rep, IP, Quantity.Mean)

# calculate IP/input ratios
casted <- cast(averaged, target + cond + rep ~ IP)
casted$IPratio <- casted$Ab/casted$input
casted$mockratio <- casted$BSA/casted$input
casted <- select(casted, cond, target, rep, IPratio, mockratio)
tempmelt <- melt(casted, id.vars = 1:3, measure.vars = c("IPratio", "mockratio"))
casted.ratios <- cast(tempmelt, cond + IP + rep ~ target)

casted.ratios$GAL1pr <- casted.ratios$GAL1pr/casted.ratios$PRM1
casted.ratios$TDA1_1 <- casted.ratios$TDA1_1/casted.ratios$PRM1
casted.ratios$TDA1_2 <- casted.ratios$TDA1_2/casted.ratios$PRM1
casted.ratios$TDA1_3 <- casted.ratios$TDA1_3/casted.ratios$PRM1
casted.ratios$TDA1pr_up <- casted.ratios$TDA1pr_up/casted.ratios$PRM1
casted.ratios$TDA1pr_dn <- casted.ratios$TDA1pr_dn/casted.ratios$PRM1
casted.ratios$PRM1 <- casted.ratios$PRM1/casted.ratios$PRM1

remelted <- melt(casted.ratios, id.vars = 1:3)

# group samples to calculate p-values
byrep <- cast(remelted, IP + target ~ cond + rep)
byrep <- byrep[byrep$target != "PRM1", ]
t.result <- apply(byrep[, -(1:2)], 1, function(x) t.test(x[1:3], x[4:6]))
byrep$pval <- unlist(lapply(t.result, function(x) x$p.value))

# calculate mean and sd
summarized <- cast(remelted, cond + IP + target ~ ., c(mean, sd))

# positions for bar plot
dodge <- position_dodge(width = 0.9)

# plot antibody data
plotdata <- subset(summarized, target != "PRM1" & IP == "Ab")
plotdata$mark <- c("*", "", "", "", "", "", "", "", "", "", "", "")
plotdata$target <- factor(plotdata$target, levels = c("GAL1pr", "TDA1pr_up", "TDA1pr_dn", "TDA1_1", "TDA1_2", "TDA1_3"))
plotdata$cond <- factor(plotdata$cond, levels = c("glucose", "galactose"))

pdf("Nup60-TAP_ChIP_normed_Ab.pdf", 2.5, 1.6)
ggplot(plotdata) + 
 geom_errorbar(aes(x = target, ymin = mean - sd/sqrt(3), ymax = mean + sd/sqrt(3), fill = cond), stat = "identity", position = dodge, width = 0.5) + 
 geom_bar(aes(x = target, y = mean, fill = cond), color = "black", stat = "identity", position = dodge) + 
 scale_fill_manual(values = cols, name = "", labels = c("Glucose", "Galactose")) + 
  theme_classic() + scale_y_continuous(limits = c(0, max(byrep$mean + byrep$sd/sqrt(3))*1.15), expand = c(0, 0)) + 
  xlab("") + ylab("Nup60-TAP IP/input\nnormalized to PRM1") + 
 scale_x_discrete(labels = c("GAL1pr", "prom 1", "prom 2", "CDS 1", "CDS 2", "CDS 3")) +
 theme(text = element_text(size = 7), legend.position = c(0.8, 0.99), axis.ticks.x = element_blank(), legend.key.size = unit(0.25, "cm")) +
 geom_text(aes(label = mark, x = target, fill = cond, y = mean + sd/sqrt(3) + 0.15), size = 6, position = dodge)
dev.off()

# plot BSA mock-IP data
plotdata <- subset(summarized, target != "PRM1" & IP == "BSA")
plotdata$mark <- c("", "", "", "", "", "", "", "", "", "", "", "")
plotdata$target <- factor(plotdata$target, levels = c("GAL1pr", "TDA1pr_up", "TDA1pr_dn", "TDA1_1", "TDA1_2", "TDA1_3"))
plotdata$cond <- factor(plotdata$cond, levels = c("glucose", "galactose"))

pdf("Nup60-TAP_ChIP_normed_BSA.pdf", 3.3, 1.6)
ggplot(plotdata) + 
 geom_errorbar(aes(x = target, ymin = mean - sd/sqrt(3), ymax = mean + sd/sqrt(3), fill = cond), stat = "identity", position = dodge, width = 0.5) + 
 geom_bar(aes(x = target, y = mean, fill = cond), color = "black", stat = "identity", position = dodge) + 
 scale_fill_manual(values = cols, name = "", labels = c("Glucose", "Galactose")) + 
  theme_classic() + scale_y_continuous(limits = c(0, max(byrep$mean + byrep$sd/sqrt(3))*1.15), expand = c(0, 0)) + 
  xlab("") + ylab("Nup60-TAP IP/input\nnormalized to PRM1") + 
 scale_x_discrete(labels = c("GAL1pr", "prom 1", "prom 2", "CDS 1", "CDS 2", "CDS 3")) +
 theme(text = element_text(size = 7), legend.position = c(0.8, 0.99), axis.ticks.x = element_blank(), legend.key.size = unit(0.25, "cm")) +
 geom_text(aes(label = mark, x = target, fill = cond, y = mean + sd/sqrt(3) + 0.15), size = 6, position = dodge)
dev.off()
