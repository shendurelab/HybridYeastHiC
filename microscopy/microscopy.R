# microscopy.R
# plotting of live-cell microscopy data
# Seungsoo Kim

# load libraries ----
library(ggplot2)
library(scales)
library(RColorBrewer)
library(dplyr)
library(data.table)
library(reshape)

# color palette and condition names ----
brewercols <- brewer.pal(4, "Set1")
cols = brewercols[c(2, 4, 1)]
conds = c("Glucose", "Galactose", "Saturated")

# HAS1 allele clustering ----

# load data
clustering <- read.csv("clustering_data.csv", header = FALSE)

# frequency polygon plot
pdf("HAS1_clustering_freqpoly.pdf", 2.2, 1.7)
ggplot(clustering) + 
  geom_freqpoly(aes(x = V1, color = factor(V2, levels = conds)), binwidth = 0.15, size = .5) + 
  scale_color_manual(values = cols, name = "", labels = conds) + 
  theme_classic() + 
  scale_x_continuous(expand = c(0, 0)) + 
  theme(legend.position = c(0.55, 0.98),
        text = element_text(size = 7),
        axis.text = element_text(size = 7),
        legend.key.size = unit(0.25, "cm")) + 
  xlab(expression(paste("Distance between ", italic("HAS1"), "-", italic("TDA1"), " alleles (", mu, "m)", sep = ""))) + ylab("Percent of cells") +
  theme(plot.margin = unit(c(0.2, 0.2, 0, 0), "cm"))
dev.off()

# calculate frequency d < 0.55 µm
clustering.dt <- data.table(clustering)
clustering.dt.out <- clustering.dt[, list(prop = sum(V1 < .55)/sum(V1 > 0)), by = V2]

# clustering bar plot
pdf("HAS1_clustering_bar.pdf", 1.3, 2)
ggplot(data.table(clustering.dt.out)) + 
  geom_bar(aes(x = factor(V2, levels = conds), y = 100*prop, fill = factor(V2, levels = conds)), stat = "identity", color = "black") + 
  scale_fill_manual(values = cols, name = "", labels = conds) +
  theme_classic() + xlab("") + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 7),
        axis.text = element_text(size = 7),
        legend.position = "none",
        legend.key.size = unit(0.25, "cm"),
        plot.margin = unit(c(0.1, 0.1, 0, 0), "cm")) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 32)) + ylab("\n") +
  scale_x_discrete(labels = conds) + 
  geom_segment(aes(x = 1, xend = 2, y = 27.5, yend = 27.5), size = 0.25) + geom_segment(aes(x = 1, xend = 3, y = 31.5, yend = 31.5), size = 0.25)
dev.off()

# haploid HAS1 peripheral localization ----

# load data, calculate summary statistics
hap.peripheral <- read.csv("haploid_peripheral.csv", header = TRUE)
hap.peripheral$peripheral <- hap.peripheral$ON/hap.peripheral$total*100
hap.peripheral <- select(hap.peripheral, strain, condition, replicate, peripheral)
hap.peripheral <- dcast(hap.peripheral, strain + condition ~ replicate)

hap.peripheral.mat <- as.matrix(hap.peripheral[, 3:5])
hap.peripheral$ave <- apply(hap.peripheral.mat, 1, mean)
hap.peripheral$sd <- apply(hap.peripheral.mat, 1, sd)
hap.peripheral$sem <- hap.peripheral$sd/sqrt(3)

# sort data by strain, then condition
strains <- c("WTalpha", "WTa", "mlp2D", "nup2D", "nup100D")
hap.peripheral$strain <- factor(hap.peripheral$strain, levels = strains)
hap.peripheral$condition <- factor(hap.peripheral$condition, levels = conds)
hap.peripheral <- hap.peripheral[order(hap.peripheral$strain, hap.peripheral$condition), ]

# strain names and labels
strainlabs <- c(expression(MAT * alpha), 
        expression("MATa"), 
        expression(mlp2 * Delta), 
        expression(nup2 * Delta), 
        expression(nup100 * Delta))

# plotting parameters
dodge <- position_dodge(width = 0.8)
wid = 0.8
spac = 0.4
hap.peripheral$x <- c(wid*(0:2), spac + wid*(3:4), 2*spac + wid*(5:6), 3*spac + wid*(7:8), 4*spac + wid*(9:10))
hap.peripheral$comb <- paste(hap.peripheral$strain, hap.peripheral$cond, sep = "-")

# statistical significance
hap.peripheral$mark <- c("", "*", "*", "", "*", "", "*", "", "", "", "*")

# peripheral localization bar plot
pdf("HAS1_peripheral.pdf", 3, 1.7)
ggplot(hap.peripheral) + 
  geom_errorbar(aes(x = x, ymin = ave-sem, ymax = ave+sem), width = 0.5) +
  geom_bar(aes(x = x, fill = factor(condition, levels = conds), y = ave), width = wid, stat = "identity", color = "black") + 
  scale_fill_manual(values = cols, name = "", labels = conds) +
  theme_classic() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 89)) + 
  ylab("Peripheral localization\n(% of cells)") + xlab("") +
  scale_x_continuous(breaks = c(wid, 3.5*wid+spac, 5.5*wid+2*spac, 7.5*wid+3*spac, 9.5*wid+4*spac), labels = strainlabs) + 
  theme(text = element_text(size = 7),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 7),
        legend.key.size = unit(0.25, "cm"),
        plot.margin = unit(c(0.2, 0, 0, 0), "cm"),
        legend.position = c(0.742, 0.9)) +
  geom_text(aes(label = mark, x = x, y = ave+sem+3), size = 6)
dev.off()

# diploid HAS1 peripheral localization ----

# diploid HAS1 joint peripheral localization and pairing analysis ----
dip.reanalysis <- read.csv("diploid_joint_peripheral_pairing.csv", header = TRUE)
dip.reanalysis$pairing <- dip.reanalysis$paired/(dip.reanalysis$paired + dip.reanalysis$unpaired)*100
dip.reanalysis$total <- dip.reanalysis$paired + dip.reanalysis$unpaired

peripheral_labels <- c("ON-ON", "ON-OFF", "OFF-OFF")

pdf("diploid_peripheral.pdf", 2.1, 1.7)
ggplot(dip.reanalysis) + geom_bar(aes(x = factor(peripheral, levels = peripheral_labels), y = total, fill = factor(condition, levels = conds)), stat = "identity", position = "dodge", color = "black") + 
 theme_classic() + 
 scale_y_continuous(expand = c(0, 0), limits = c(0, 55)) + 
 scale_fill_manual(values = cols, name = "") +
 theme(text = element_text(size = 7),
       axis.ticks.x = element_blank(),
       axis.text = element_text(size = 7),
       legend.key.size = unit(0.25, "cm"),
       legend.position = c(0.6, 0.96),
       plot.margin = unit(c(0.2, 0, 0, 0), "cm")) +
 xlab("Peripheral localization") + ylab("Percent of cells")
dev.off()

# Chi-squared test
chisq.test(dip.reanalysis[, -1])
dip.reanalysis.table <- dcast(select(dip.reanalysis, condition, peripheral, total), condition ~ peripheral)
chisq.test(dip.reanalysis.table[, -1])

# joint analysis
dip.reanalysis.pairing <- melt(dip.reanalysis[, c("condition", "peripheral", "pairing")], value.name = "pairing")

pdf("diploid_peripheral_pairing.pdf", 2.1, 1.7)
ggplot(dip.reanalysis.pairing) + 
  geom_bar(aes(x = factor(peripheral, levels = peripheral_labels), y = value, fill = factor(condition, levels = conds)), stat = "identity", position = "dodge", color = "black") + 
  theme_classic() + 
  scale_x_discrete() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 55)) + 
  scale_fill_manual(values = cols, name = "") +
  theme(text = element_text(size = 7),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 7),
        legend.key.size = unit(0.25, "cm"),
        legend.position = c(0.6, 0.96),
        plot.margin = unit(c(0.2, 0, 0, 0), "cm")) +
  xlab("Peripheral localization") + ylab("Percent clustered (< 0.55 um)")
dev.off()

# cell cycle peripheral localization ----

# load data
cell.cyc.peripheral <- read.csv("cell_cycle_peripheral.csv", header = TRUE)

# calculate percentages of cells with peripheral localization
cell.cyc.peripheral$percent <- cell.cyc.peripheral$peripheral/cell.cyc.peripheral$total*100
cell.cyc.peripheral <- select(cell.cyc.peripheral, cell.cycle, condition, replicate, percent)

# group by replicate
cell.cyc.peripheral <- dcast(cell.cyc.peripheral, cell.cycle + condition ~ replicate)

# calculate summary statistics
cell.cyc.peripheral.mat <- as.matrix(cell.cyc.peripheral[, 3:5])
cell.cyc.peripheral$ave <- apply(cell.cyc.peripheral.mat, 1, mean)
cell.cyc.peripheral$sd <- apply(cell.cyc.peripheral.mat, 1, sd)
cell.cyc.peripheral$sem <- cell.cyc.peripheral$sd/sqrt(3)

# sort by cell cycle and condition
phases = c("G1", "S", "G2/M")
cell.cyc.peripheral$cell.cycle <- factor(cell.cyc.peripheral$cell.cycle, levels = phases)
cell.cyc.peripheral$condition <- factor(cell.cyc.peripheral$condition, levels = conds)
cell.cyc.peripheral <- cell.cyc.peripheral[order(cell.cyc.peripheral$cell.cycle, cell.cyc.peripheral$condition), ]

# asterisks for statistical significance
cell.cyc.peripheral$mark <- c("", "*", "", "*", "", "*")

dodge = position_dodge(width = 0.9)

pdf("cell_cycle_peripheral.pdf", 2, 1.6)
ggplot(cell.cyc.peripheral) + 
  geom_errorbar(aes(x = cell.cycle, fill = condition, ymin = ave - sem, ymax = ave + sem), width = 0.5, position = dodge) +
  geom_bar(aes(x = cell.cycle, fill = condition, y = ave), width = 0.9, stat = "identity", color = "black", position = dodge) + 
  scale_fill_manual(values = cols, name = "", labels = conds) +
  theme_classic() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 89)) + 
  ylab("Peripheral localization\n(% of cells)") + xlab("") +
  theme(text = element_text(size = 7),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 7),
        legend.key.size = unit(0.25, "cm"),
        legend.position = c(0.8, 0.99),
        plot.margin = unit(c(0.2, 0, 0, 0), "cm")) +
  geom_text(aes(label = mark, x = cell.cycle, fill = condition, y = ave+sem+3), size = 6, position = dodge)
dev.off()

# calculate p-values
cell.cyc.peripheral.array <- array(data = as.vector(melt(cell.cyc.peripheral[3:5])$value), dim = c(2, 3, 3))
ttest.pval <- function(x) {
 t.test(x[1, ], x[2, ])$p.value
}
apply(cell.cyc.peripheral.array, 2, ttest.pval)
