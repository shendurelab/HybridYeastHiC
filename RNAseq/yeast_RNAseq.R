# yeast_RNAseq.R
# plotting and DESeq2 analysis of RNA-seq data
# Seungsoo Kim

# load libraries ----
library(ggplot2)
library(dplyr)
library(reshape)
library(RColorBrewer)
library(grid)
library(gplots)
library(scales)
library(DESeq2)

# get gene lengths
genelens <- read.table("nobackup/references/sacCer3_genes_length.txt")
colnames(genelens) <- c("gene", "length")

# samples
samps <- c("asy_r1", "asy_r2", "asy_r3", "gal_r1", "gal_r2", "gal_r3", "sat_r1", "sat_r2", "sat_r3")

# load data
counts <- read.table("nobackup/counts/all.counts.txt")
colnames(counts) <- c("gene", samps)

# add gene lengths to data, and exclude non-uniquely assigned read counts
data <- merge(genelens, counts)

# extract only raw counts for normalization
rawcounts <- data[, -(1:2)]
rownames(rawcounts) <- data$gene

# normalize counts
normed <- sweep(rawcounts, 2, colSums(rawcounts), "/") # by total read count per sample
normed <- sweep(normed, 1, genelens$length, "/") # by gene length (in bp)
normed <- normed*1000000000 # multiply by 1 billion to make FPKM
normwgenes <- cbind(data$gene, normed) # add back gene annotations
colnames(normwgenes) <- c("gene", samps)

# colors
brewercols <- brewer.pal(4, "Set1")
cols=brewercols[c(2, 4, 1)]

# condition names
conds <- c("Glucose", "Galactose", "Saturated")

# function for calculating number of * to add based on p-value
# takes a vector of increasingly stringent (lower) p-value cutoffs and 
# outputs a vector of strings, each with the appropriate number of asterisks
stars <- function(thresh, pval) {
 n = 0
 for (t in thresh) {
  if (pval < t) {
   n = n + 1
  }
 }
 return(paste(rep("*", n), collapse = ""))
}

thresh = c(0.05, 0.01, 0.001, 0.0001) # p-value thresholds for asterisks

# loop through genes
genesofinterest <- c("YMR290C", "YMR291W")
genelabs <- c("YMR290C (HAS1)", "YMR291W (TDA1)")
for (i in 1:length(genesofinterest)) {
 geneofinterest <- genesofinterest[i]
 genelab <- genelabs[i]
 
 # subset relevant data
 sub <- subset(normwgenes, gene == geneofinterest)
 table <- data.frame
 table <- sub[, 2:4]
 table[2, ] <- sub[, 5:7]
 table[3, ] <- sub[, 8:10]
 colnames(table) <- c("r1", "r2", "r3")
 table$ave <- (table$r1 + table$r2 + table$r3)/3
 table$sd <- apply(table[, 1:3], 1, sd)
 rownames(table) <- conds
 table$cond <- factor(conds, levels = conds)

 # calculate p-values 
 gal.test <- t.test(table[1, 1:3], table[2, 1:3])
 gal.stars <- stars(thresh, gal.test$p.value)
 sat.test <- t.test(table[1, 1:3], table[3, 1:3])
 sat.stars <- stars(thresh, sat.test$p.value)
 table$stars <- c("", gal.stars, sat.stars)
 
 # bar plot of FPKM of one gene in each condition, with error bars (SEM) and asterisks
 pdf(paste("2017-04-12_fpkm_", geneofinterest, ".pdf", sep=""), 1.8, 1.5)
 print(ggplot(table) + 
     geom_text(aes(x = cond, y = (ave + sd) + max(table$ave)*.05, label = stars)) + 
     geom_errorbar(aes(x = cond, ymin = ave - sd/sqrt(3), ymax = ave + sd/sqrt(3)), width = .2) + 
     geom_bar(aes(x = cond, fill = cond, y = ave), stat = "identity", color = "black") + 
     theme_classic() + scale_y_continuous(limits = c(0, 1.2*max(table$ave + table$sd)), expand = c(0, 0)) + 
     scale_fill_manual(values = cols) + 
     theme(plot.margin = unit(c(0.2, 0, 0, 0), "cm"), 
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5, face = "italic"), 
        text = element_text(size = 7)) + 
     xlab("") + ylab("FPKM") + ggtitle(genelab))
 dev.off()
}

# DESeq2 ----

# sample table
samples <- read.table("sample_table.txt", header=TRUE, stringsAsFactors = TRUE)

# run DESeq2
dds <- DESeqDataSetFromMatrix(countData = rawcounts, 
               colData = samples, 
               design = ~ condition)
dds <- dds[ rowSums(counts(dds)) > 1, ] # filter out genes with <= 1 fragment
dds <- DESeq(dds)

# calculate p-values and log2 Fold Change for galactose vs. glucose
gal <- data.frame(results(dds, contrast = c("condition", "gal", "asy")))
gal$gene <- rownames(gal)

# calculate p-values and log2 Fold Change for saturated vs. glucose
sat <- data.frame(results(dds, contrast = c("condition", "sat", "asy")))
sat$gene <- rownames(sat)

# histograms of log2 Fold Change, with GAL1 (for galactose only), HAS1, and TDA1 labeled

pdf("2017-04-12_gal_foldChange_hist.pdf", 5, 1.5)
ggplot(gal) + geom_histogram(aes(x = log2FoldChange), binwidth = 0.2, fill = cols[2]) + 
 geom_vline(aes(xintercept = gal[gal$gene == "YBR020W", ]$log2FoldChange)) + 
 geom_vline(aes(xintercept = gal[gal$gene == "YMR290C", ]$log2FoldChange)) + 
 geom_vline(aes(xintercept = gal[gal$gene == "YMR291W", ]$log2FoldChange)) + 
 theme_classic() + theme(text = element_text(size = 7), legend.position = "none") + 
 scale_y_continuous(expand = c(0, 0)) + 
 ylab("Genes") + xlab(expression(paste(log[2], " Fold Change in Galactose vs. Glucose")))
dev.off()

pdf("2017-04-12_sat_foldChange_hist.pdf", 5, 1.5)
ggplot(sat) + geom_histogram(aes(x = log2FoldChange), binwidth = 0.2, fill = cols[3]) + 
 geom_vline(aes(xintercept = sat[sat$gene == "YMR290C", ]$log2FoldChange)) + 
 geom_vline(aes(xintercept = sat[sat$gene == "YMR291W", ]$log2FoldChange)) + 
 theme_classic() + theme(text = element_text(size = 7), legend.position = "none") + 
 scale_y_continuous(expand = c(0, 0)) + 
 ylab("Genes") + xlab(expression(paste(log[2], " Fold Change in Saturated vs. Glucose")))
dev.off()
