library(vcfR)
library(reshape2)
library(pinfsc50)
library(reshape2)
library(ggplot2)
library(ape)
library(stringr)
vcf_file <- "data/final_all_merged.vcf"
# read vcf
nes_vcf <- read.vcfR(vcf_file, verbose = FALSE )
dp <- extract.gt(nes_vcf, element = "DP", as.numeric = TRUE)

# Reorganize and render violin plots.
dpf <- melt(dp, varnames=c("Index", "Sample"), value.name = "Depth", na.rm=TRUE)
dpf <- dpf[ dpf$Depth > 0,]
p <- ggplot(dpf, aes(x=Sample, y=Depth)) + geom_violin(fill="#C0C0C0", adjust=1.0,scale = "count", trim=TRUE)

p <- p + theme_bw()
p <- p + ylab("Read Depth (DP)")
p <- p + theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 60, hjust = 1))
p <- p + stat_summary(fun.data=mean_sdl, geom="pointrange", color="black")
p <- p + scale_y_continuous(trans=scales::log2_trans(), breaks=c(1, 10, 100, 1000))
p
# heatmap
heatmap.bp(dp, rlabel = FALSE)



gt <- extract.gt(nes_vcf, element = "gt")
# transpose and data.frame
#gt <- as.data.frame(t(gt), stringsAsFactors = FALSE)
gt <- as.data.frame(gt, stringsAsFactors = FALSE)
ind_names <- names(gt)
gt <- t(gt)
row.names(gt) <- ind_names
# NA handling
gt[gt == "./."] <- NA
# split columns
snp_geno <- do.call(cbind, apply(gt, 2, function(x) colsplit(x, "/", c("a","b"))))
# convert
snp_genotypes <- inbreedR::convert_raw(snp_geno)
# check data
inbreedR::check_data(snp_genotypes)

library(inbreedR)
nes <- snp_genotypes

# some checks
typed <- rowSums(!is.na(nes))
hist(typed)
het <- MLH(nes)
plot(het, typed)


hist(het)
g2_snps(snp_genotypes)

# mouse
library(inbreedR)
data("mouse_snps")
mouse_snps
typed <- rowSums(!is.na(mouse_snps))
het <- MLH(mouse_snps)
plot(het, typed)


# filter SNPs that were typed in at least 90% of individuals
test <- which(colSums(!is.na(nes)) > 90)

nes1 <- nes[test]

typed <- rowSums(!is.na(nes1))
het <- MLH(nes1)
hist(het)
g2_snps(nes1, nboot = 100)
plot(het, typed)


# total number of sequences per individual against heterozygosity
num_reads <- read.table("data/number_of_reads_per_ind.txt")

plot(num_reads$V1/4, het)

