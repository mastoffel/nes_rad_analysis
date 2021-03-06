# nes rad vcf file last filtering steps
library(vcfR)
library(reshape2)
library(ape)
library(stringr)
# vcf_file <- "data/filtered_snps.recode.vcf"
vcf_file <- "data/final_all_merged.vcf"
# read vcf
nes <- read.vcfR(vcf_file, verbose = FALSE )
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

# nes2 <- read.vcfR(vcf_file2, verbose = FALSE )
# summary
nes
head(nes)
plot(nes)

# chromR format
chrom <- create.chromR(name = 'supernes', vcf=nes)
head(chrom)
plot(chrom)
# fixed <- chrom@vcf@fix
chromoqc(chrom, dp.alpha = 22)

# AN total number of alleles in called genotypes
all_AN <- as.numeric(extract.info(chrom, "AN"))
all_AF <- as.numeric(extract.info(chrom, "AF"))
all_DP <- as.numeric(extract.info(chrom, "DP"))
hist(all_AN)
hist(all_AF)
hist(all_DP, breaks = 200)

gts <- chrom@vcf@gt
# snps genotyped
snps_per_ind <- rowSums(!is.na(gts))
hist(snps_per_ind, breaks=96, xlab = "number of individuals", ylab = "snps typed", main = "")

# get one snp per contig, the one with the highest AN
# contig id
head(nes)
# extract fixed part of vcf
fixed <- as.data.frame(chrom@vcf@fix, stringsAsFactors = FALSE)

extract_info <- function(element, vcfR){
    element <- extract.info(vcfR, element = element, as.numeric = TRUE)
}

info_df <- lapply(c("AF", "AN", "DP", "MQ"), extract_info, nes)
info_df <- data.frame(do.call(cbind, info_df), stringsAsFactors = FALSE)
names(info_df) <- c("AF", "AN", "DP", "MQ")
info_all_num <- info_df
# split up information in INFO part into columns
# splitted_info <- str_split(fixed$INFO, ";")
# 
# extract_INFO <- function(x) {
#     df <- as.data.frame(str_split(x, "="), stringsAsFactors = FALSE)
#     names(df) <- df[1, ]
#     out <- df[-1, ]
# }
# 
# info_df <- lapply(splitted_info, extract_INFO)
# info_col_filt <- lapply(info_df, function(x) x[c("AF", "AN", "DP", "MQ")])
# 
# info_all <- do.call(rbind, info_col_filt)
# info_all_num <- as.data.frame(apply(info_all, 2, as.numeric))

# create df with INFO in seperate columns
snp_filter_df <- cbind(fixed, info_all_num)
snp_filter_df <- snp_filter_df[-8]
head(snp_filter_df)
str(snp_filter_df)
snp_filter_df$POS <- as.numeric(snp_filter_df$POS)

# filter indels and multiallelic SNPs (where REF or ALT is longer than 1)
no_multi_allel <- !((str_count(snp_filter_df$REF) > 1) | (str_count(snp_filter_df$ALT) > 1))
snps_filtered1 <- snp_filter_df[no_multi_allel, ]

# filter AF 0 or 1
# calculate DP CI 95% and filter DP according to that
# also filter MAF < 0.025
library(dplyr)
CI <- 0.95
CI_DP <- stats::quantile(snps_filtered1$DP, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
CI_AF <- stats::quantile(snps_filtered1$AF, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)

snps_filtered2 <- snps_filtered1 %>% 
    filter(DP > CI_DP[1]) %>%
    filter(DP < CI_DP[2]) %>%
    filter(AF > 0.025) %>%
    filter(AF < 0.975)

# ###filter to get one snp per contig, the snp with the highest AN
head(snps_filtered2)
# filter for AN first, then for DP and then for the most left (min) position
snps_filtered3 <- snps_filtered2 %>% 
                        group_by(CHROM) %>% 
                        filter(AN==max(AN)) %>%
                        filter(DP==max(DP)) %>%
                        filter(POS==min(POS))
plot(snps_filtered3$DP)


whitelist_snps <- snps_filtered3[c("CHROM", "POS")]

write.table(whitelist_snps, file = "data/whitelist_snps.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

## command line filter initial vcf from whitelist
# vcftools --vcf data/final_all_merged.vcf --out data/filtered_snps --positions data/whitelist_snps.txt --recode --recode-INFO-all



### now filter for DP per sample with matrix

# vcf_file <- "data/filtered_snps.recode.vcf"
vcf_file <- "data/filtered_snps.recode.vcf"
# read vcf
nes <- read.vcfR(vcf_file, verbose = FALSE)
dp <- extract.gt(nes, element = "DP", as.numeric = TRUE)

dp_under_20 <- t(dp < 20)
dp_under_10 <- t(dp < 10)
dp_under_30 <- t(dp < 30)
dp_under_5 <- t(dp < 5)
sum(dp_under_10)

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

# create genotypes

gt <- extract.gt(nes, element = "GT")
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



library(inbreedR)
nes <- snp_genotypes
nes[dp_under_5] <- NA
# some checks
typed <- rowSums(!is.na(nes))
hist(typed)
het <- MLH(nes)
plot(het, typed)

snp_geno_dp5 <- as.matrix(snp_genotypes)
log_mat <- !(unname(as.matrix(dp_under_5)))
rowSums(log_mat)
# snp_geno_dp20[!log_mat] <- NA

rowSums(!is.na(snp_genotypes))



#snps_filtered3 <- snps_filtered2 %>% group_by(CHROM) %>% slice(which.max(AN))

gts <- as.data.frame(chrom@vcf@gt)

contig <- sapply(fixed[, 1], function(x) str_split(x, "\\|")[[1]][1])
df_filter <- data.frame("contig" = contig, "AN" = all_AN)

# filter for highest AN
# library(dplyr)

# unique contig names
unique_contigs <- unique(contig)

library(plyr)
ddply(df_filter, "contig", subset, AN == max(AN))

sel <- ave(df_filter$contig, df_filter$AN, FUN = max) == df_filter$contig

rows_to_keep <- sapply(unique_contigs, function(x) which.max(df_filter[df_filter$contig == x, ]$AN))





# extract dp information from genotype slot
dp <- extract.gt(nes, element='DP', as.numeric=TRUE)
# boxplots of sequence depth
par(mar=c(8,4,1,1))
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Depth", las=2)
abline(h=seq(0,1e4, by=100), col="#C0C0C088")

# violin plots
if( require(reshape2) & require(ggplot2) ){
    dpf <- melt(dp, varnames=c('Index', 'Sample'), value.name = 'Depth', na.rm=TRUE)
    dpf <- dpf[ dpf$Depth > 0,]
    p <- ggplot(dpf, aes(x=Sample, y=Depth)) + geom_violin(fill="#C0C0C0", adjust=1.0,
        scale = "count", trim=TRUE)
    p <- p + theme_bw()
    p <- p + theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 60, hjust = 1))
    #  p <- p + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="black")
    p <- p + scale_y_continuous(trans=scales::log2_trans(), breaks=c(1, 10, 100, 800))
    p
} else {
    message("The packages reshape2 and ggplot2 are required for this example but do not appear
          to be installed.  Please use install.packages(c('reshape2', 'ggplot2', 'scales')) if you would
          like to install them.")
}

##### filtering on sequence depth ######
# get confidence interval for sequence depth per sample
sums <- apply(dp, MARGIN=2, quantile, probs=c(0.05, 0.95), na.rm=TRUE)
# substract lower border from the values in dp
dp2 <- sweep(dp, MARGIN=2, FUN = "-", sums[1,])
# where substraction is 0 or smaller, make it NA
dp[dp2 < 0] <- NA
# substract upper border from the values in dp
dp2 <- sweep(dp, MARGIN=2, FUN = "-", sums[2,])
# where value is bigger, make it NA
dp[dp2 > 0] <- NA
# filter for depth of coverage of at least 4
dp[dp < 4] <- NA

## plot everything
par(mar=c(8,4,1,1))
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Depth")
abline(h=seq(0,200, by=20), col="#C0C0C088")

# some numbers
dim(dp)
plot(colMeans(dp, na.rm = TRUE))

# chromR format
chrom <- create.chromR(name = 'supernes', vcf=nes)
fixed <- chrom@vcf@fix
fixed[1:10, 8]

test <- extract.info(chrom, "AN")

mean(chrom@var.info$DP)
plot(chrom)
chromoqc(chrom)


# filtering
chrom2 <- masker(chrom, min_MQ = 20, min_DP = 300)
plot(chrom2)
# extract genotypes
gt <- extract.gt(nes)
# transpose and data.frame
gt <- as.data.frame(t(gt), stringsAsFactors = FALSE)
# NA handling
gt[gt == "."] <- NA
# split columns
snp_geno <- do.call(cbind, apply(gt, 2, function(x) colsplit(x, "/", c("a","b"))))
# convert
mouse_snp_genotypes <- inbreedR::convert_raw(snp_geno)
# check data
check_data(mouse_snp_genotypes)