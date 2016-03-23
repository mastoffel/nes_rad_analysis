# raw reads before denovo map


all_reads <- read.table("rad_reads.txt")
all_reads <- all_reads[-nrow(all_reads), ]
libs <- factor(rep(1:6, each = 32))
all_reads$libs <- libs
all_reads$reads <- all_reads$V1 / 4
library(ggplot2)

p <- ggplot(all_reads, aes(x = 1:nrow(all_reads), y = reads))
p + geom_point(aes(colour = libs), size = 4)

