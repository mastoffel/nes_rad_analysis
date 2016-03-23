# construct denovomap

all_reads <- read.table("rad_reads.txt")
all_reads <- all_reads[-nrow(all_reads), ]

fileConn<-file("output.txt")

for (i in seq(from = 1, to = nrow(all_reads), by = 2)) {
    if (i == 1) {
    text <- paste0("-r ", as.character(all_reads$V2)[i], " \\") # in write.table will be one backslash
    } else {
        text <- paste0(text, "\n", paste0("-r ", as.character(all_reads$V2)[i], " \\"))
    }
}

write.table(text, file = "test.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
