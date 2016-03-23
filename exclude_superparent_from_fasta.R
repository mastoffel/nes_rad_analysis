# exclude superparents from fasta files
args<-commandArgs(TRUE)
# read in fasta
fa_file <- readLines(args[1])
# search lines with superparent
all_sp <- grep("superparent", fa_file)
# delete sequences
fa_file <- fa_file[-c(all_sp, all_sp + 1)]
# overwrite old file
write(fa_file, file = args)
