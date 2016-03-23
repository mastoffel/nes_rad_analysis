# creates a whitelist from files with individual names

ind_names <- read.table("ind_names.txt")
ind_names <- as.character(unlist(ind_names))

library(stringr)
all_ids <- ind_names[str_detect(ind_names, "fq")]
all_ids <- all_ids[!(str_detect(all_ids, "superparent"))]
all_ids <- all_ids[str_detect(all_ids, ".2.fq")]

#all_ids <- gsub( "\\..*$", "", ind_names )
write(all_ids, file = "whitelist")
