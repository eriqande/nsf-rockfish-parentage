library(tidyverse)
library(readxl)
library(stringr)


source("R/rockfish-funcs.R")


#### Call genos from the microhaplot RDS files ####

# this next step uses a little data frame that has the names and GTSeq run numbers
# I got that file like this:
# 2017-01-27 16:26 /haplot-rds/--% (master) pwd
# /Users/eriq/Documents/git-repos/nsf-rockfish-parentage/extdata/haplot-rds
# 2017-01-27 16:26 /haplot-rds/--% (master) ls -l | awk 'BEGIN {print "gtseq_run", "file"} NR >1 && !/posinfo/ {num = $NF; gsub(/gtseq/, "", num); gsub(/[._].*$/, "", num);  print num, $NF}' > ../../data/rds-file-list.txt 

# get the names of the files
fdf <- read.table("data/rds-file-list.txt", stringsAsFactors = FALSE, header = TRUE) %>%
  tbl_df()
dir <- "extdata/haplot-rds"


# cycle over them, read them and add the gtseq_run column on each.
# at the end, bind them together.
genos_long <- lapply(1:nrow(fdf), function(i) {
  message("Working on ", fdf$file[i])
  call_genos_from_haplotRDS(path = file.path(dir, fdf$file[i])) %>%
    mutate(gtseq_run = fdf$gtseq_run[i]) %>%
    select(gtseq_run, everything())
}) %>%
  bind_rows()


# that gives us about 2.8 million records, and we go ahead and save
# it in extdata/processed, with xz compression in 5.5 Mb
saveRDS(genos_long, file = "extdata/processed/called_genos.rds", compress = "xz")


#### Read the sample sheets from all the excel files ####

# Same drill here.  First we make a file that holds a data frame of file names:
# 2017-01-28 15:36 /sample_sheets/--% (master) pwd
# /Users/eriq/Documents/git-repos/nsf-rockfish-parentage/data/sample_sheets
# 2017-01-28 15:36 /sample_sheets/--% (master) ls -l | awk 'BEGIN {print "gtseq_run", "file"} NR > 1 {num = $NF; gsub(/GTseq/, "", num); gsub(/_.*$/, "", num);  print num, $NF}' > ../../data/sample-sheet-file-list.txt

fdf <- read.table("data/sample-sheet-file-list.txt", stringsAsFactors = FALSE, header = TRUE) %>%
  tbl_df()
dir <- "data/sample_sheets"

sample_sheets <- lapply(1:nrow(fdf), function(i) {
  message("Working on ", fdf$file[i])
  read_excel(path = file.path(dir, fdf$file[i]), sheet = "sample_sheet", skip = 18) %>%
    tidyr::separate(Sample_Plate, into = c("NMFS_DNA_ID", "ssBOX_ID", "ssBOX_POSITION")) %>%
    mutate(id = str_replace(Sample_ID, "satrovirens_0*", "s")) %>%
    mutate(gtseq_run = fdf$gtseq_run[i]) %>%
    select(gtseq_run, id, everything())
}) %>%
  bind_rows()

saveRDS(sample_sheets, "data/sample-sheet-tibble.rds", compress = "xz")



#### And finally, let's get the meta-data read in and cleaned up (if it needs it) ####
meta <- read_excel("data/nsf-rockfish-metadata.xlsx", sheet = 1)

# when we read that in, we lose the "None."s in the LEFTOVER_SAMPLE fields.  That 
# is OK for now.  
names(meta) <- str_replace(names(meta), "^Marine::", "marine_")  # get the colons out of the column names

saveRDS(meta, "data/meta-data-tibble.rds", compress = "xz")




