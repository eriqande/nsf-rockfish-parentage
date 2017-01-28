

# this is script that cycles over all the RDS files.  

# it uses a little data frame that has the names and GTSeq run numbers
# I got that file like this:
# 2017-01-27 16:26 /haplot-rds/--% (master) pwd
# /Users/eriq/Documents/git-repos/nsf-rockfish-parentage/extdata/haplot-rds
# 2017-01-27 16:26 /haplot-rds/--% (master) ls -l | awk 'BEGIN {print "gtseq_run", "file"} NR >1 && !/posinfo/ {num = $NF; gsub(/gtseq/, "", num); gsub(/[._].*$/, "", num);  print num, $NF}' > ../../data/rds-file-list.txt 

files <- read.table("data/rds-file-list.txt", stringsAsFactors = FALSE, header = TRUE)
dir <- "extdata/haplot-rds/"
