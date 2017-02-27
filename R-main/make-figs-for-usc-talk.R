

library(ggmap)


# get the meta data
meta <- readRDS("./data/processed/meta-data-tibble.rds") %>%
  filter(PROJECT_NAME == "NSF Rockfish Dispersal Project") %>%
  filter(!is.na(LATITUDE_M) & !is.na(LONGITUDE_M)) %>%
  filter(REPORTED_LIFE_STAGE %in% c("ADULT", "JUVENILE"))




#### Make some sampling maps ####
nums <- meta %>%
  group_by(REPORTED_LIFE_STAGE) %>% 
  tally()

latmean <- mean(meta$LATITUDE_M)
longmean <- mean(meta$LONGITUDE_M)

cb_map <- get_map(location = c(-121.92, 36.56), source = "google", maptype = "satellite", zoom = 11)

ggmap(cb_map) +
  geom_jitter(data = meta, mapping = aes(x = LONGITUDE_M, y = LATITUDE_M), colour = "yellow", alpha = 0.5) +
  facet_wrap(~ REPORTED_LIFE_STAGE) +
  geom_text(data = nums, aes(label = n), x = -122.1, y = 36.7, colour = "yellow")


ggsave("sampling_map1.pdf", width = 8, height = 6)




#### Make some plots of numbers of alleles ####

genos <- readRDS("Rmd/rds_outputs/kelp_genos_used_long.rds")

num_alleles <- genos %>%
  filter(!is.na(Allele)) %>%
  group_by(Locus) %>%
  summarise(num_alleles = n_distinct(Allele))

ggplot(num_alleles, aes(x = num_alleles)) + 
  geom_histogram(fill = "blue") +
  xlab("Number of Alleles (Haplotypes) at a Locus") +
  ylab("Number of Loci")

ggsave("alle-count-histo.pdf", width = 5, height = 3)



#### Plots of read depth ####
read_genos <- readRDS("./extdata/processed/genos-aggregated-and-no-hi-missers.rds")

reads <- read_genos %>% 
  filter(!is.na(allele.balance))

ggplot(reads, aes(x = depth)) + 
  geom_histogram(fill = "violet", alpha = 0.5)

depth_bins <- reads %>% 
  mutate(bins = cut(depth, breaks = c(0,10,20,30,40,50, 75, 100, 250, 500, 1000, 10^6))) %>%
  group_by(bins) %>% 
  tally() %>%
  ungroup() %>%
  mutate(percentage = 100 * n/sum(n),
         cumul_perc = cumsum(percentage))

depth_bins %>%
  mutate(Percentage = sprintf("%.1f", percentage),
         Cumulative = sprintf("%.1f", cumul_perc) ) %>%
  rename(`Read Depth Bin` = bins) %>%
  select(-percentage, -cumul_perc) %>%
  write.table(sep = "  &  ", eol = "\\\\\n", quote = F, row.names = F)



#### Get plots of homozygosities

# i did this in the Rmarkdown document

#### Read depths of mismatching loci

# read in the output from the Rmarkdown
matchers <- readRDS("Rmd/rds_outputs/matchers.rds") %>%
  tbl_df()

# pick out just the ones we want
fish1 <- read_genos %>%
  filter(NMFS_DNA_ID %in% matchers$NMFS_DNA_ID_1) %>%
  select(NMFS_DNA_ID, locus, gene_copy, allele, depth)

fish2 <- read_genos %>%
  filter(NMFS_DNA_ID %in% matchers$NMFS_DNA_ID_2) %>%
  select(NMFS_DNA_ID, locus, gene_copy, allele, depth)


together <- left_join(matchers, fish1, by = c("NMFS_DNA_ID_1" = "NMFS_DNA_ID")) %>%
  left_join(., fish2, by = c("NMFS_DNA_ID_2" = "NMFS_DNA_ID", "locus"))

# now we group by pair and locus and test whether the genotypes are concordant or not
tog_conc <- together %>%
  group_by(NMFS_DNA_ID_1, NMFS_DNA_ID_2, locus) %>%
  mutate(concordant = sort(allele.x) == sort(allele.y))

# check that we have the right results:
tog_conc %>%
  group_by(NMFS_DNA_ID_1, NMFS_DNA_ID_2) %>%
  summarise(num_disc = sum(!concordant, na.rm=TRUE) / 2) %>%
  arrange(desc(num_disc))

# Yep!! that is right.  So, now, let's look at the read-depths there
disco <- tog_conc %>%
  filter(concordant == FALSE) %>%
  group_by(NMFS_DNA_ID_1, NMFS_DNA_ID_2, locus) %>%
  mutate(homoz.x = allele.x[1] == allele.x[2],
         homoz.y = allele.y[1] == allele.y[2])

disco %>%
  select(NMFS_DNA_ID_1, NMFS_DNA_ID_2, locus, allele.x, allele.y, depth.x, depth.y)

# too much work for now, gonna do a more complete analysis later.  