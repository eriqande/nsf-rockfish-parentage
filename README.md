nsf-rockfish-parentage
================
31 January, 2017

-   [Data Files](#data-files)
    -   [Other data files](#other-data-files)
-   [Early maneuvers. Data Checking, etc.](#early-maneuvers.-data-checking-etc.)

<!-- README.md is generated from README.Rmd. Please edit that file -->
This is an Rstudio project housing our analyses of parentage (and sibling inference) of rockfish using our GTseq data off the MiSeq.

The full genotype data don't live in this repo, since that would be too much for GitHub, but, the code for the early data processing steps is set up as if you have them in a directory called `extdata` in this repository. The `data` directory here holds the meta data information for the individuals, and the plate maps that we use to associate the IDs from the genotyping data with all the metadata.

Data Files
----------

Check out `./R-main/01-compile-and-tidy-data.R` to see what we went through to tidy data up into a reasonably good long format.

Those steps produce these key files. The first two are part of this repository on GitHub. The third is not, at the present time.

1.  `./data/meta-data-tibble.rds` an xz-compressed tibble of all the meta data from the SWFSC tissue repository. Open it in R with `readRDS()`.
2.  `./data/sample-sheet-tibble.rds` an xz-compressed tibble of all the sample sheet information. This lets us associate `gtseq_run` number and genotyping `id` with a `NMFS_DNA_ID` so we can associate meta-data with the genotyped individuals.
3.  `./extdata/processed/called_genos.rds` an xz-compressed tibble of all the called genotypes in a long format. Individuals are identified by the `gtseq_run` and `id`, which can be used to get back to their `NMFS_DNA_ID`.
4.  `./extdata/processed/called_genos_na_explicit.rds` an xz-compressed tibble of all the genotypes in long format. This one has explicitly put NAs in for every individual listed in the sample sheet (in the previous file, individuals with no genotypes scored at all don't appear!). It also includes the `NMFS_DNA_ID`. This is the one to work from.

All subsequent analyses start from these files.

### Other data files

There are a few more intermediates that I am generating as I work through this. Two very important ones are:

1.  `./data/processed/juvie-gsi-sim-assignments.rds` a data frame of GSI-assignments of juvenile fish to a baseline that included adults of a few different species.
2.  `./extdata/processed/genos-aggregated-and-no-hi-missers.rds` a data frame of all the genotypes for individuals that are typed at &gt;= 85 loci. This has the NMFS\_DNA\_ID in it and that is what should be used, as multiply-genotyped individuals have had the highest read-depth result put in there. Hence NMFS\_DNA\_IDs might have genotypes taken from more than gtseq\_run.

Early maneuvers. Data Checking, etc.
------------------------------------

1.  First things first, verify who it is we have (and have not) genotyped, etc. This is just a general overview of the data. R-markdown notebook: `./Rmd/general-data-overview.Rmd`
2.  Check a few things for genotyping accuracy in `./Rmd/genotyping-accuracy-assessments.Rmd`
3.  Aggregate individual genotypes, toss hi-missers, and do GSI to identify species of juvies in `/Rmd/sorting-to-species.Rmd`.
