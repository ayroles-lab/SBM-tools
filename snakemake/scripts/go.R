my_logfile = snakemake@log[["log"]]
#my_logfile = "SBM/snakemake/logs/GO/fdr-1e-3/hs.log"
snakemake@source("logger.R")
#source(here::here("SBM/snakemake/scripts/logger.R"))
log4r_info("Starting.")
print = log4r_info

print("Loading packages and functions")

org.db = snakemake@params[["orgdb"]]
#org.db = "org.Dm.eg.db"
#pak::pkg_install(org.db)
library(org.db, character.only=TRUE)

snakemake@source("go_functions.R")
#source(here::here("SBM/snakemake/scripts/go_functions.R"))

block_dir = snakemake@params["blockDir"]
#block_dir = here::here("SBM/snakemake/cache/blockSummary/fdr-1e-3/hs")
print(paste("Data folder:", block_dir))

print("Running GO analysis")
go_set = makeEnrichment(block_path = block_dir)
#go_set = makeEnrichment(block_path = "cache/blockSummary/fdr-1e-3/LUNG")

print("Writting GO object")
saveRDS(go_set, file = snakemake@output[["GO"]])

print("Writting block summary table")
write.csv(go_set$summary, file = snakemake@output[["blockSummary"]])

go_set_table = ldply(go_set$CP, function(x) x@result) %>%
  filter(p.adjust < 0.05, Count >= 4) %>% rename(Name = .id)

print("Writting GO table")
write.csv(go_set_table, file = snakemake@output[["GOcsv"]])

print("Done!")
