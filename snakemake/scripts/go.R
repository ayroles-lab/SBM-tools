my_logfile = snakemake@log[["log"]]
snakemake@source("logger.R")
log4r_info("Starting.")
print = log4r_info

print("Loading packages and functions")
snakemake@source("go_functions_gtex.R")

print(paste("Data folder:", snakemake@params["blockDir"]))

print("Running GO analysis")
go_set = makeEnrichment(block_path = snakemake@params["blockDir"])
#go_set = makeEnrichment(block_path = "../data/output/SBM/gtex/blockSummary/fdr-1e-3/LUNG")

print("Writting GO object")
saveRDS(go_set, file = snakemake@output[["GO"]])

print("Writting block summary table")
write.csv(go_set$summary, file = snakemake@output[["blockSummary"]])

go_set_table = ldply(go_set$CP, function(x) x@result) %>%
  filter(p.adjust < 0.05, Count >= 4) %>% rename(Name = .id)

print("Writting GO table")
write.csv(go_set_table, file = snakemake@output[["GOcsv"]])

print("Done!")
