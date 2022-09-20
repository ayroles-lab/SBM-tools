library(clusterProfiler)
library(enrichplot)
library(readr)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(cowplot)

tissue = "HEART"

go_folder = paste0(here::here(), "/data/output/SBM/gtex/GO/fdr-1e-3/")
go_plots = paste0(here::here(), "/data/output/SBM/gtex/GO_plots/fdr-1e-3/")
dir.create(go_plots, recursive = TRUE, showWarnings = FALSE)
tissues = dir(go_folder)
tissue = tissues[1]

enGoSummary_list = llply(tissues, function(tissue) read_csv(file.path(go_folder, tissue, "block_summary.csv")) %>%
  mutate(tissue = tissue))
  names(enGoSummary_list) = tissues
enGoSummary_all = ldply(enGoSummary_list)

p = ggplot(enGoSummary_all, aes(Nested_Level, Assortativity, 
           group = interaction(tissue, Nested_Level), 
           fill = tissue, 
           color = tissue)) + 
           geom_boxplot() +geom_line() 
save_plot("test.png", p)

ldply(tissues, function(tissue) ldply(1:5, function(level, tissue) {
table_en = table(enGoSummary_list[[tissue]]$n_enrich[enGoSummary_list[[tissue]]$Nested_Level==level]!=0)
data.frame(tissue = tissue,
           n_blocks = sum(table_en), 
          "n_enriched" = table_en["TRUE"], 
          "proportion" = table_en["TRUE"]/sum(table_en),
           Level = level, row.names = "")
}, tissue)) %>% arrange(Level)

library(doMC)
registerDoMC(10)
llply(tissues, function(tissue){
    enGo = readRDS(file.path(go_folder, tissue, "GO.rds"))
    en_s = enGo$summary %>%
        filter(n_enrich > 1, N_genes > 4)
    enriched_blocks = enGo$CP[en_s$Name]
    plot_path = paste0(go_plots, tissue)
    dir.create(paste0(plot_path), recursive = TRUE, showWarnings = FALSE)
    for(block in en_s$Name){
        p1 = dotplot(enriched_blocks[[block]], showCategory=30) + ggtitle(block)
        #save_plot("test.png", p1)
        level = en_s[en_s$Name==block,]$Nested_Level
        save_plot(file.path(plot_path, paste0("Level_", level, "_", block, ".png")), p1, base_height = 8)
    }
}, .parallel = TRUE)


    # pw_upper <- pairwise_termsim(enriched_blocks[[block]])
    # p1 <- emapplot(pw_upper, showCategory = 15) + 
    #             theme(legend.position = "none") +
    #             theme(plot.title = element_text(size=28)) +
    #             ggtitle(block)
    # p1$data$color = mean(p1$data$color)
