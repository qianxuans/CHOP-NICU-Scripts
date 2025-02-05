library(ggtree)
library(treeio)
library(ggplot2)
#other
library(here)
# Input ----
cluster_summary <- readRDS("cluster_summary.RDS")
tree_beast_dir <- "nexus/"
tree_beast_list <- list.files(path = tree_beast_dir,
                              pattern = "*best.trees",
                              all.files = TRUE,
                              full.names = TRUE,
                              recursive = TRUE)

for(tree_path in tree_beast_list){
  
  cluster_id <- gsub("\\_best.trees|cluster","",basename(tree_path))
  cluster_entry <- which(sapply(cluster_summary[["clusters"]],function(cluster) cluster$cluster == cluster_id))
  cluster_mrsd <- paste0("20",cluster_summary[["clusters"]][[cluster_entry]]$last_seen)
  
  tree_beast <- read.beast(tree_path)
  
  tree_plot <- ggtree(tree_beast,
                      mrsd = as.Date(cluster_mrsd)) + 
    theme_tree2() + 
    geom_tiplab(size = 2.5) + 
    geom_range("CAheight_0.95_HPD",
               color = "#6FB3E1",
               center = "auto")+
    scale_x_continuous(expand = expansion(mult = c(0, 0.5)))
  
  pdf(file=paste0("BEAST2_tree_Cluster",
                  cluster_id,
                  ".pdf"),
      width=10,
      height=7.5)
  print(tree_plot)
  dev.off()
}

