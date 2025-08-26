# Transmission Trees Inference using TransPhylo
library(TransPhylo)
library(ape)
library(treeio)
library(lubridate)
library(coda)
library(here)
# Helper Function ----
# Convert the beast tree to newick format
convert2transmission <- function(beast_tree_path,
                                 cluster_date_dir,
                                 surveillance_end = 2024.495,
                                 mcmc_iteration = 1000000,
                                 analysis_seed = 1205){
  
  set.seed(analysis_seed)
  # cluster ID
  cluster_id <- gsub("_best.trees","",basename(beast_tree_path))
  
  cluster_name <- gsub("cluster","Cluster",cluster_id)
  
  # read the beast tree
  beast_tree <- read.nexus(beast_tree_path)
  # convert to newick format tree with the unit of branch length being time (years)
  # read the date for inferring transmission tree 
  cluster_date <- setNames(read.table(file.path(cluster_date_dir,paste0(cluster_id,"_date.txt")),
                                      header = FALSE),
                           c("genome","collection_date"))
  
  cluster_last_date <- decimal_date(max(as.Date(cluster_date$collection_date)))
  # latest date used in inferring transmission tree
  
  cluster_ptree <- ptreeFromPhylo(beast_tree,
                                  dateLastSample=cluster_last_date)
  
  
  cluster_res <- inferTTree(cluster_ptree,
                            mcmcIterations=mcmc_iteration,
                            # shape and scale calculated using NICU cluster data
                            w.shape=0.588,
                            w.scale=0.461,
                            # end of surveillance date in this study is when observation stopped
                            dateT=surveillance_end,
                            updateOff.p = TRUE)
  
  # check MCMC coverage
  cluster_mcmc=convertToCoda(cluster_res)
  
  cluster_mcmc_EffectiveSize <- as.data.frame(t(effectiveSize(cluster_mcmc)))
  
  # Discard the first 10% results as burinin
  
  cluster_med=medTTree(cluster_res,burnin = 0.1)
  
  cluster_index_case <- min(cluster_med$ctree[,1])
  
  cluster_index_transmission <- sort(cluster_med$ctree[,1])[2]
  
  cluster_ttree <- extractTTree(cluster_med)
  
  # Visualize the transmission tree
  pdf(file = file.path("./output/transmission_tree/",paste0(cluster_name,"_transmission_tree.pdf")),width = 15,height = 10)
  plot(output_RDS$transmission_tree)
  title(main = paste0(cluster_name, " Transmission Tree"))
  dev.off()
  
  pdf(file = file.path("./output/colored_phylogeny/",paste0(cluster_name,"_colored_phylogenetic_tree.pdf")),width = 15,height = 10)
  plot(output_RDS$colored_phylo_tree)
  title(main = paste0(cluster_name, " Colored Phylogeny"))
  dev.off()
  
  cluster_res_post_burnin <- cluster_res[(0.1*length(cluster_res)+1):length(cluster_res)]
  
  # Visualize the unsampled cases probabilities
  pdf(file = file.path("./output/incidentcases/",paste0(cluster_name,"_IncidentCases.pdf")),width = 15,height = 10)
  cluster_ic <- getIncidentCases(cluster_res_post_burnin,show.plot = TRUE)
  title(main = paste0(cluster_name, " Sampled and Unsampled Cases"))
  dev.off()
  
  cluster_index_case_range <- sapply(cluster_res_post_burnin, function(res_entry) min(res_entry$ctree$ctree[,1]))
  
  cluster_index_case_ci <- quantile(as.numeric(cluster_index_case_range),c(0.025,0.975))
  
  return(
    list(
      sum_df = data.frame(
        cluster = cluster_id,
        transm_index_mean = cluster_index_case,
        transm_index_lower = as.numeric(cluster_index_case_ci[1]),
        transm_index_upper = as.numeric(cluster_index_case_ci[2])
      ),
      mcmc_check = cluster_mcmc_EffectiveSize,
      transmission_tree = cluster_ttree,
      colored_phylo_tree = cluster_med,
      raw_res = cluster_res,
      raw_res_post_burnin = cluster_res_post_burnin,
      index_case_range = cluster_index_case_range
    )
  )
  }

# Perform analysis ----
## Input ----

setwd(dir = file.path(here(),"output"))
beast_tree_dir <- "./input/beast_tree/"
beast_tree_list <- list.files(path = beast_tree_dir,
                              pattern = "cluster*",
                              recursive = FALSE,
                              all.files = TRUE,
                              full.names = TRUE)

## Execute function ----
for(tree_path in beast_tree_list){
  
  cluster_id <- gsub("_best.trees","",basename(tree_path))
  
  cluster_sum <- convert2transmission(beast_tree_path = tree_path,
                       cluster_date_dir = "./input/cluster_date/",
                       surveillance_end = 2024.495,
                       mcmc_iteration = 1000000,
                       analysis_seed = 1205)
  
  
  saveRDS(cluster_sum,
          paste0(cluster_id,"_TransPhylo.RDS"))
}
