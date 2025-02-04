# Libraries ----
#use parallel package to use multiple cores
library(parallel)
#use compiler package to compile the function and improve the performance
library(compiler)
#packages to analyze the trees
library(TreeTools)
library(ape)
#clustering tool
library(cluster)
#others 
library(dplyr)
library(purrr)
library(here)
# Input ----
output_path <- "./"
everything_tree_path <- "./input/NICU_1446_ufboot_bnni_alrt.treefile"
study_snp_matrix_path <- "./input/NICU_NICU_sum_snp_unique.RDS"
study_global_snp_matrix_path <- "./input/NICU_NCBI_sum_snp_unique.RDS"
group_tree_path <- "./input/group_trees/"
setwd(dir = output_path)

# Hierarchical Clustering ----
#A CC-free way to group the study genomes in the everything-tree using the pairwise distances
hierarchical_clustering <- function(everything_tree_path,
                                    study_snp_matrix_path){
  # read the tree
  everything_tree <- ape::read.tree(everything_tree_path)
  # mid-point rooting of the tree
  everything_tree <- phytools::midpoint_root(everything_tree)
  # pre-order the tree
  everything_tree <- Preorder(everything_tree)
  # summarize the tree by iterating through nodes in the tree
  all_nodes <- (length(everything_tree$tip.label) + 1):(length(everything_tree$tip.label) + everything_tree$Nnode)
  
  everything_tree_sum <- lapply(seq_along(all_nodes),
                                function(node_entry){
                                  sub_everything_tree <- Subtree(everything_tree,all_nodes[node_entry])
                                  
                                  if(everything_tree$node.label[node_entry] == "Root"){
                                    
                                    return(list(
                                      node = all_nodes[node_entry],
                                      SHaLRT_support = "Root",
                                      bootstrap_support = "Root",
                                      genomes = sub_everything_tree$tip.label
                                      ))
                                    
                                  }else{
                                    
                                    return(list(
                                      node = all_nodes[node_entry],
                                      SHaLRT_support = as.numeric(strsplit(everything_tree$node.label[node_entry],
                                                                split = "/")[[1]][1]),
                                      bootstrap_support = as.numeric(strsplit(everything_tree$node.label[node_entry],
                                                                   split = "/")[[1]][2]),
                                      genomes = sub_everything_tree$tip.label
                                    ))
                                  }
                                })
  
  #Matrix of pairwise distances from tree
  distance_matrix <-  cophenetic.phylo(everything_tree)
  #the various methods can become options in the future
  hclust_result <- hclust(as.dist(distance_matrix),
                          method = "complete")
  #hierarchical clustering
  #We test the range from 2 to (number of genomes in the everything - 1) clustering patterns 
  hclust_result_silhouette_scores <- as.numeric(unlist(sapply((2:(length(everything_tree$tip.label)-1)),
                                                              function(clustering_pattern_entry){
                                                                clusters <- cutree(hclust_result, k = clustering_pattern_entry)
                                                                return(mean(silhouette(clusters,distance_matrix)[, "sil_width"]))
                                                              })))
  
  final_group <- cutree(hclust_result,
                        k = (which.max(hclust_result_silhouette_scores) + 1))

  #Read the SNP matrix of study genomes comparing to study genomes to perform the QC
  study_snp_matrix <- readRDS(study_snp_matrix_path)
  #Look at each group to check the pair-wise SNP distance
  #If the genome in the group has no link with other genomes within 1000 gsnps (can be changed)
  #The genome is a singletons
  
  hc_groups_qc <- lapply(sort(unique(final_group)),
                         function(hc_group){
                           
                           group_genomes <- names(final_group)[final_group == hc_group]
                           
                           if(length(group_genomes) > 1){
                             
                             # Find the SHaLRT/bootstrap support of this HC group
                             sum_entry <- which(sapply(everything_tree_sum,
                                                       function(summary){identical(sort(summary$genomes),
                                                                                   sort(group_genomes))}))
                             
                             group_SHaLRT_support <- everything_tree_sum[[sum_entry]]$SHaLRT_support
                             group_bootstrap_support <- everything_tree_sum[[sum_entry]]$bootstrap_support
                             
                             # If all genomes over 1000gsnp limit? 
                             group_study_snp_matrix <- study_snp_matrix[(study_snp_matrix$subject %in% group_genomes) &
                                                                          (study_snp_matrix$query %in% group_genomes),]
                             
                             if(min(group_study_snp_matrix$gsnp) >= 1000){
                               
                               group_all_overlimit <- TRUE
                               genomes_overlimit <- group_genomes
                               
                             }else{
                               
                               group_all_overlimit <- FALSE
                               
                               # Iterate through group_genomes to find which genomes has >= 1000 gsnp
                               genomes_overlimit <- as.character(unlist(sapply(group_genomes,
                                                                  function(genome){
                                                                    group_study_snp_matrix_genome <- group_study_snp_matrix[group_study_snp_matrix$subject == genome |
                                                                                                                              group_study_snp_matrix$query == genome,]
                                                                    if(min(group_study_snp_matrix_genome$gsnp) >= 1000){
                                                                      return(genome)
                                                                    }
                                                                  })))
                             }
                             
                             
                             return(list(
                               hc_group = hc_group,
                               genomes = group_genomes,
                               all_overlimit = group_all_overlimit,
                               genomes_overlimit = genomes_overlimit,
                               SHaLRT_support = group_SHaLRT_support,
                               bootstrap_support = group_bootstrap_support
                             ))
                             
                           }else{
                             # If there is only 1 group genome: singleton
                             return(list(
                               hc_group = hc_group,
                               genomes = group_genomes,
                               all_overlimit = FALSE,
                               genomes_overlimit = NULL,
                               SHaLRT_support = "N/A",
                               bootstrap_support = "N/A"
                             ))
                           }
                         })
  
  #return hc_groups_qc
  return(hc_groups_qc)
}

hierarchical_clustering <- cmpfun(hierarchical_clustering)

hierarchical_clustering_groups <- hierarchical_clustering(everything_tree_path,
                                                          study_snp_matrix_path) 

# Intermediate Steps ----
## Assembly-Scan to find N50 ----
## WhatsGNU introduces the closest global genomes to the study genomes ----
## Snippy uses highest N50 within the group as the reference, and generates the msa of the group ----
## IQtree uses the msa from snippy to build the trees ----
# Input for determining strains ----
# Output from iqtree
group_trees <- list.files(path = group_tree_path,
                          pattern = "\\.contree.treefile$",
                          all.files = TRUE,
                          full.names = FALSE,
                          recursive = FALSE)

get_determine_strains_input <- function(hierarchical_clustering_groups,
                                        group_tree_path,
                                        group_tree,
                                        study_snp_matrix_path,
                                        study_global_snp_matrix_path){
  ## Group details ----
  group_id <- as.integer(gsub("Group|\\.contree.treefile","",group_tree))
  qc_entry <- which(unlist(sapply(hierarchical_clustering_groups,
                                  function(group){
                                    group$hc_group
                                  })) == group_id)
  group_genomes <- hierarchical_clustering_groups[[qc_entry]]$genomes
  # Exclude genomes over gsnp limit 
  group_genomes <- group_genomes[!(group_genomes %in% hierarchical_clustering_groups[[qc_entry]]$genomes_overlimit)]
  ##  SNP-distance matrix ----
  study_snp_matrix <- readRDS(study_snp_matrix_path)
  study_snp_matrix <- study_snp_matrix[study_snp_matrix$subject %in% group_genomes &
                                         study_snp_matrix$query %in% group_genomes,]
  
  ## Summarize group tree ----
  group_tree_newick <- read.tree(file.path(group_tree_path,group_tree))
  group_tree_newick <- Preorder(group_tree_newick)
  group_tree_newick_all_nodes <- (length(group_tree_newick$tip.label) + 1):(length(group_tree_newick$tip.label) + group_tree_newick$Nnode)
  group_tree_newick_sum <- lapply(seq_along(group_tree_newick_all_nodes),
                                  function(node_entry){
                                    sub_group_tree_newick <- Subtree(group_tree_newick,
                                                                     group_tree_newick_all_nodes[node_entry])
                                    return(list(
                                      node = group_tree_newick_all_nodes[node_entry],
                                      bootstrap_support = as.numeric(group_tree_newick$node.label[node_entry]),
                                      genomes = sub_group_tree_newick$tip.label
                                    ))
                                  })
  
  ## Test gsnp from min to 500 ----
  lapply((round(min(study_snp_matrix$gsnp)+0.1):500),
         function(gsnp_cutoff){
           ### SNP only ----
           #Find all single links among the study genomes at this gsnp cutoff
           study_snp_matrix_cutoff <- study_snp_matrix[study_snp_matrix$gsnp <= gsnp_cutoff,]
           # Within those links, the genomes in study_snp_matrix_cutoff would be linked to at least one another genome and considered part of a strain
           # Which means that those genomes will not be singletons
           strain_genomes <- unique(c(study_snp_matrix_cutoff$subject,
                                      study_snp_matrix_cutoff$query))
           
           # Find singleton genomes
           singleton_genomes <- c(setdiff(group_genomes,strain_genomes),
                                  hierarchical_clustering_groups[[qc_entry]]$genomes_overlimit)
           
           # Summrize the single-linkage 
           redundant_strain_list <- lapply(seq_along(strain_genomes),
                                           function(genome_entry){
                                             genome <- strain_genomes[genome_entry]
                                             # Find the links associated with this genome
                                             # Every genomes involved in the links associated with this genome will be considered from the same strain
                                             genome_link <- study_snp_matrix_cutoff[study_snp_matrix_cutoff$subject == genome | 
                                                                                      study_snp_matrix_cutoff$query == genome,]
                                             strain_linked_genomes <- unique(c(genome_link$subject,
                                                                               genome_link$query))
                                             return(
                                               list(id = genome_entry,
                                                    genomes = strain_linked_genomes)
                                             )
                                           })
           
           
           
           # Perform single-linkage clustering with redundant_strain_list
           # Strain_id_index is used to track the strain ID, starting with 1
           
           strain_id_index <- 1
           
           # Use both unique_strain_df to perform single-linkage clustering
           
           unique_strain_df <- data.frame(genome = strain_genomes,
                                          category = "clone",
                                          strain_id = NA)
           
           invisible(sapply(strain_genomes,
                            function(genome){
                              # Find all redundant strains involving the genome in redundant_strain_list
                              redundant_entry <- which(sapply(redundant_strain_list,function(redundant_list) genome %in% redundant_list$genomes))
                              # the genomes in redundant strains with redundant_entry are now considered the same strain 
                              unique_strain_genomes <- unique(as.character(unlist(sapply(redundant_entry,
                                                                                         function(entry) redundant_strain_list[[entry]]$genomes))))
                              # Check if any of unique_strain_genomes is already labeled with strain_id
                              # if all are labeled, no action needed 
                              not_labeled_unique_strain_genomes <- unique_strain_df$genome[unique_strain_df$genome %in% unique_strain_genomes & 
                                                                                             is.na(unique_strain_df$strain_id)]
                              
                              if(length(not_labeled_unique_strain_genomes) > 0){
                                if(identical(sort(unique_strain_genomes),
                                             sort(not_labeled_unique_strain_genomes))){
                                  # if none are labeled, label all genomes in unique_strain_genomes with strain_id_index, and add 1 to strain_id_index
                                  unique_strain_df$strain_id[unique_strain_df$genome %in% unique_strain_genomes] <<- strain_id_index
                                  strain_id_index <<- strain_id_index + 1
                                }else{
                                  # If any one of the genome is already labelled a with strain_id
                                  # the rest of the unlabeled genomes will be labeled the same strain_id
                                  # If there are multiple strain_id, use the smallest id
                                  smallest_strain_id <- min(na.omit(unique_strain_df$strain_id[unique_strain_df$genome %in% unique_strain_genomes]))
                                  unique_strain_df$strain_id[unique_strain_df$genome %in% unique_strain_genomes] <<- smallest_strain_id
                                }
                              }
                            }))
           
           if(length(singleton_genomes) > 0){
             # Add singletons to unique_strain_df if there is any
             unique_strain_df <- rbind(unique_strain_df, 
                                       data.frame(genome = singleton_genomes, 
                                                  category = "singleton",
                                                  strain_id = seq_along(singleton_genomes) + max(unique_strain_df$strain_id)))
           }
           
           # Count before-correction singleton/clones
           before_correction_clones <- length(unique(unique_strain_df$strain_id[unique_strain_df$category == "clone"]))
           before_correction_singletons <- length(unique(unique_strain_df$strain_id[unique_strain_df$category == "singleton"]))
           
           ### Phylogenetic Tree ----
           
           # Only process clones (strains that have no less than 2 genomes, not singleton strains)
           # Find the strain number of the strain with no less than 2 genomes (clones)
           # Add bootstrap support for the strains
           unique_strain_df <- unique_strain_df %>% mutate(correction = NA,
                                                           bootstrap_support = NA)
           
           unique_strain_df$bootstrap_support[unique_strain_df$category == "singleton"] <- 0
           
           unique_strain_tree <- lapply(sort(unique(unique_strain_df$strain_id[unique_strain_df$category == "clone"])),
                                        function(strain_id){
                                          # strain genomes determined only with gsnp
                                          strain_snp_genomes <- unique_strain_df$genome[unique_strain_df$strain_id == strain_id]
                                          #iterate group_tree_newick_sum to find the clades covering all strain_snp_genomes
                                          strain_clade <- group_tree_newick_sum[sapply(group_tree_newick_sum, function(clade) 
                                            all(strain_snp_genomes %in% clade$genomes))] %>%
                                            #the best final clade is the one with least number of genomes 
                                            {.[which.min(sapply(., function(clade) length(clade$genomes)))]}
                                          
                                          #return result to the unique_strain_tree
                                          list(
                                            strain_id = strain_id,
                                            snp_only_genomes = sort(strain_snp_genomes),
                                            snp_tree_genomes = sort(strain_clade[[1]]$genomes),
                                            bootstrap_support = strain_clade[[1]]$bootstrap_support,
                                            node = strain_clade[[1]]$node
                                          )
                                        })
           
           ### Corrected for phylogenetic structure -----
           
           correction_strain_tree <- unique_strain_tree[sapply(unique_strain_tree, function(strain) !identical(strain$snp_only_genomes, strain$snp_tree_genomes))]
           # Count discrepancy
           discrepancy <- sum(unlist(sapply(correction_strain_tree,
                                            function(strain){
                                              length(setdiff(strain$snp_tree_genomes[!grepl("GCA_|GCF_",strain$snp_tree_genomes)],
                                                             strain$snp_only_genomes))
                                            })))
           # Only perform correction when there is discrepancy
           if(length(correction_strain_tree) > 0){
             # Filter the list to keep only the strains with discrepancy between tree and gsnp
             # And only accept the correction only if the bootstrap support >= 70 (https://doi.org/10.1093/sysbio/42.2.182)
             correction_strain_tree <- correction_strain_tree[sapply(correction_strain_tree, function(strain) ((strain$bootstrap_support >= 70) | (strain$bootstrap_support == "Root")))]
             clones_corrected <- length(correction_strain_tree)
             # The function to perform Cladebreaker to return the list of Cladebreaker-divided strains
             correction_strain_list <- lapply(correction_strain_tree,
                                              function(strain_clade){
                                                # 2 situations are expected:
                                                if(!any(grepl("GCA_|GCF_",strain_clade$snp_tree_genomes))){
                                                  # 1st situation: only discrepancy genomes with no breaker genomes
                                                  # use the minimal ID for all snp_tree_genomes in unique_strain_df
                                                  return(
                                                    list(
                                                      list(original_strain = strain_clade$strain_id,
                                                           corrected_genomes =  strain_clade$snp_tree_genomes,
                                                           category = "clone",
                                                           bootstrap_support = as.integer(strain_clade$bootstrap_support)
                                                      ))
                                                  )
                                                }else{
                                                  # 2nd situation: breaker genomes were introduced (with or without discrepancy genomes)
                                                  
                                                  strain_clade_newick <- Subtree(group_tree_newick,
                                                                                 strain_clade$node)
                                                  
                                                  strain_clade_newick_all_nodes <- (length(strain_clade_newick$tip.label) + 1):(length(strain_clade_newick$tip.label) + strain_clade_newick$Nnode)
                                                  # From strain_clade_newick_all_nodes
                                                  # We need to find the biggest clades in strain_clade_newick without being broken by breaker genomes
                                                  # And the genomes in those biggest pure-study-genomes clades are from the same strains
                                                  tmp_list <- lapply(seq_along(strain_clade_newick_all_nodes),
                                                                     function(node_entry){
                                                                       
                                                                       node <- strain_clade_newick_all_nodes[node_entry]
                                                                       sub_strain_clade_newick <- Subtree(strain_clade_newick,
                                                                                                          node)
                                                                       #Check if this clade has only study genomes
                                                                       if(!any(grepl("GCA_|GCF_",sub_strain_clade_newick$tip.label))){
                                                                         # If there are only study genomes
                                                                         # proceed to check if this is the biggest pure-study-genome clade (There are breaker genomes in the parent clade)
                                                                         # Parent node:
                                                                         parent_node <- strain_clade_newick$edge[which(strain_clade_newick$edge[,2] == node), 1]
                                                                         parent_sub_strain_clade_newick <- Subtree(strain_clade_newick,
                                                                                                                   parent_node)
                                                                         if(any(grepl("GCA_|GCF_",parent_sub_strain_clade_newick$tip.label))){
                                                                           
                                                                           return(
                                                                             list(
                                                                               original_strain = strain_clade$strain_id,
                                                                               corrected_genomes = sub_strain_clade_newick$tip.label,
                                                                               category = "clone",
                                                                               bootstrap_support = as.integer(strain_clade_newick$node.label[node_entry])
                                                                             )
                                                                           )
                                                                         }
                                                                       }
                                                                     }) %>% purrr::compact()
                                                  
                                                  #Add the singletons to the list
                                                  singletons <- setdiff(strain_clade$snp_tree_genomes[!grepl("GCA_|GCF_", strain_clade$snp_tree_genomes)],
                                                                        unlist(sapply(tmp_list, `[[`, "corrected_genomes")))
                                                  
                                                  #Add the singleton list to the tmp_list
                                                  if(length(singletons) > 0){
                                                    tmp_list <- c(tmp_list,lapply(singletons, function(singleton) list(
                                                      original_strain = strain_clade$strain_id,
                                                      corrected_genomes = singleton,
                                                      category = "singleton",
                                                      bootstrap_support = 0
                                                    )))
                                                  }
                                                  return(tmp_list)
                                                  
                                                }
                                              })
             
             
             # Use correction_strain_list to update info in unique_strain_df
             invisible(sapply(correction_strain_list,
                              function(correction_strain){
                                for(i in seq_along(correction_strain)){
                                  if(length(correction_strain) == 1){
                                    unique_strain_df$strain_id[unique_strain_df$genome %in% correction_strain[[i]]$corrected_genomes] <<- correction_strain[[i]]$original_strain
                                  }else{
                                    unique_strain_df$strain_id[unique_strain_df$genome %in% correction_strain[[i]]$corrected_genomes] <<- paste0(correction_strain[[i]]$original_strain,
                                                                                                                                                 "_",
                                                                                                                                                 i)
                                  }
                                  
                                  unique_strain_df$category[unique_strain_df$genome %in% correction_strain[[i]]$corrected_genomes] <<- correction_strain[[i]]$category
                                  unique_strain_df$correction[unique_strain_df$genome %in% correction_strain[[i]]$corrected_genomes] <<- TRUE
                                  unique_strain_df$bootstrap_support[unique_strain_df$genome %in% correction_strain[[i]]$corrected_genomes] <<- correction_strain[[i]]$bootstrap_support
                                }
                              }))
           }else{
             clones_corrected <- 0
           }
           
           
           unique_strain_df$correction[is.na(unique_strain_df$correction)] <- FALSE
           
           # Fill in the bootstrap support
           invisible(sapply(unique(unique_strain_df$strain_id[is.na(unique_strain_df$bootstrap_support)]),
                            function(id){
                              unique_strain_df$bootstrap_support[unique_strain_df$strain_id == id] <<- unique_strain_tree[[as.integer(id)]]$bootstrap_support 
                            }))
           
           # Count before-correction singleton/clones
           after_correction_clones <- length(unique(unique_strain_df$strain_id[unique_strain_df$category == "clone"]))
           after_correction_singletons <- length(unique(unique_strain_df$strain_id[unique_strain_df$category == "singleton"]))
           
           ### Return result ----
           # Sort the output using unique_strain_df
           strain_composition <- lapply(sort(unique(unique_strain_df$strain_id)),
                                        function(id){
                                          list(
                                            strain_id = id,
                                            category = unique(unique_strain_df$category[unique_strain_df$strain_id == id]),
                                            genome = unique_strain_df$genome[unique_strain_df$strain_id == id],
                                            correction = unique(unique_strain_df$correction[unique_strain_df$strain_id == id]),
                                            bootstrap_support = unique(unique_strain_df$bootstrap_support[unique_strain_df$strain_id == id])
                                          )
                                        })
           
           # Calculate the mean/median bootstrap support for strains
           
           
           
           median_strain_bootstrap_support <- unlist(median(sapply(strain_composition[sapply(strain_composition, function(x) x$category == "clone")], `[[`, "bootstrap_support")))
           
           median_strain_bootstrap_support <- ifelse(is.null(median_strain_bootstrap_support), 
                                                     0,
                                                     ifelse(is.na(as.numeric(median_strain_bootstrap_support)),
                                                            0,
                                                            as.numeric(median_strain_bootstrap_support)))
           
           mean_strain_bootstrap_support <- mean(sapply(strain_composition[sapply(strain_composition, function(x) x$category == "clone")], `[[`, "bootstrap_support"))
           
           mean_strain_bootstrap_support <- ifelse(is.null(mean_strain_bootstrap_support), 
                                                   0,
                                                   ifelse(is.na(as.numeric(mean_strain_bootstrap_support)),
                                                          0,
                                                          as.numeric(mean_strain_bootstrap_support)))
           return(
             list(
               HC_group = group_id,
               cutoff = gsnp_cutoff,
               discrepancy = discrepancy,
               before_correction_singletons = before_correction_singletons,
               before_correction_clones = before_correction_clones,
               clones_corrected = clones_corrected,
               after_correction_singletons = after_correction_singletons,
               after_correction_clones = after_correction_clones,
               median_strain_bootstrap_support = median_strain_bootstrap_support,
               mean_strain_bootstrap_support = mean_strain_bootstrap_support,
               strain_composition = strain_composition
             )
           )
         })
}


get_determine_strains_input <- cmpfun(get_determine_strains_input)



# Setup parallel processing
cl <- makeCluster(detectCores()-2)

clusterExport(cl,
              c("hierarchical_clustering_groups",
                "group_tree_path",
                "group_trees",
                "study_snp_matrix_path",
                "study_global_snp_matrix_path",
                "get_determine_strains_input"))

clusterEvalQ(cl, {
  library(TreeTools)
  library(ape)
  library(cluster)
  library(dplyr)
  library(purrr)
})

determine_strains_input <- parLapplyLB(cl,
                                       group_trees,
                                       function(group_tree){
                                         get_determine_strains_input(hierarchical_clustering_groups,
                                                                     group_tree_path,
                                                                     group_tree,
                                                                     study_snp_matrix_path,
                                                                     study_global_snp_matrix_path)
                                       })
stopCluster(cl)
rm(cl)

# Endpoints for determining strains -----
## Plateau ----
get_plateau_strains <- function(determine_strains_input,
                                hierarchical_clustering_groups,
                                plateau_length){
  
  group_output <- lapply(determine_strains_input, function(group_input){
    
    group_id <- group_input[[1]]$HC_group
    
    # Compare the strain compositions and get an array showing if identical with previous cutoff
    
    rle_array <- sapply(group_input[-1], function(cutoff_composition) {
      
      cutoff <- cutoff_composition$cutoff
      
      get_sorted_composition <- function(target_cutoff) {
        strain_data <- group_input[sapply(group_input, `[[`, "cutoff") == target_cutoff]
        composition <- lapply(strain_data[[1]]$strain_composition, `[[`, "genome")
        composition[order(sapply(composition, length))]
      }
      
      current_composition <- get_sorted_composition(cutoff)
      previous_composition <- get_sorted_composition(cutoff - 1)
      
      if (length(current_composition) != length(previous_composition)) {
        return(0)
      }
      
      current_sorted <- lapply(current_composition, sort)
      previous_sorted <- lapply(previous_composition, sort)
      
      return(if(identical(current_sorted, previous_sorted)) 1 else 0)
    })
    
    names(rle_array) <- sapply(group_input[-1], `[[`, "cutoff")
    
    
    # Find the earliest plateau
    
    plateau <- rle(rle_array)
    plateau_pos <- which(plateau$values == 1 & plateau$lengths >= plateau_length)[1]
    
    plateau_cutoff <- if(!is.na(plateau_pos)){
      if (plateau_pos == 1) as.integer(names(rle_array[1]))
      else as.integer(names(rle_array)[sum(plateau$lengths[1:(plateau_pos-1)]) + 1])
    }else "No Plateau Found"
    
    # Pull out the plateau strains
    
    plateau_strains <- if(plateau_cutoff != "No Plateau Found"){
      group_input[which(sapply(group_input, `[[`, "cutoff") == plateau_cutoff)][[1]]
    }else NULL
    
    return(list(
      group = group_id,
      plateau = plateau_cutoff,
      plateau_length = plateau_length,
      composition = plateau_strains
    ))
  })
  
  # A data frame describing the plateaus for groups
  
  plateau_df <- do.call(rbind,
                        lapply(group_output,
                        function(group){
                          data.frame(
                            group = group$group,
                            plateau = group$plateau,
                            plateau_length = plateau_length
                          )
                        })) %>%
    arrange(group)
  
  # Use group_output to generate the final strain data frame
  
  # Groups >= 2 genomes 
  strain_df <- do.call(rbind,
                       lapply(group_output,
                              function(group){
                                
                                group_strain_composition <- lapply(group$composition$strain_composition,
                                                                   `[[`,
                                                                   "genome")
                                return(do.call(rbind,
                                               lapply(seq_along(group_strain_composition),
                                                      function(strain_entry){
                                                        return(
                                                          data.frame(
                                                            strain_id = paste0(group$group,
                                                                               "_",
                                                                               strain_entry),
                                                            genome = unlist(group_strain_composition[[strain_entry]])
                                                          )
                                                        )
                                                      })))
                                
                              }))
  
  # Groups with only 1 genome
  singleton_groups <-  hierarchical_clustering_groups[sapply(hierarchical_clustering_groups,
                                                                           function(group){
                                                                             length(group$genomes) == 1
                                                                           })]
  strain_df <- rbind(strain_df,
                     do.call(rbind,
                       lapply(singleton_groups,
                              function(singleton_group){
                                return(
                                  data.frame(
                                    strain_id = paste0(singleton_group$hc_group,
                                                       "_1"),
                                    genome = unlist(singleton_group$genomes)
                                  )
                                )
                              })))
    
  
    
  return(list(
    strains = strain_df,
    plateaus = plateau_df,
    plateau_length = plateau_length,
    composition_details = group_output
  ))
  
}

final_strains <- get_plateau_strains(determine_strains_input,
                                     hierarchical_clustering_groups,
                                     plateau_length)


