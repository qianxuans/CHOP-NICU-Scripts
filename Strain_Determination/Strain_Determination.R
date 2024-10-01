#!/usr/bin/env Rscript
# Function to check and install required packages ####
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    cat("Installing missing packages:", paste(new_packages, collapse=", "), "\n")
    install.packages(new_packages, repos="https://cran.rstudio.com/")
  }
}
# List of required packages
required_packages <- c("TreeTools",
                       "ape",
                       "parallel",
                       "dplyr",
                       "cluster",
                       "ggplot2",
                       "scales")
#load packages
library(TreeTools)
library(parallel)
library(ape)
library(dplyr)
library(cluster)
library(ggplot2)
library(scales)
# Access command-line arguments ####
args <- commandArgs(trailingOnly = TRUE)
# Initialize variables to store flag values
sub_trees_dir <- NULL
matrix_path <- NULL
num_cores <- NULL
max_threshold <- NULL
plateau_length <- NULL
groups <- NULL
output_dir <- NULL
# Loop through arguments and assign values to the appropriate variables
for (i in seq_along(args)) {
  if (args[i] == "--subtrees") {
    sub_trees_dir <- as.character(args[i + 1])
  } else if (args[i] == "--matrix") {
    matrix_path <- as.character(args[i + 1])
  } else if (args[i] == "--cores") {
    num_cores <- as.integer(args[i + 1])
  } else if (args[i] == "--max_threshold") {
    max_threshold <- as.integer(args[i + 1])
  } else if (args[i] == "--plateau_length") {
    plateau_length <- as.integer(args[i + 1])
  } else if (args[i] == "--groups") {
    groups <- as.character(args[i + 1])
  } else if (args[i] == "--output") {
    output_dir <- as.character(args[i + 1])
  }
}

# Check if both arguments are provided
if (is.null(sub_trees_dir) || is.null(matrix_path) || is.null(max_threshold) || is.null(plateau_length) || is.null(groups)) {
  warning("Missing required arguments.--groups, --subtrees, --matrix , --max_threshold, --plateau_length, and --output must all be provided.")
  quit(status = 1)
}

#the list of the sub-trees
group_tree_list <- data.frame(path = list.files(path = sub_trees_dir,
                                                pattern = "*contree",
                                                all.files = TRUE,
                                                full.names = TRUE,
                                                recursive = FALSE),
                              group = NA)

for(n in seq_len(nrow(group_tree_list))){
  group_tree_list$group[n] <- gsub("Group|.contree",
                                   "",
                                   strsplit(group_tree_list$path[n],
                                            split = "/")[[1]][length(strsplit(group_tree_list$path[n],
                                                                              split = "/")[[1]])])
}
remove(n)

#Read the SNP matrix
study_snp_df <- readRDS(matrix_path)
#groups determined by "Hierarchical_Clustering.R"
hc_groups_df <- read.csv(groups,
                         header = TRUE)
hc_groups_stats_df <- as.data.frame(table(hc_groups_df$group))

colnames(hc_groups_stats_df) <- c("group",
                                  "genome")
#For the genome(s) that is(are) not in any group divided using hierarchical clustering
#They are definitely singletons
hc_singetons <- hc_groups_df$genome[hc_groups_df$group %in% hc_groups_stats_df$group[hc_groups_stats_df$genome == 1]]

for(i in seq_len(nrow(hc_groups_stats_df))){
  group = hc_groups_stats_df$group[i]
  group_isolates <- hc_groups_df$genome[hc_groups_df$group == group]
  #Study-genome SNP distance matrix within this group
  if(hc_groups_stats_df$genome[i] > 1){
    group_study_snp_df <- study_snp_df %>% filter(subject %in% group_isolates &
                                                    query %in% group_isolates)
    if(unique(min(group_study_snp_df$snp)) >= 500){
      hc_groups_stats_df$over_snp_limit[i] <- 1
    }else{
      hc_groups_stats_df$over_snp_limit[i] <- 0
    }
    remove(group_study_snp_df)
  }else{
    hc_groups_stats_df$over_snp_limit[i] <- 0
  }
  remove(group)
  remove(group_isolates)
}
remove(i)
#For the isolates in the groups with over_snp_limit = 1
#They are considered singletons at this step 
hc_singetons <- c(hc_singetons,
                  hc_groups_df$genome[hc_groups_df$group %in% hc_groups_stats_df$group[hc_groups_stats_df$over_snp_limit == 1]])
#summarize the strains
plateau_strain_sum_df <- data.frame()
group_snp_cutoff <- data.frame()
# Analysis####
group_tree_list <- group_tree_list %>% filter(group %in% hc_groups_stats_df$group[hc_groups_stats_df$over_snp_limit == 0 & hc_groups_stats_df$genome >= 2])
remove(hc_groups_stats_df)
#use parLapplyLB to perform parallel
#make the cluster object with all available cores
#but leave one core for the system
if(is.null(num_cores)){
  cl <- makeCluster(detectCores())
}else{
  cl <- makeCluster(num_cores)
}

#export the dplyr library to the cluster object
clusterEvalQ(cl, {
  library(TreeTools)
  library(ape)
  library(dplyr)
  library(cluster)
  library(ggplot2)
  library(scales)
})
#export the data needed for the parallel process
clusterExport(cl,
              c("group_tree_list",
                "hc_groups_df",
                "study_snp_df",
                "max_threshold",
                "plateau_length",
                "output_dir"),
              envir = environment())
#Function of determining the strains
#Used in the parallel process by parLapplyLB
DetermineStrains <- function(hci) {
  #Phase0: Summarize sub-tree and import SNP-distance matrix####
  group <- group_tree_list$group[hci]
  #There are global genomes in the sub-tree of this group as we introduced top 3 of each study genome
  #Here we are just referring group_isolates to the study genomes
  group_isolates <- hc_groups_df$genome[hc_groups_df$group == group]
  
  #SNP-distance matrix
  #Study-genome SNP distance matrix within this group
  group_study_snp_df <- study_snp_df %>% filter(subject %in% group_isolates &
                                                  query %in% group_isolates)
  #Summarize the nodes of sub-tree of this group
  group_tree_newick <- read.tree(group_tree_list$path[hci])
  group_tree_newick <- phytools::midpoint_root(group_tree_newick)
  group_tree_newick <- Preorder(group_tree_newick)
  group_tree_newick_sum <- as.data.frame(matrix(nrow = length(group_tree_newick$node.label),
                                                ncol = 4))
  
  colnames(group_tree_newick_sum) <- c("clade",
                                       "node",
                                       "num_taxa",
                                       "taxa")
  
  group_tree_newick_sum$clade <- seq_len(group_tree_newick$Nnode)
  
  node_list <- (length(group_tree_newick$tip.label)+1) : (length(group_tree_newick$tip.label) + group_tree_newick$Nnode)
  
  for(i in seq_len(length(node_list))){
    sub_group_tree_newick <- Subtree(group_tree_newick,node_list[i])
    group_tree_newick_sum$node[i] <- node_list[i]
    group_tree_newick_sum$num_taxa[i] <- length(sub_group_tree_newick$tip.label)
    group_tree_newick_sum$taxa[i] <- paste(sub_group_tree_newick$tip.label,
                                           collapse = ",")
    remove(sub_group_tree_newick)
  }
  remove(i)
  remove(node_list)
  #The SNP cutoff tested
  strain_cutoff <- data.frame(snp_cutoff = round(min(group_study_snp_df$snp)+0.1):max_threshold,
                              snp_only_clones = NA,
                              snp_only_singletons = NA,
                              phylogeny_discrepancy = NA,
                              phylogeny_breaker = NA,
                              phylogeny_divided_strains = NA,
                              phylogeny_divided_isolates = NA)
  
  sum_clone_num <- data.frame(snp_cutoff = round(min(group_study_snp_df$snp)+0.1):max_threshold,
                              before_tree_singletons = NA,
                              before_tree_clones = NA,
                              after_tree_singletons = NA,
                              after_tree_clones = NA,
                              clones_broken = NA)
  
  #get a list to store the strain information(what genomes belong to what strains at different snp cutoff)
  #the name of the elements(sub-list) within strain_genome_list is the snp cutoff
  strain_genome_list <- NULL
  
  for(snp_cutoff in round(min(group_study_snp_df$snp)+0.1):max_threshold){
    #Phase 1: SNP only####
    #in the redundant_strain_df, the strain isolates(non-singleton) might appear more than once 
    #with different temporary strain id
    #this will be handled later
    redundant_strain_df <- as.data.frame(matrix(nrow = 0,
                                                ncol = 2))
    
    colnames(redundant_strain_df) <- c("strain",
                                       "isolate")
    
    #All the possible links among the study genomes at this SNP cutoff
    group_study_snp_df_cutoff <- group_study_snp_df %>% filter(snp <= snp_cutoff)
    #Within those possible links, the isolates below would be connected and considered the strain isolates
    #Which means that those the isolates in strain_isolate will not be singletons
    strain_isolate <- data.frame(isolate = unique(c(group_study_snp_df_cutoff$subject,
                                                    group_study_snp_df_cutoff$query)))
    
    #If the logical below is true (length(setdiff(group_isolates,strain_isolate$isolate) > 0)), 
    #this means that at this snp cutoff, 
    #there are singleton isolates and strain isolates
    #if false,
    #at this snp, all isolates are identified as strain isolate 
    #and no singleton isolate exists at this point
    #so the data frame singleton_isolate will be one with 0 row as a placeholder
    
    if(length(setdiff(group_isolates,strain_isolate$isolate) > 0)){
      singleton_isolate <- data.frame(isolate = setdiff(group_isolates,
                                                        strain_isolate$isolate))
    }else{
      singleton_isolate <- as.data.frame(matrix(nrow = 0,
                                                ncol = 1))
      colnames(singleton_isolate) <- "isolate"
    }
    
    #Set a character variable as a place holder for strain number(tmp_strain_num)
    #which starts from 1
    #Every iteration, there will be 1 added to the tmp_strain_num
    tmp_strain_num <- 1
    for(i in 1:nrow(strain_isolate)){
      #For this isolate (isolate[i])
      isolate <- strain_isolate$isolate[i]
      #Find the links associated with this isolate
      isolate_links <- group_study_snp_df_cutoff %>% filter(group_study_snp_df_cutoff$subject == isolate | group_study_snp_df_cutoff$query == isolate)
      #Every genomes involved in those links associated with this isolate
      #will be given tmp_strain_num as the temporary ID
      redundant_strain_df_tmp <- data.frame(strain = tmp_strain_num,
                                            isolate = unique(c(isolate_links$subject,
                                                               isolate_links$query)))
      #rbind to the redundant_strain_df
      redundant_strain_df <- rbind(redundant_strain_df,
                                   redundant_strain_df_tmp)
      
      #add 1 to the temporary ID
      tmp_strain_num <- tmp_strain_num + 1
      
      remove(isolate)
      remove(isolate_links)
      remove(redundant_strain_df_tmp)
    }
    
    remove(tmp_strain_num)
    remove(group_study_snp_df_cutoff)
    
    #Check if there is any overlapping of isolates in redundant_strain_df
    #for example, one isolate might have multiple IDs, 
    #if so, the link is found to merge multiple strains into the same strain
    #Note: 
    #The isolates in the unique_strain_df$isolate are the strain isolates
    #The singleton isolates will be added once the unique strain IDs are assigned
    
    #The real strain ID also starts at 1
    unique_strain_df <- data.frame(strain = NA,
                                   isolate = unique(redundant_strain_df$isolate))
    real_strain_id <- 1
    
    #loop from the first isolate to the last in the list of strain_isolate
    #to find the potential link to merge multiple strains into one
    
    for(i in 1:nrow(strain_isolate)){
      #the name of the isolate
      isolate <- strain_isolate$isolate[i]
      #the isolate is involved in what redundant strains
      redundant_strain <- redundant_strain_df$strain[redundant_strain_df$isolate == isolate]
      #any isolates in the isolate_redundant_strain are considered from the same strain
      #the isolates from the redundant_strain:
      redundant_strain_isolate <- unique(redundant_strain_df$isolate[redundant_strain_df$strain %in% redundant_strain])
      
      #see if the strain of those isolates is already labeled with real_strain_id, if not, label the strain with real_strain_id, and add 1 to real_strain_id
      #if so, skip to next unique(redundant_strain_df$strain)
      #if any one of the isolates is already labelled a read_strain_id, the rest of the isolates, which are not labeled, are labeled the same read_strain_id
      if(sum(is.na(unique_strain_df$strain[unique_strain_df$isolate %in% redundant_strain_isolate])) == 
         length(unique_strain_df$strain[unique_strain_df$isolate %in% redundant_strain_isolate])){
        #situation 1: none of the isolates are labeled with real_strain_id (all in strain column in unique_strain_df is "NA")
        unique_strain_df$strain[unique_strain_df$isolate %in% redundant_strain_isolate] <- real_strain_id
        real_strain_id <- real_strain_id + 1
      }else if(sum(is.na(unique_strain_df$strain[unique_strain_df$isolate %in% redundant_strain_isolate])) < 
               length(unique_strain_df$strain[unique_strain_df$isolate %in% redundant_strain_isolate]) &
               sum(is.na(unique_strain_df$strain[unique_strain_df$isolate %in% redundant_strain_isolate])) > 0){
        #situation 2: some of them are labeled(some in strain column in unique_strain_df is "NA", while some are already labeled)
        #if there are multiple strain numbers in the unique_strain_df, use the smallest one.
        unique_strain_df$strain[unique_strain_df$isolate %in% redundant_strain_isolate] <- min(unique(na.omit(unique_strain_df$strain[unique_strain_df$isolate %in% redundant_strain_isolate])))
      }else if(sum(is.na(unique_strain_df$strain[unique_strain_df$isolate %in% redundant_strain_isolate])) == 0){
        #situation 3: all are labeled(no NA in strain column in unique_strain_df)
        #do nothing
      }
      remove(isolate)
    }
    remove(redundant_strain_df)
    remove(redundant_strain)
    remove(redundant_strain_isolate)
    remove(i)
    remove(real_strain_id)
    
    #The singleton isolates exist, add the singletons to the unique_strain_df
    #But for now, strain ID would not be assigned to singleton isolates
    if(nrow(singleton_isolate) > 0){
      unique_strain_df_singleton <- data.frame(strain = "Singleton",
                                               isolate = singleton_isolate$isolate)
      unique_strain_df <- rbind(unique_strain_df,
                                unique_strain_df_singleton)
      remove(unique_strain_df_singleton)
    }
    remove(strain_isolate)
    remove(singleton_isolate)
    
    #record the numbers of singletons/ clones and clusters
    sum_clone_num$before_tree_singletons[sum_clone_num$snp_cutoff == snp_cutoff] <- sum(unique_strain_df$strain == "Singleton")
    sum_clone_num$before_tree_clones[sum_clone_num$snp_cutoff == snp_cutoff] <- max(as.numeric(unique_strain_df$strain[unique_strain_df$strain != "Singleton"]))
    strain_cutoff$snp_only_singletons[strain_cutoff$snp_cutoff == snp_cutoff] <- sum(unique_strain_df$strain == "Singleton")
    strain_cutoff$snp_only_clones[strain_cutoff$snp_cutoff == snp_cutoff] <- max(as.numeric(unique_strain_df$strain[unique_strain_df$strain != "Singleton"]))

    #Phase 2: Phylogeny####
    #Only for clones(strains that have no less than 2 isolates, which is not singleton strains)
    #Find the strain number of the strain with no less than 2 isolates (clones)
    clone_phylogeny <- data.frame(strain = unique(unique_strain_df$strain[unique_strain_df$strain != "Singleton"]),
                                  num_strain_isolates = NA,
                                  strain_isolates = NA,
                                  best_clade = NA,
                                  best_clade_isolates = NA,
                                  num_discrepancy = NA,
                                  discrepancy_isolates = NA,
                                  num_breaker = NA,
                                  breakers = NA)
    
    for(i in seq_len(nrow(clone_phylogeny))){
      
      strain_num <- clone_phylogeny$strain[i]
      
      strain_isolates <- unique_strain_df$isolate[unique_strain_df$strain == strain_num]
      
      clone_phylogeny$strain_isolates[i] <- paste(strain_isolates,
                                                  collapse = ",")
      
      clone_phylogeny$num_strain_isolates[i] <- length(strain_isolates)
      
      clade_cover <- as.data.frame(matrix(nrow = 0,
                                          ncol = 9))
      
      colnames(clade_cover) <- c("clade",
                                 "num_taxa_in_clade",
                                 "taxa_in_clade",
                                 "num_breaker_included_taxa_in_clade",
                                 "breaker_included_taxa_in_clade",
                                 "num_breaker_taxa_in_clade",
                                 "breaker_taxa_in_clade",
                                 "strain0_clade1_num",
                                 "strain0_clade1")
      
      for(n in seq_len(nrow(group_tree_newick_sum))){
        clade_num <- group_tree_newick_sum$clade[n]
        
        clade_isolates <- strsplit(group_tree_newick_sum$taxa[n],
                                   split = ",")[[1]]
        
        #If length(setdiff(isolates,clade_isolate)) == 0,
        #this means that all the study isolates in this strain are covered in the clade
        #only when all the study isolates in this strain are covered, the clade is counted
        if(length(setdiff(strain_isolates,clade_isolates)) == 0){
          
          clade_cover_tmp <- as.data.frame(matrix(nrow = 1,
                                                  ncol = 9))
          
          colnames(clade_cover_tmp) <- c("clade",
                                         "num_taxa_in_clade",
                                         "taxa_in_clade",
                                         "num_breaker_included_taxa_in_clade",
                                         "breaker_included_taxa_in_clade",
                                         "num_breaker_taxa_in_clade",
                                         "breaker_taxa_in_clade",
                                         "strain0_clade1_num",
                                         "strain0_clade1")
          
          clade_cover_tmp$clade <- clade_num
          #this step counts both study and brekaer(global) genomes
          clade_cover_tmp$num_breaker_included_taxa_in_clade <- length(clade_isolates)
          clade_cover_tmp$breaker_included_taxa_in_clade <- group_tree_newick_sum$taxa[n]
          #this only counts the study genomes
          #no breaker(global) genomes
          breaker_excluded_taxa_in_clade <- strsplit(group_tree_newick_sum$taxa[n],
                                                     split = ",")[[1]][!grepl("GCA",
                                                                              strsplit(group_tree_newick_sum$taxa[n],
                                                                                       split = ",")[[1]])]
          clade_cover_tmp$taxa_in_clade <- paste(breaker_excluded_taxa_in_clade,
                                                 collapse = ",")
          
          clade_cover_tmp$num_taxa_in_clade <- length(breaker_excluded_taxa_in_clade)
          #only count the breaker genomes in the clade
          clade_cover_tmp$num_breaker_taxa_in_clade <- length(strsplit(group_tree_newick_sum$taxa[n],
                                                                       split = ",")[[1]][grepl("GCA",
                                                                                               strsplit(group_tree_newick_sum$taxa[n],
                                                                                                        split = ",")[[1]])])
          
          clade_cover_tmp$breaker_taxa_in_clade <- paste(strsplit(group_tree_newick_sum$taxa[n],
                                                                  split = ",")[[1]][grepl("GCA",
                                                                                          strsplit(group_tree_newick_sum$taxa[n],
                                                                                                   split = ",")[[1]])],
                                                         collapse = ",")
          
          #Counts the study genomes that are in this clade but NOT in strain defined by SNP only cutoff
          clade_cover_tmp$strain0_clade1_num <- length(setdiff(breaker_excluded_taxa_in_clade,
                                                               strain_isolates))
          
          clade_cover_tmp$strain0_clade1 <- paste(setdiff(breaker_excluded_taxa_in_clade,
                                                          strain_isolates),
                                                  collapse = ",")
          
          clade_cover <- rbind(clade_cover,
                               clade_cover_tmp)
          
          remove(clade_cover_tmp)
          remove(breaker_excluded_taxa_in_clade)
        }
        remove(clade_num)
        remove(clade_isolates)
      }
      remove(n)
      clone_phylogeny$best_clade[i] <- clade_cover$clade[which.min(clade_cover$num_breaker_included_taxa_in_clade)]
      #in best_clade_isolate, only the NICU genome(s) being counted
      clone_phylogeny$best_clade_isolates[i] <- paste(strsplit(clade_cover$taxa_in_clade[which.min(clade_cover$num_breaker_included_taxa_in_clade)],
                                                               split = ",")[[1]],
                                                      collapse = ",")
      
      clone_phylogeny$num_discrepancy[i] <- clade_cover$strain0_clade1_num[which.min(clade_cover$num_breaker_included_taxa_in_clade)]
      clone_phylogeny$discrepancy_isolates[i] <- clade_cover$strain0_clade1[which.min(clade_cover$num_breaker_included_taxa_in_clade)]
      
      #count the number of breaker genomes in the clades
      clone_phylogeny$num_breaker[i] <- clade_cover$num_breaker_taxa_in_clade[which.min(clade_cover$num_breaker_included_taxa_in_clade)]
      clone_phylogeny$breakers[i] <- clade_cover$breaker_taxa_in_clade[which.min(clade_cover$num_breaker_included_taxa_in_clade)]
      #clean
      remove(clade_cover)
      remove(strain_num)
      remove(strain_isolates)
    }
    remove(i)
    
    #Count the number of discrepancy genomes and the breakers
    discrepancy_genomes <- NULL
    breakers <- NULL
    for(i in seq_len(nrow(clone_phylogeny))){
      #count discrepancy genomes 
      discrepancy_genomes_tmp <- strsplit(clone_phylogeny$discrepancy_isolates[i],
                                          split = ",")[[1]]
      
      if(length(discrepancy_genomes_tmp) > 0){
        discrepancy_genomes <- unique(c(discrepancy_genomes,
                                        discrepancy_genomes_tmp))
      }
      remove(discrepancy_genomes_tmp)
      
      #count breaker genomes 
      breakers_tmp <- strsplit(clone_phylogeny$breakers[i],
                               split = ",")[[1]]
      
      if(length(breakers_tmp) > 0){
        breakers <- unique(c(breakers,
                             breakers_tmp))
      }
      remove(breakers_tmp)
    }
    
    strain_cutoff$phylogeny_discrepancy[strain_cutoff$snp_cutoff == snp_cutoff] <- length(discrepancy_genomes)
    strain_cutoff$phylogeny_breaker[strain_cutoff$snp_cutoff == snp_cutoff] <- length(breakers)
    
    remove(i)
    remove(discrepancy_genomes)
    remove(breakers)
    
    #Phase 3: SNP + Phylogeny ####
    #use the discrepancy found in clone_phylogeny to add the discrepancy genomes back to the unique_strains_df
    #also use the breaker genomes to decide what genomes are not supposed to be the same strain because of setting SNP cutoff too high
    #This step only fixes strains with either discrepancy or breaker genomes in the clade inferred by the tree
    #For the strains filtered by this step, the SNP-only identification of strain is in full compliance with the sub-tree
    clone_phylogeny <- clone_phylogeny %>% filter(num_discrepancy > 0 |
                                                    num_breaker > 0)
    
    clone_phylogeny <- clone_phylogeny %>% arrange(num_strain_isolates)
    
    for(i in seq_len(nrow(clone_phylogeny))){
      #In this iteration, 2 situations are expected:
      if(clone_phylogeny$num_breaker[i] == 0){
        #1st situation: only discrepancy genomes with no breaker genomes
        #get a data frame of the new isolates of the new strain,
        #This step will also check if any discrepancy isolates are already labeled another strain
        new_strain_isolates_df <- data.frame(strain = NA,
                                             new_isolates = strsplit(clone_phylogeny$best_clade_isolate[i],
                                                                     split = ",")[[1]])
        for(s in seq_len(nrow(new_strain_isolates_df))){
          new_strain_isolates_df$strain[s] <- unique_strain_df$strain[unique_strain_df$isolate == new_strain_isolates_df$new_isolates[s]]
        }
        remove(s)
        #pick the smallest id as the new ID for all isolates in new_strain_isolates_df
        new_id <- min(as.numeric(new_strain_isolates_df$strain[new_strain_isolates_df$strain != "Singleton"]))
        #use the new_id to replace all 
        unique_strain_df$strain[unique_strain_df$isolate %in% new_strain_isolates_df$new_isolates] <- new_id
        #clean
        remove(new_strain_isolates_df)
        remove(new_id)
        
      }else if(clone_phylogeny$num_breaker[i] > 0){
        #2nd situation: breaker genomes were introduced (with or without discrepancy genomes)
        #Summarize the sub-tree of this clone defined by SNP cutoff
        clone_phylogeny_breaker_subtree <- Subtree(group_tree_newick,
                                                   group_tree_newick_sum$node[group_tree_newick_sum$clade == clone_phylogeny$best_clade[i]])
        
        clone_phylogeny_breaker_subtree_sum <- as.data.frame(matrix(nrow = length(clone_phylogeny_breaker_subtree$node.label),
                                                                    ncol = 6))
        
        colnames(clone_phylogeny_breaker_subtree_sum) <- c("clade",
                                                           "node",
                                                           "num_taxa",
                                                           "taxa",
                                                           "if_pure_study",
                                                           "if_biggest_pure_study")
        #Keep track of study genomes that are counted
        #What's left are the isolated study genomes divided by the cladebraker genomes
        processed_study_genomes <- c()
        
        clone_phylogeny_breaker_subtree_sum$clade <- seq_len(clone_phylogeny_breaker_subtree$Nnode)
        
        clone_phylogeny_breaker_subtree_node_list <- (length(clone_phylogeny_breaker_subtree$tip.label)+1) : (length(clone_phylogeny_breaker_subtree$tip.label) + clone_phylogeny_breaker_subtree$Nnode)
        #prune the clone_phylogeny_breaker_subtree to any possible sub trees
        for(nli in seq_len(length(clone_phylogeny_breaker_subtree_node_list))){
          
          clone_phylogeny_breaker_subtree_subtree <- Subtree(clone_phylogeny_breaker_subtree,clone_phylogeny_breaker_subtree_node_list[nli])
          
          clone_phylogeny_breaker_subtree_sum$node[nli] <- clone_phylogeny_breaker_subtree_node_list[nli]
          
          clone_phylogeny_breaker_subtree_sum$num_taxa[nli] <- length(clone_phylogeny_breaker_subtree_subtree$tip.label)
          
          clone_phylogeny_breaker_subtree_sum$taxa[nli] <- paste(clone_phylogeny_breaker_subtree_subtree$tip.label,
                                                                 collapse = ",")
          #check if this sub-clade contains purely study genomes
          #if this clade has purely study genomes that are considered from the same strain by using SNP + Phylogeny method
          if(all(clone_phylogeny_breaker_subtree_subtree$tip.labe %in% strsplit(clone_phylogeny$best_clade_isolates[i],split = ",")[[1]])){
            #if yes(1),
            #this clade is a pure study genome clade
            clone_phylogeny_breaker_subtree_sum$if_pure_study[nli] <- 1
            #but if this is the biggest one without being broken by cladebreaker genomes?
            #to answer this question, go to the parent clade to check if the parent clade is pure study-genome clade
            #if the parent clade is still pure study-genome clade, this clade is not the biggest study-only clade 
            parent_node <- clone_phylogeny_breaker_subtree$edge[which(clone_phylogeny_breaker_subtree$edge[,2] == clone_phylogeny_breaker_subtree_node_list[nli]), 1]
            clone_phylogeny_breaker_subtree_subtree_parent <- Subtree(clone_phylogeny_breaker_subtree,parent_node)
            if(all(clone_phylogeny_breaker_subtree_subtree_parent$tip.label %in% strsplit(clone_phylogeny$best_clade_isolates[i],split = ",")[[1]])){
              clone_phylogeny_breaker_subtree_sum$if_biggest_pure_study[nli] <- 0
            }else{
              #now we find the one of the biggest study-genome clade broken by the cladebreaker genomes 
              clone_phylogeny_breaker_subtree_sum$if_biggest_pure_study[nli] <- 1
              #put the study genomes in this clade to processed_study_genomes 
              processed_study_genomes <- c(processed_study_genomes,
                                           clone_phylogeny_breaker_subtree_subtree$tip.label)
            }
            remove(clone_phylogeny_breaker_subtree_subtree_parent) 
            
            remove(parent_node)
          }else{
            #if no(0), 
            #this is not a pure study genome clade
            clone_phylogeny_breaker_subtree_sum$if_pure_study[nli] <- 0
            clone_phylogeny_breaker_subtree_sum$if_biggest_pure_study[nli] <- 0
          }
          #clean
          remove(clone_phylogeny_breaker_subtree_subtree)
        }
        remove(nli)
        remove(clone_phylogeny_breaker_subtree_node_list)
        #only keep the biggest pure-study-genome clades
        clone_phylogeny_breaker_subtree_sum <- clone_phylogeny_breaker_subtree_sum %>% 
          filter(if_biggest_pure_study == 1)
        #sort the results of cladebreaker genomes
        #new_strain_isolates_df is used to give the ID to the broken strains
        #e.g. strain1-divided1, strain1-divided2, strain1-divided3, ... 
        
        new_strain_isolates_df <- data.frame()
        
        for(nli in seq_len(nrow(clone_phylogeny_breaker_subtree_sum))){
          
          new_strain_isolates_df_tmp <- data.frame(strain = paste0(clone_phylogeny$strain[i],
                                                                   "-divided-",
                                                                   nli),
                                                   new_isolates = strsplit(clone_phylogeny_breaker_subtree_sum$taxa[nli],
                                                                           split = ",")[[1]])
          new_strain_isolates_df <- rbind(new_strain_isolates_df,
                                          new_strain_isolates_df_tmp)
          remove(new_strain_isolates_df_tmp)
        }
        #look at the isolated genomes
        clone_phylogeny_breaker_isolated_genomes <- setdiff(strsplit(clone_phylogeny$best_clade_isolates[i],
                                                                     split = ",")[[1]],
                                                            processed_study_genomes)
        
        for(nls in seq_len(length(clone_phylogeny_breaker_isolated_genomes))){
          id <- nls + nli
          
          new_strain_isolates_df_tmp <- data.frame(strain = paste0(clone_phylogeny$strain[i],
                                                                   "-divided-",
                                                                   id,
                                                                   "-singleton"),
                                                   new_isolates = clone_phylogeny_breaker_isolated_genomes[nls])
          
          new_strain_isolates_df <- rbind(new_strain_isolates_df,
                                          new_strain_isolates_df_tmp)
          remove(new_strain_isolates_df_tmp)
          remove(id)
        }
        #clean
        remove(processed_study_genomes)
        remove(clone_phylogeny_breaker_isolated_genomes)
        remove(nli)
        remove(nls)
        #use the info in the new_strain_isolates_df to modify unique_strain_df
        for(s in seq_len(nrow(new_strain_isolates_df))){
          unique_strain_df$strain[unique_strain_df$isolate == new_strain_isolates_df$new_isolates[s]] <- new_strain_isolates_df$strain[s]
        }
        #clean
        remove(s)
        remove(new_strain_isolates_df)
        remove(clone_phylogeny_breaker_subtree)
        remove(clone_phylogeny_breaker_subtree_sum)
      }
    }
    remove(clone_phylogeny)
    remove(i)
    #Record after-tree singletons:
    #For the singletons inferred by the cladebreaker,
    #Count them as singleton
    sum_clone_num$after_tree_singletons[sum_clone_num$snp_cutoff == snp_cutoff] <- sum(unique_strain_df$strain == "Singleton") + sum(grepl("-singleton",unique_strain_df$strain))
    #record after-tree clones:
    clone_count <- unique(unique_strain_df$strain)
    #Count the clones broken by the cladebreaker genomes
    clone_broken <- c()
    for(id in seq_len(length(clone_count))){
      if(grepl("-divided-",clone_count[id])){
        clone_broken <- c(clone_broken,
                          strsplit(clone_count[id],split = "-")[[1]][1])
        
        clone_broken <- unique(clone_broken)
      }
    }
    remove(id)
    #Remove Singletons from the strain ID to count the number of clones
    clone_count <- clone_count[clone_count != "Singleton"]
    clone_count <- clone_count[!grepl("-singleton",clone_count)]
    
    sum_clone_num$after_tree_clones[sum_clone_num$snp_cutoff == snp_cutoff] <- length(unique(clone_count))
    sum_clone_num$clones_broken[sum_clone_num$snp_cutoff == snp_cutoff] <- length(unique(clone_broken))
    strain_cutoff$phylogeny_divided_strains[strain_cutoff$snp_cutoff == snp_cutoff] <- length(unique(clone_broken))
    remove(clone_count)
    remove(clone_broken)
    #add the strain genome information to strain_genome_list 
    strain_genome_list_tmp <- NULL
    #get a data frame to sort the info to be added to strain_genome_list 
    
    for(tmp_strain in unique(unique_strain_df$strain)){
      strain_isolates <- unique_strain_df$isolate[unique_strain_df$strain == tmp_strain]
      strain_isolates <- strain_isolates[order(strain_isolates)]
      strain_genome_list_tmp[[length(strain_genome_list_tmp)+1]] <- list(strain_isolates)
      names(strain_genome_list_tmp)[[length(strain_genome_list_tmp)]] <- tmp_strain
      #clean
      remove(strain_isolates)
    }
    strain_cutoff$phylogeny_divided_isolates[strain_cutoff$snp_cutoff == snp_cutoff] <- nrow(unique_strain_df %>% filter(grepl("-divided-",strain)))
    remove(tmp_strain)
    remove(unique_strain_df)
    #Add the strain_genome_list_tmp to strain_genome_list
    strain_genome_list[[length(strain_genome_list)+1]] <- strain_genome_list_tmp
    names(strain_genome_list)[[length(strain_genome_list)]] <- snp_cutoff
    remove(strain_genome_list_tmp)
  }
  remove(snp_cutoff)
  
  
  group_snp_cutoff_tmp <- data.frame(group = group,
                                     num_genomes = length(group_isolates),
                                     plateau = NA)
  
  #compare the list to see if the i is the same as previous 
  strain_genome_list_compare <- data.frame(cutoff = names(strain_genome_list[1:length(strain_genome_list)]),
                                           same_as_previous = NA)
  
  for(i in 2:nrow(strain_genome_list_compare)){
    if(length(setdiff(strain_genome_list[[i]], strain_genome_list[[i-1]])) == 0 & 
       length(setdiff(strain_genome_list[[i-1]], strain_genome_list[[i]])) == 0){
      strain_genome_list_compare$same_as_previous[i] <- 1
    }else{
      strain_genome_list_compare$same_as_previous[i] <- 0
    }
  }
  remove(i)
  strain_genome_list_compare$cutoff <- as.numeric(strain_genome_list_compare$cutoff)
  plateau <- rle(strain_genome_list_compare$same_as_previous)
  plateau_run <- which(plateau$values == 1 & plateau$lengths >= plateau_length)[1]
  
  if(!is.na(plateau_run)){
    if(plateau_run == 1){
      group_snp_cutoff_tmp$plateau <- strain_genome_list_compare$cutoff[1]
    }else if(plateau_run > 1){
      group_snp_cutoff_tmp$plateau <- strain_genome_list_compare$cutoff[sum(plateau$lengths[1:(plateau_run-1)]) + 1]
    }
  }else{
    group_snp_cutoff_tmp$plateau <- "No Plateau Found"
  }
  #clean
  remove(plateau)
  remove(plateau_run)
  plateau_strain_sum_df_tmp <- data.frame()
  
  
  if(group_snp_cutoff_tmp$plateau != "No Plateau Found"){
    #if the plateau is found:
    #Use the plateau to define the strains:
    plateau_strain_list_defined <- strain_genome_list[names(strain_genome_list) == as.character(group_snp_cutoff_tmp$plateau)]
  }else{
    #if no plateau is found
    #use the max_threshold as the threshold
    plateau_strain_list_defined <- strain_genome_list[names(strain_genome_list) == as.character(max_threshold)]
  }
  for(i in seq_len(length(plateau_strain_list_defined[[1]]))){
    
    plateau_strain_sum_df_tmp_tmp <- data.frame(strain = paste0(group,
                                                                "_",
                                                                names(plateau_strain_list_defined[[1]][i])),
                                                isolate = unlist(plateau_strain_list_defined[[1]][i]))
    
    plateau_strain_sum_df_tmp <- rbind(plateau_strain_sum_df_tmp,
                                       plateau_strain_sum_df_tmp_tmp)
    
    remove(plateau_strain_sum_df_tmp_tmp)
    
  }
  remove(i)
  remove(plateau_strain_list_defined)
  remove(strain_genome_list_compare)
  
  #export the number
  write.csv(sum_clone_num,
             file = paste0(output_dir,
                           "/Group",
                    group,
                    "_clones_number.csv"),
            quote = FALSE,
            row.names = FALSE)
  
  write.csv(strain_cutoff,
             file = paste0(output_dir,
                           "/Group",
                    group,
                    "_stats.csv"),
            quote = FALSE,
            row.names = FALSE)
  
  
  #visualize the number of discrepancy genomes and clones
  line_color <- c("#78C52C",
                  "#E39223")
  names(line_color) <- c("discrepancy",
                         "clone")
  plot_df_tmp1 <- strain_cutoff %>% mutate(group = "discrepancy") %>% select(snp_cutoff,
                                                                             phylogeny_discrepancy,
                                                                             group)
  colnames(plot_df_tmp1) <- c("snp",
                              "number",
                              "category")
  plot_df_tmp2 <- sum_clone_num %>% mutate(group = "clone") %>% select(snp_cutoff,
                                                                       after_tree_clones,
                                                                       group)
  
  colnames(plot_df_tmp2) <- c("snp",
                              "number",
                              "category")
  
  plot_df <- rbind(plot_df_tmp1,
                   plot_df_tmp2)
  remove(plot_df_tmp1)
  remove(plot_df_tmp2)
  plot <- ggplot(plot_df,
                 aes(x = snp,
                     y = number,
                     color = category)) +
    scale_color_manual(values = line_color,
                       name = "Number",
                       labels = c("Clones",
                                  "Discrepancy"))+
    geom_line(alpha = 0.75,
              linewidth = 2) +
    scale_x_continuous(name = "SNP Threshold",
                       breaks = seq(0, max(plot_df$snp),50),
                       limits = c(min(plot_df$snp), max(plot_df$snp))) + 
    scale_y_continuous(
      name = "Number",
      breaks = pretty_breaks(),
    ) + 
    theme(
      axis.title.y.left = element_text(colour = "black",
                                       size = 30,
                                       face = "bold"),
      axis.text.y.left = element_text(colour = "black",
                                      size = 20),
      axis.title.x = element_text(color = "black",
                                  size = 30,
                                  face = "bold"),
      axis.text.x = element_text(colour = "black",
                                 size = 20),
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      axis.ticks.length=unit(.25, "cm"),
      plot.background = element_rect(fill = "transparent"),
      panel.background = element_rect(fill = "transparent"),
      legend.title = element_text(size=20),
      legend.text = element_text(size=15),
      legend.key.size = unit(1.5, 'cm')
    )
  
  pdf(file=paste0(output_dir,
                  "/Group",
                  group,
                  "_discrepancy clones.pdf"),
      width=12.5,
      height=7.5)
  print(plot)
  dev.off()
  
  #clean
  remove(group)
  remove(group_isolates)
  remove(group_study_snp_df)
  remove(group_tree_newick_sum)
  remove(group_tree_newick)
  remove(sum_clone_num)
  remove(strain_cutoff)
  remove(plot_df)
  remove(plot)
  remove(line_color)
  #return the results
  return(list(cutoffs = group_snp_cutoff_tmp,
              plateau = plateau_strain_sum_df_tmp))
  
}
#Execute the parallel process of identifying the strains of different groups
DetermineStrains_results <- parLapplyLB(cl,
                                        seq_len(nrow(group_tree_list)),
                                        DetermineStrains)
#Export the results from the parallel process
plateau_strain_sum_df <- do.call(rbind, lapply(DetermineStrains_results, function(x) x$plateau))
group_snp_cutoff <- do.call(rbind, lapply(DetermineStrains_results, function(x) x$cutoffs))
write.csv(group_snp_cutoff,
          file = paste0(output_dir,"All_Groups_Plateau.csv"),
          quote = FALSE,
          row.names = FALSE)
#Stop the parallel and release resources
stopCluster(cl)
rm(cl)
remove(DetermineStrains_results)
gc()
# Add the singletons inferred by hierarchical clustering to the *_strain_sum_df####
plateau_strain_sum_df <- rbind(plateau_strain_sum_df,
                               data.frame(strain = "HC_Singleton",
                                          isolate = hc_singetons))
# Assign new strain IDs####
if(length(setdiff(plateau_strain_sum_df$isolate,hc_groups_df$genome)) == 0 &
   length(setdiff(hc_groups_df$genome,plateau_strain_sum_df$isolate)) == 0 ){
  
  plateau_strain_sum_df <- plateau_strain_sum_df %>% mutate(new_strain_id = NA)
  
  new_id_df <- data.frame(old_id = unique(plateau_strain_sum_df$strain[!grepl("Singleton|singleton",plateau_strain_sum_df$strain)]),
                          new_id = NA)
  new_id_df$new_id <- seq_len(nrow(new_id_df))
  
  new_id_df <- rbind(new_id_df,
                     data.frame(old_id = unique(plateau_strain_sum_df$isolate[grepl("Singleton|singleton",
                                                                                    plateau_strain_sum_df$strain)]),
                                new_id =  (nrow(new_id_df) + 1) : (nrow(new_id_df) + length(unique(plateau_strain_sum_df$isolate[grepl("Singleton|singleton",
                                                                                                                                       plateau_strain_sum_df$strain)])))))
  
  for(i in seq_len(nrow(plateau_strain_sum_df))){
    if(grepl("Singleton|singleton",plateau_strain_sum_df$strain[i])){
      plateau_strain_sum_df$new_strain_id[i] <- new_id_df$new_id[new_id_df$old_id == plateau_strain_sum_df$isolate[i]] 
    }else{
      plateau_strain_sum_df$new_strain_id[i] <- new_id_df$new_id[new_id_df$old_id == plateau_strain_sum_df$strain[i]] 
    }
  }
  remove(i)
  remove(new_id_df)
}

write.csv(plateau_strain_sum_df,
          file = paste0(output_dir,"/Strains_summary.csv"),
          quote = FALSE,
          row.names = FALSE)
