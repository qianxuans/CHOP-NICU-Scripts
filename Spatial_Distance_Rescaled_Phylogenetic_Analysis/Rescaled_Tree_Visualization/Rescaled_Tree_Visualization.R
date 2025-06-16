# Libraries ----
library(ggtree)
library(phangorn)
library(ape)
library(ggplot2)
library(patchwork)
library(lme4)
library(lmerTest)

## Input ----
# ml_tree_dir NOT PROVIDED
# bed_coord_path NOT PROVIDED
# cluster_summary NOT PROVIDED
# genome_summary NOT PROVIDED
# output_dir NOT PROVIDED
# pixel_scale NOT PROVIDED

# Function to nnls and visualize the trees ----
nnls_tree <- function(output_dir,
                      ml_tree_dir,
                      bed_coord_path,
                      genome_summary,
                      cluster_summary,
                      pixel_scale){
  ## Input ----
  setwd(dir = output_dir)
  # get the list of trees
  ml_tree_list <- list.files(ml_tree_dir,
                             pattern = "*contree",
                             all.files = TRUE,
                             full.names = TRUE,
                             recursive = FALSE)
  
  # bed coordinates
  bed_coord <- read.delim(bed_coord_path)
  
  for(cluster_index in seq_along(cluster_summary)){
    
    # Get a matrix of pair-wise geographic distance 
    cluster_genomes <- cluster_summary[[cluster_index]]$genomes
    
    cluster_dist_matrix <- matrix(0, 
                                  nrow = length(cluster_genomes),
                                  ncol = length(cluster_genomes),
                                  dimnames = list(cluster_genomes, cluster_genomes))
    
    for (subject_index in 1:length(cluster_genomes)) {
      for (query_index in 1:length(cluster_genomes)) {
        if (subject_index != query_index) { 
          
          subject_genome <- cluster_genomes[subject_index]
          subject_genome_bed <- genome_summary$bed[genome_summary$genome == subject_genome]
          query_genome <- cluster_genomes[query_index]
          query_genome_bed <- genome_summary$bed[genome_summary$genome == query_genome]
          
          # Get coordinates for subject genome
          subject_coord_x <- bed_coord$X[bed_coord$bed == subject_genome_bed]
          subject_coord_y <- bed_coord$Y[bed_coord$bed == subject_genome_bed]
          
          # Get coordinates for query genome
          query_coord_x <- bed_coord$X[bed_coord$bed == query_genome_bed]
          query_coord_y <- bed_coord$Y[bed_coord$bed == query_genome_bed]
          
          if(length(na.omit(subject_coord_x)) == 1 & length(na.omit(query_coord_x)) == 1){
            # Calculate Euclidean distance
            cluster_dist_matrix[subject_index, query_index] <- sqrt((subject_coord_x - query_coord_x)^2 + (subject_coord_y - query_coord_y)^2) * pixel_scale
          }else{
            cluster_dist_matrix[subject_index, query_index] <- 0
          }
          
        }
      }
    }
    
    # Read the ML tree for this cluster
    cluster_ml_tree_path <- grep(paste0("iqtree_cluster",cluster_summary[[cluster_index]]$cluster,".contree"),
                                 ml_tree_list,
                                 value = TRUE)
    
    if(length(cluster_ml_tree_path)  == 1){
     
      
      cluster_ml_tree <- read.tree(cluster_ml_tree_path) 
      
      # Root at the oldest isolate/genome
      cluster_ml_tree <- root(cluster_ml_tree,
                              outgroup = intersect(cluster_genomes,
                                                   genome_summary$genome[genome_summary$collection_date == cluster_summary[[cluster_index]]$first_seen])[1],
                              resolve.root = TRUE)
      # Mid-point rooting is not used for now
      # cluster_ml_tree <- phytools::midpoint_root(cluster_ml_tree)
      
      spatial_ml_tree <- nnls.tree(cluster_dist_matrix,
                                   cluster_ml_tree,
                                   balanced = TRUE)
      
      # If after correction, all branch lengths are 0
      # skip the output
      if(!all(spatial_ml_tree$edge.length == 0)){
        
        ml_tree <- ggtree(cluster_ml_tree,
                          layout = "rectangular",
                          size=0.75) + 
          theme(legend.position = "none") + 
          ggtitle(paste0("Genomic Maximum Likelihood Phylogeny"))+
          theme(
            plot.title = element_text(
              hjust = 0.5,
              size = 20
            )
          )+
          coord_cartesian(clip = 'off')
        
        ml_tree <- ml_tree +
          geom_treescale(fontsize = 8.5,
                         linesize = 1.5,
                         label = "substitutions/site")
        
        spatial_tree <- ggtree(spatial_ml_tree,
                               layout = "rectangular",
                               size=0.75) + 
          theme(legend.position = "none") + 
          ggtitle(paste0("Spatial Distance Rescaled Phylogeny"))+
          theme(
            plot.title = element_text(
              hjust = 0.5,
              size = 20
            )
          ) + 
          coord_cartesian(clip = 'off')
        
        spatial_tree <- spatial_tree +
          geom_treescale(fontsize = 8.5,
                         linesize = 1.5,
                         label = "meters")
        
        
        combined_plot <- ml_tree + spatial_tree + 
          plot_layout(ncol=2) +
          plot_annotation(
            title = paste0("Cluster",cluster_summary[[cluster_index]]$cluster),
            theme = theme(
              plot.title = element_text(
                hjust = 0.5,
                size = 30,
                face = "bold"
              )
            )
          ) + 
          coord_cartesian(clip = 'off')
        
        pdf(file= paste0("Cluster",
                         cluster_summary[[cluster_index]]$cluster,
                         "_Trees.pdf"),
            width=12.5,
            height=10)
        print(combined_plot)
        dev.off() 
      }
    }
  }
}

# Execute ----

nnls_tree(
  output_dir,
  ml_tree_dir,
  bed_coord_path,
  genome_summary,
  cluster_summary,
  pixel_scale
)
