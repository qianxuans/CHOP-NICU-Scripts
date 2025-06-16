# Generate the floor plan showing the collection locations of cluster genomes within NICU
# The color gradient of the dot represents the day post-initial detection

# Libraries ----
library(ggplot2)
library(ggmap)
library(scales)
library(colorspace)
library(dplyr)
library(ggrepel)
library(openxlsx2)
library(magick)
library(phangorn)
library(ape)
library(patchwork)
library(lme4)
library(lmerTest)

# Input ----
#output_dir NOT PROVIDED
#summary_path NOT PROVIDED
#bed_coord_path NOT PROVIDED
#floorplan_png_path NOT PROVIDED
#ml_tree_dir NOT PROVIDED
#pixel_scale NOT PROVIDED

# Function to generate the floor plan ----

floorplan_cluster_distance <- function(output_dir,
                                       summary_path,
                                       bed_coor_path,
                                       ml_tree_list,
                                       floorplan_png_path,
                                       pixel_scale){
  
  ## Universal variables / Input ----
  setwd(dir = output_dir)
  ### floorplan png
  map_plot <- image_read(floorplan_png_path)
  map_plot <- image_ggplot(map_plot)
  ### Bed X Y coordinates
  bed_coord <- read.table(bed_coord_path,
                          header = TRUE)
  ### Get the list of trees
  
  ml_tree_list <- list.files(ml_tree_dir,
                             pattern = "*contree",
                             all.files = TRUE,
                             full.names = TRUE,
                             recursive = FALSE)
  ### Cluster summary
  cluster_summary <- readRDS(summary_path)[["clusters"]]
  ### Genome summary
  genome_summary <- readRDS(summary_path)[["genomes"]]
  ## Iterate through clusters ----
  
  for(cluster_index in seq_along(cluster_summary)){
    
    # Read the ML tree for this cluster
    cluster_ml_tree_path <- grep(paste0("iqtree_cluster",cluster_summary[[cluster_index]]$cluster,".contree"),
                                 ml_tree_list,
                                 value = TRUE)
    
    if(length(cluster_ml_tree_path)  == 1){
      
      cluster_ml_tree <- read.tree(cluster_ml_tree_path) 
      
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
              cluster_dist_matrix[subject_index, query_index] <- pixel_scale * sqrt((subject_coord_x - query_coord_x)^2 + (subject_coord_y - query_coord_y)^2)
            }else{
              cluster_dist_matrix[subject_index, query_index] <- 0
            }
            
          }
        }
      }
      
      # The first day detected:
      first_day <- as.Date(paste0("20",cluster_summary[[cluster_index]]$first_seen))
      # The last day detected: 
      last_day <- as.Date(paste0("20",cluster_summary[[cluster_index]]$last_seen))
      
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
      # Convert it back to distance matrix
      phylogentic_corrected_distance_matrix <- cophenetic.phylo(spatial_ml_tree)
      
      # Two data frame will be generated as the input for 2 visualizations
      ### The first is the line plot showing the accumuative and new distance this cluster spread in the NICU  ----
      if(first_day == last_day){
        distance_by_day <- rbind(
          data.frame(
          cluster = cluster_summary[[cluster_index]]$cluster,
          invasive_status =cluster_summary[[cluster_index]]$invasive_status,
          date = first_day,
          day_post_initial_detection = 0,
          distance = 0,
          distance_category = "New"
        ),
        data.frame(
          cluster = cluster_summary[[cluster_index]]$cluster,
          invasive_status =cluster_summary[[cluster_index]]$invasive_status,
          date = first_day,
          day_post_initial_detection = 0,
          distance = 0,
          distance_category = "Accumulated"
        ))
      }else{
        
        persistent_period <- c(0,seq_len(last_day - first_day))
        
        names(persistent_period) <- as.character(as.Date(first_day:last_day))
        
        distance_by_day <- data.frame()
        
        accumulated_distance <- 0
        
        for(day_entry in persistent_period){
          # Calculate the spatial distance the cluster spread this day
          day_new_genome <- intersect(cluster_genomes,
                                      genome_summary$genome[as.Date(paste0("20",genome_summary$collection_date)) == (first_day + day_entry)])
          
          day_existing_genome <- intersect(cluster_genomes,
                                           genome_summary$genome[as.Date(paste0("20",genome_summary$collection_date)) < (first_day + day_entry)])
          
          
          if(length(day_new_genome) > 0){
            
            new_distance <- sum(unlist(sapply(day_new_genome,
                                              function(new_genome_entry){
                                                # every new genome is compared to each of the existing genomes
                                                if(length(day_existing_genome) > 0){
                                                  
                                                  
                                                  new_genome_to_existing_distance <- unlist(sapply(day_existing_genome,
                                                                                                   function(existing_genome_entry){
                                                                                                     return(phylogentic_corrected_distance_matrix[new_genome_entry,
                                                                                                                                                  existing_genome_entry])
                                                                                                   }))
                                                  
                                                  return(
                                                    if(min(new_genome_to_existing_distance) < 0){
                                                      0
                                                    }else{
                                                      min(new_genome_to_existing_distance)
                                                    }
                                                    
                                                  )
                                                  
                                                }else{
                                                  return(0)
                                                }
                                              })))
          }else{
            new_distance <- 0
          }
          
          # add the new distance to accumulated_distance
          accumulated_distance <- accumulated_distance + new_distance
          
          distance_by_day <- rbind(distance_by_day,
                                   data.frame(
                                     cluster = cluster_summary[[cluster_index]]$cluster,
                                     invasive_status =cluster_summary[[cluster_index]]$invasive_status,
                                     date = names(persistent_period[persistent_period == day_entry]),
                                     day_post_initial_detection = day_entry,
                                     distance = new_distance,
                                     distance_category = "New"
                                   ),
                                   data.frame(
                                     cluster = cluster_summary[[cluster_index]]$cluster,
                                     invasive_status =cluster_summary[[cluster_index]]$invasive_status,
                                     date = names(persistent_period[persistent_period == day_entry]),
                                     day_post_initial_detection = day_entry,
                                     distance = accumulated_distance,
                                     distance_category = "Accumulated"
                                   ))
          
          
        }
        
      }
      
      ### The second is the floor plan showing the collection locations and Days post-initial detection of the cluster ----
      floorplan_df <- do.call(rbind,
                              lapply(persistent_period,
                                     function(day){
                                
                                # Genome collected at this day 
                                genome_this_day <- unique(genome_summary$genome[genome_summary$strain == cluster_summary[[cluster_index]]$strain & 
                                                               genome_summary$genome != "Environmental" &
                                                               as.Date(paste0("20",genome_summary$collection_date)) == as.Date(paste0("20",cluster_summary[[cluster_index]]$first_seen)) + day])
                                
                                  
                                if(length(genome_this_day) > 0){
                                  return(
                                    data.frame(
                                      genome = genome_this_day,
                                      day_post_inital_detection = day,
                                      bed = sapply(genome_this_day,
                                                   function(genome){
                                                     genome_summary$bed[genome_summary$genome == genome]
                                                   }),
                                      X = sapply(genome_this_day,
                                                 function(genome){
                                                   genome_bed <- genome_summary$bed[genome_summary$genome == genome]
                                                   genome_X <- bed_coord$X[bed_coord$bed == genome_bed]
                                                   if(length(genome_X) == 1){
                                                     genome_X
                                                   }else{
                                                     NA
                                                   }
                                                 }),
                                      Y = sapply(genome_this_day,
                                                 function(genome){
                                                   genome_bed <- genome_summary$bed[genome_summary$genome == genome]
                                                   genome_Y <- bed_coord$Y[bed_coord$bed == genome_bed]
                                                   if(length(genome_Y) == 1){
                                                     genome_Y
                                                   }else{
                                                     NA
                                                   }
                                                 })
                                    )
                                  )
                                }
                              }))
      
      distance_by_day$distance_category <- factor(distance_by_day$distance_category,
                                                  levels = c("New","Accumulated"))
      ### Visualization ----
      # y scale based on the max distance
      if (max(distance_by_day$distance) >= 500) {
        max_y <- ceiling(max(distance_by_day$distance) / 500) * 500
        plot_ratio <- 0.075
      } else if (max(distance_by_day$distance) >= 100) {
        max_y <- 500
        plot_ratio <- 0.075 * 3
      } else {
        max_y <- 100
        plot_ratio <- 0.075 * 15
      }
        
      
      
      # Line
      distance_plot <- ggplot(data = distance_by_day) +
        geom_line(aes(x = day_post_initial_detection, 
                      y = distance,
                      group = distance_category,
                      color = distance_category),
                  linewidth = 1.5,
                  alpha = 1) + 
        scale_x_continuous(limits = c(0, 1000)) +
        scale_y_continuous(limits = c(0, max_y)) +
        scale_color_manual(values = setNames(c("#DC8EB7","#8DA131"),
                                             c("New","Accumulated")),
                           name = "Category") + 
        labs(x = "Days Post Initial Detection",
             y = "Distance (m)") +
        theme(
          axis.line.x = element_line(linewidth = 1),
          axis.line.y = element_line(linewidth = 1),
          axis.text.x = element_text(face = "bold",
                                     angle = 0,
                                     size = 15),
          legend.key.size = unit(0.2,"cm"),
          legend.text = element_text(size = 12.5),
          legend.title = element_text(size = 12.5,
                                      face = "bold"),
          axis.text.y = element_text(face = "bold",
                                     size = 10),
          #axis.title.x = element_text(face = "bold", size = 20),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face = "bold", size = 15),
          plot.background = element_rect(fill = "transparent"),
          legend.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"),
          panel.border = element_rect(colour = "black",
                                      fill = NA,
                                      linewidth = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "inside",
          legend.position.inside = c(0.10, 0.6),
        ) + 
        coord_fixed(ratio = plot_ratio)
      
      
      # Get a bar indicating Days post-initial detection
      get_dpid_bar <- function(low_color = "#4575b4",
                               mid_color = "#ffeb3b",
                               high_color = "#d73027",
                               persistent_period = 0:1000){
        
        color_fn <- scales::col_numeric(
          palette = c(low_color, mid_color, high_color),
          domain = persistent_period
        )
        
        color_df <- data.frame(
          dpid = persistent_period,
          color = color_fn(persistent_period)
        )
        
        
        
        dpid_bar <- ggplot(color_df, aes(x = dpid,
                                         y = 1,
                                         fill = color)) +
          geom_tile() +
          labs(x = "Days Post Initial Detection") +
          scale_fill_identity() +
          theme_minimal() +
          scale_x_continuous(limits = c(0, 1000)) +
          theme(
            axis.title.x = element_text(face = "bold", size = 20),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank()
          ) +
          coord_fixed(ratio = 40)
        
        return(dpid_bar)
      }
      
      dpid_bar <- get_dpid_bar(low_color = "#4575b4",
                               mid_color = "#ffeb3b",
                               high_color = "#d73027",
                               persistent_period)
      
      # Floor plan
      floorplan <- map_plot +
        geom_point(data = floorplan_df,
                   aes(x=X,
                       y=Y,
                       colour = day_post_inital_detection),
                   alpha= 1,
                   size = 5,
                   shape=16)+
        scale_colour_gradient2(low = "#4575b4",
                               mid = "#ffeb3b",
                               high = "#d73027",
                               midpoint = median(persistent_period),
                               name = "Days Post Initial Detection") +
        geom_text_repel(data= unique(floorplan_df[,2:5]),
                        segment.linetype = "dotted",
                        segment.size = unit(0.75, "cm"),
                        min.segment.length = 0.5,
                        aes(x=X,
                            y=Y,
                            label = day_post_inital_detection,
                            colour = day_post_inital_detection),
                        size = 5.5,
                        box.padding = 0.35,
                        max.overlaps = 100,
                        show.legend=F,
                        direction = "both",
                        max.iter = 100000) + 
        theme(legend.text = element_text(size = 15),
              legend.title = element_text(size = 15),
              legend.key.size = unit(0.75,"cm"),
              legend.position = "none")
      
      # Combine 3 plots
      
      combined_plot <- distance_plot + dpid_bar +floorplan +
        plot_layout(ncol=1) +
        plot_annotation(
          title = paste0("Cluster",cluster_summary[[cluster_index]]$cluster),
          theme = theme(
            plot.title = element_text(
              hjust = 0.5,
              size = 25,
              face = "bold"
            )
          )
        )
      
      
      
      pdf(file= paste0("Cluster",
                       cluster_summary[[cluster_index]]$cluster,
                       "_Accumulated_Distance.pdf"),
          width=10,
          height=11.5)
      print(combined_plot)
      dev.off()
    }
  }
}

# Execute the function ----
floorplan_cluster_distance(output_dir,
                           summary_path,
                           bed_coor_path,
                           ml_tree_list,
                           floorplan_png_path,
                           pixel_scale)
