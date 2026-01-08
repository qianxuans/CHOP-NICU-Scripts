# Persistence plot: 
#x axis indicates the date
#y axis indicates the clusters
#circles indicates the genomes
#colors of the circles show the invasive status 
# Libraries ----
library(ggplot2)
library(dplyr)
library(ggrepel)
# Input ----
persistence_plot_input <- readRDS("input/persistence_plot_input.RDS")
# Function to generate plots 
persistence_plot <- function(output_dir,
                             summary_path,
                             environmental_metadata_path){
  setwd(dir = output_dir)
  cluster_summary <- persistence_plot_input$cluster_summary
  genomes_summary <- persistence_plot_input$genomes_summary
  environmental_metadata <-  persistence_plot_input$environmental_metadata
  # Get the data frame as the input for the plot 
  all_plot_df <- do.call(rbind,
                     lapply(cluster_summary,
                            function(cluster){
                              
                              data.frame(
                                date = sapply(cluster$genomes,
                                              function(genome) paste0("20",genomes_summary$collection_date[genomes_summary$genome == genome])),
                                cluster = cluster$cluster,
                                status = sapply(cluster$genomes,
                                                function(genome) if(grepl("marc.",genome)) "Invasive" else if(genomes_summary$patientID[genomes_summary$genome == genome] == "Environmental") "Environmental" else "Colonizing"),
                                patient = sapply(cluster$genomes,
                                                 function(genome){
                                                   if(genomes_summary$patientID[genomes_summary$genome == genome] == "Environmental"){
                                                     
                                                     site <- environmental_metadata$site[environmental_metadata$sampleID == strsplit(genome,
                                                                                                                                     split = "\\.")[[1]][2]]
                                                     paste0(toupper(substr(site, 1, 1)), tolower(substr(site, 2, nchar(site))))
                                                     
                                                   }else{
                                                     genomes_summary$patientID[genomes_summary$genome == genome]
                                                   }
                                                   
                                                 })
                              )
                     }))
  # Get the data frame as the input for label
  all_label_df <- unique(all_plot_df)
  
  invasive_plot_df <- do.call(rbind,
                              lapply(cluster_summary,
                                     function(cluster){
                                       if(cluster$invasive_status == "Invasive"){
                                         data.frame(
                                           date = sapply(cluster$genomes,
                                                         function(genome) paste0("20",genomes_summary$collection_date[genomes_summary$genome == genome])),
                                           cluster = cluster$cluster,
                                           status = sapply(cluster$genomes,
                                                           function(genome) if(grepl("marc.",genome)) "Invasive" else if(genomes_summary$patientID[genomes_summary$genome == genome] == "Environmental") "Environmental" else "Colonizing"),
                                           patient = sapply(cluster$genomes,
                                                            function(genome){
                                                              if(genomes_summary$patientID[genomes_summary$genome == genome] == "Environmental"){
                                                                
                                                                site <- environmental_metadata$site[environmental_metadata$sampleID == strsplit(genome,
                                                                                                                                                split = "\\.")[[1]][2]]
                                                                paste0(toupper(substr(site, 1, 1)), tolower(substr(site, 2, nchar(site))))
                                                                
                                                              }else{
                                                                genomes_summary$patientID[genomes_summary$genome == genome]
                                                              }
                                                              
                                                            })
                                         )
                                       }
                                     }))
  invasive_label_df <- unique(invasive_plot_df)
    
  
  #Color used in the plot
  
  label_color <- c("#FFAFAF",
                   "#DAEAF8",
                   "#C7C0BB")
  names(label_color) <- c("Invasive",
                          "Colonizing",
                          "Environmental")
  
  label_border_color <- c("#CD6868",
                          "#25537D",
                          "#9B938E")
  names(label_border_color) <- c("Invasive",
                                 "Colonizing",
                                 "Environmental")
  
  # Change the class 
  all_plot_df$date <- as.Date(all_plot_df$date)
  invasive_plot_df$date <- as.Date(invasive_plot_df$date)
  all_label_df$date <- as.Date(all_label_df$date)
  invasive_label_df$date <- as.Date(invasive_label_df$date)
  all_plot_df$cluster <- factor(all_plot_df$cluster,
                                levels =  sort(unique(as.integer(all_plot_df$cluster)),decreasing = TRUE))
  invasive_plot_df$cluster <- factor(invasive_plot_df$cluster,
                                levels =  sort(unique(as.integer(invasive_plot_df$cluster)),decreasing = TRUE))
  invasive_plot_df$status <- factor(invasive_plot_df$status, 
                                    levels = c("Colonizing", "Invasive", "Environmental"))
  all_plot_df$status <- factor(all_plot_df$status,
                                    levels = c("Colonizing", "Invasive", "Environmental"))
  all_label_df$status <- factor(all_label_df$status,
                               levels = c("Colonizing", "Invasive", "Environmental"))
  invasive_label_df$status <- factor(invasive_label_df$status,
                               levels = c("Colonizing", "Invasive", "Environmental"))
  all_label_df$cluster <- factor(all_label_df$cluster,
                                levels =  sort(unique(as.integer(all_label_df$cluster)),decreasing = TRUE))
  invasive_label_df$cluster <- factor(invasive_label_df$cluster,
                                     levels =  sort(unique(as.integer(invasive_label_df$cluster)),decreasing = TRUE))
  
  # invasive cluster plot 
  
  invasive_persistence_plot <- ggplot(invasive_plot_df,
         aes(x = date,
             y = cluster)) +
    geom_line(linewidth = 1,
              aes(group = cluster),
              color = "black") +
    geom_point(size = 15,
               aes(color = status,
                   shape = status),
               alpha = 0.75) +
    scale_shape_manual(values = setNames(c(16,16,18),
                                         c("Colonizing","Invasive","Environmental")),
                       name = "Isolate Source",
                       guide = "none") + 
    scale_color_manual(values = setNames(c("#011F5B","#990000","#5C5F61"),
                                         c("Colonizing","Invasive","Environmental")),
                       name = "Isolate Source",
                       guide = guide_legend(override.aes = list(size = 20, 
                                                                shape = c(16, 16, 18)))) +
    scale_x_date(date_breaks = "6 month", 
                date_labels = "%b %Y") + 
    geom_label_repel(data = invasive_label_df,
                     aes(x = date,
                         y = cluster,
                         label = patient,
                         fill = status,
                         color = status),
                     direction = "both",
                     label.size = 0.7,
                     nudge_y = 0.175,
                     show.legend = FALSE,
                     segment.alpha = 0.35,
                     max.overlaps = 1000) +
    scale_fill_manual(values = label_color) +
    labs(x = "Date",
         y = "Cluster") +
    theme(
      axis.line.x = element_line(linewidth = 2),
      axis.ticks.y = element_blank(),
      axis.ticks.length = unit(0.5,"cm"),
      axis.text.x = element_text(face = "bold",
                                 angle = 45,
                                 hjust = 1,
                                 size = 40),
      axis.text.y = element_text(face = "bold",
                                 size = 40),
      axis.title.x = element_text(face = "bold", size = 50),
      axis.title.y = element_text(face = "bold", size = 50),
      plot.background = element_rect(fill = "transparent"),
      panel.background = element_rect(fill = "transparent"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = "transparent",
                                fill = "transparent"),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 25),
    )
  
  pdf(file="Invasive_clusters_persistence.pdf",
      width=17.5,
      height=17.5)
  print(invasive_persistence_plot)
  dev.off()
  
  
  #For all clusters
  
  all_persistence_plot <- ggplot(all_plot_df,
                                      aes(x = date,
                                          y = cluster)) +
    geom_line(linewidth = 1,
              aes(group = cluster),
              color = "black") +
    geom_point(size = 15,
               aes(color = status,
                   shape = status),
               alpha = 0.75) +
    scale_shape_manual(values = setNames(c(16,16,18),
                                         c("Colonizing","Invasive","Environmental")),
                       name = "Isolate Source",
                       guide = "none") + 
    scale_color_manual(values = setNames(c("#011F5B","#990000","#5C5F61"),
                                         c("Colonizing","Invasive","Environmental")),
                       name = "Isolate Source",
                       guide = guide_legend(override.aes = list(size = 20, 
                                                                shape = c(16, 16, 18)))) +
    scale_x_date(date_breaks = "6 month", 
                 date_labels = "%b %Y") + 
    geom_label_repel(data = all_label_df,
                     aes(x = date,
                         y = cluster,
                         label = patient,
                         fill = status,
                         color = status),
                     direction = "both",
                     label.size = 0.7,
                     nudge_y = 0.175,
                     show.legend = FALSE,
                     segment.alpha = 0.35,
                     max.overlaps = 1000) +
    scale_fill_manual(values = label_color) +
    labs(x = "Date",
         y = "Cluster") +
    theme(
      axis.line.x = element_line(linewidth = 2),
      axis.ticks.y = element_blank(),
      axis.ticks.length = unit(0.5,"cm"),
      axis.text.x = element_text(face = "bold",
                                 angle = 45,
                                 hjust = 1,
                                 size = 40),
      axis.text.y = element_text(face = "bold",
                                 size = 40),
      axis.title.x = element_text(face = "bold", size = 50),
      axis.title.y = element_text(face = "bold", size = 50),
      plot.background = element_rect(fill = "transparent"),
      panel.background = element_rect(fill = "transparent"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = "transparent",
                                fill = "transparent"),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 25),
    )
  
  pdf(file="All_clusters_persistence.pdf",
      width=17.5,
      height=50)
  print(all_persistence_plot)
  dev.off()
  
  # For only clusters with environmental samples 
  
  environmental_plot_df <- all_plot_df[all_plot_df$cluster %in% unique(all_plot_df$cluster[grepl("Environmental",all_plot_df$status)]),]
  environmental_label_df <- all_label_df[all_label_df$cluster %in% unique(all_plot_df$cluster[grepl("Environmental",all_plot_df$status)]),]
  # Change the class 
  environmental_plot_df$date <- as.Date(environmental_plot_df$date)
  environmental_label_df$date <- as.Date(environmental_label_df$date)
  environmental_plot_df$cluster <- factor(environmental_plot_df$cluster,
                                levels =  sort(unique(as.integer(environmental_plot_df$cluster)),decreasing = TRUE))
  environmental_plot_df$status <- factor(environmental_plot_df$status, 
                                    levels = c("Colonizing", "Invasive", "Environmental"))
  
  
  env_persistence_plot <- ggplot(environmental_plot_df,
                                 aes(x = date,
                                     y = cluster)) +
    geom_line(linewidth = 1,
              aes(group = cluster),
              color = "black") +
    geom_point(size = 20,
               aes(color = status,
                   shape = status),
               alpha = 0.75) +
    scale_shape_manual(values = setNames(c(16,16,18),
                                         c("Colonizing","Invasive","Environmental")),
                       name = "Isolate Source",
                       guide = "none") + 
    scale_color_manual(values = setNames(c("#011F5B","#990000","#5C5F61"),
                                         c("Colonizing","Invasive","Environmental")),
                       name = "Isolate Source",
                       guide = guide_legend(override.aes = list(size = 25, 
                                                                shape = c(16, 16, 18)))) +
    scale_x_date(date_breaks = "6 month", 
                 date_labels = "%b %Y") + 
    geom_label_repel(data = environmental_label_df,
                     aes(x = date,
                         y = cluster,
                         label = patient,
                         fill = status,
                         color = status),
                     direction = "both",
                     label.size = 0.75,
                     size = 5,
                     nudge_y = 0.175,
                     show.legend = FALSE,
                     segment.alpha = 0.35,
                     max.overlaps = 1000) +
    scale_fill_manual(values = label_color) +
    labs(x = "Date",
         y = "Cluster") +
    theme(
      axis.line.x = element_line(linewidth = 2),
      axis.ticks.y = element_blank(),
      axis.ticks.length = unit(0.5,"cm"),
      axis.text.x = element_text(face = "bold",
                                 angle = 45,
                                 hjust = 1,
                                 size = 30),
      axis.text.y = element_text(face = "bold",
                                 size = 30),
      axis.title.x = element_text(face = "bold", size = 35),
      axis.title.y = element_text(face = "bold", size = 35),
      plot.background = element_rect(fill = "transparent"),
      panel.background = element_rect(fill = "transparent"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = "transparent",
                                fill = "transparent"),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 25),
    )
  
  pdf(file="Environmental_clusters_persistence.pdf",
      width=12.5,
      height=10)
  print(env_persistence_plot)
  dev.off()
}