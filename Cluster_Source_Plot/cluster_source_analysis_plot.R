#make the plot where x axis indicates the date,  y axis indicates the strains
#the triangle indicates the admission date of the index patient
#the rectangle indicates the 95% HPD of the tMRCA
#the solid box indicates the mean of the tMRCA
#the dot indicates the first time the isolate of the strains was spotted
# Libraries ---- 
library(ggplot2)
library(dplyr)
library(swimplot)
library(dplyr)

# Input ----
#Note: The input data is not provided as it contains sensitive patient information that cannot be shared to ensure privacy and confidentiality.

# Helper functions ---
# Function to convert tMRCA to Unix timestamp (seconds since January 1, 1970) ----
convert_to_unix <- function(x) {
  date <- as.Date(paste0(floor(x), "-01-01")) + (x - floor(x)) * 365.25
  as.numeric(as.POSIXct(as.Date(date)))
}

# Function to analyze the transphylo output ----
analyze_transphylo <- function(transphylo_plot_dir,
                               transphylo_output_dir){
  
  setwd(dir = transphylo_plot_dir)
  
  # make the directories for 3 kinds of plots 
  sapply(c("transmission_tree", "colored_phylogeny", "incidentcases"), 
         dir.create, showWarnings = FALSE, recursive = TRUE)
  
  
  transphylo_output <- list.files(path = transphylo_output_dir,
                                  pattern = "*_TransPhylo.RDS",
                                  all.files = TRUE,
                                  full.names = TRUE,
                                  recursive = FALSE)
  
  
  transphylo_sum <- mclapply(transphylo_output,
                             function(output_idx){
                               
                               output_RDS <- readRDS(output_idx)
                               
                               # Get the range of cluster index case
                               
                               index_range <- output_RDS$sum_df
                               
                               # Check MCMC effective size
                               mcmc_check <- cbind(
                                 cluster = output_RDS$sum_df$cluster,
                                 output_RDS$mcmc_check
                               )
                               
                               
                               # Export the trees and inferred sampled and unsampled cases
                               
                               cluster_name <- gsub("cluster","Cluster",output_RDS$sum_df$cluster)
                               
                               pdf(file = file.path("./transmission_tree/",paste0(cluster_name,"_transmission_tree.pdf")),width = 15,height = 10)
                               plot(output_RDS$transmission_tree)
                               title(main = paste0(cluster_name, " Transmission Tree"))
                               dev.off()
                               
                               
                               pdf(file = file.path("./colored_phylogeny/",paste0(cluster_name,"_colored_phylogeny.pdf")),width = 15,height = 10)
                               plot(output_RDS$colored_phylo_tree)
                               title(main = paste0(cluster_name, " Colored Phylogeny"))
                               dev.off()
                               
                               
                               pdf(file = file.path("./incidentcases/",paste0(cluster_name,"_IncidentCases.pdf")),width = 15,height = 10)
                               cluster_ic <- getIncidentCases(output_RDS$raw_res_post_burnin,show.plot = TRUE)
                               title(main = paste0(cluster_name, " Sampled and Unsampled Cases"))
                               dev.off()
                               
                               # analysis of the Sampled and Unsampled Cases
                               totalcases <- sum(cluster_ic$sampledCases) + sum(cluster_ic$unsampCases)
                               sampledcases_pct <- 100 * cluster_ic$sampledCases/totalcases
                               unsampledcases_pct <- 100 * cluster_ic$unsampCases/totalcases
                               
                               return(
                                 list(
                                   index = index_range,
                                   mcmc_es = mcmc_check,
                                   sampledCases = list(
                                     raw = list(
                                       sampled_pct = sampledcases_pct,
                                       unsampled_pct = unsampledcases_pct
                                     ),
                                     stats_df = data.frame(
                                       mean_sampled_pct = mean(sampledcases_pct),
                                       sd_sampled_pct = sd(sampledcases_pct),
                                       mean_unsampled_pct = mean(unsampledcases_pct),
                                       sd_unsampled_pct = sd(unsampledcases_pct)
                                     )
                                   )
                                 )
                               )
                             },
                             mc.cores = 3)
  
  saveRDS(transphylo_sum,
          "transphylo_sum.RDS")
  
  return(transphylo_sum)
}

# Function to get the plot showing the probability of unsampled cases ----
get_unsampledcases_plot <- function(output_dir,
                                    transphylo_output_dir,
                                    summary_path){
  
  # Set working directory
  setwd(dir = output_dir)
  # Read Summary RDS 
  cluster_summary <- readRDS(summary_path)[["clusters"]]
  
  # If the sum of TransPhylo exists, skip this step
  if(file.exists(file.path("~/lab/cluster/NICU/TransPhylo/plots/transphylo_sum.RDS"))){
    transphylo_sum <- readRDS("~/lab/cluster/NICU/TransPhylo/plots/transphylo_sum.RDS")
  }else{
    transphylo_sum <- analyze_transphylo(transphylo_plot_dir = "~/lab/cluster/NICU/TransPhylo/plots",
                                         transphylo_output_dir = transphylo_output_dir)
  }
  
  # Get the data frame for plot
  
  plot_df <- do.call(rbind,
                     mclapply(transphylo_sum,
                              function(sum_entry){
                                
                                cluster_id <- gsub("cluster","",sum_entry$index$cluster)
                                
                                
                                cluster_sum_entry <- sapply(cluster_summary,
                                                            function(sum_entry){
                                                              if(sum_entry$cluster == cluster_id){
                                                                return(sum_entry)
                                                              }
                                                            }) %>%
                                  compact()
                                
                                
                                
                                
                                
                                return(
                                  cbind(
                                    cluster_id = cluster_id,
                                    invasive_status = cluster_sum_entry[[1]]$invasive_status,
                                    persistence = cluster_sum_entry[[1]]$persistence,
                                    sum_entry$sampledCases$stats_df
                                  )
                                )
                              },
                              mc.cores = detectCores()-4))
  
  # Stats
  # Colonizing Clusters n = 30 
  # Invasive Clusters n = 18
  # Colonizing Cluster Unsampled Cases Probabilities
  # Mean 3.76
  # Median 3.31
  # IQR [2.52,4.71]
  # Invasive Clusters Unsampled Cases Probabilities
  # Mean 4.64
  # Median 4.58
  # IQR [2.76,6.09]
  
  # Make the plots
  unsampledcases_plot <- ggplot(plot_df,
                                aes(x = invasive_status,
                                    y = mean_unsampled_pct,
                                    color = persistence)) + 
    geom_boxplot(alpha = 0.5,
                 size = 1.5,
                 outliers = FALSE) +
    geom_point(size = 5,
               alpha = 0.75,
               show.legend = FALSE,
               position = position_jitterdodge(jitter.width = 0.5,
                                               dodge.width = 0.5)) +
    geom_signif(comparisons = list(c(unique(plot_df$invasive_status))),
                map_signif_level = TRUE,
                test = "wilcox.test",
                test.args = list(exact = FALSE,
                                 correct = TRUE),
                size = 1.5,
                vjust = 0.5,
                textsize = 10,
                color = "black",
                show.legend = FALSE) + 
    scale_color_viridis_c(name = "Persistence (Days)") +
    labs(x = "Invasive Status",
         y = "% Unsampled Cases Probablity") +
    theme(axis.text.x = element_text(face = "bold",
                                     size = 15),
          axis.text.y = element_text(face = "bold",
                                     size = 15),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face = "bold",
                                      size = 17.5),
          plot.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.position = "right",
          axis.line = element_line(colour = "black",
                                   linewidth = 1),
          panel.border = element_rect(colour = "black",
                                      fill = NA,
                                      linewidth = 1)
    )
  
  # error bar code not used
  #geom_linerange(aes(ymin = mean_unsampled_pct - sd_unsampled_pct, 
  #                   ymax = mean_unsampled_pct + sd_unsampled_pct,
  #                   color = persistence), 
  #               linewidth = 0.5,
  #               alpha = 0.5,
  #               position = position_jitterdodge(jitter.width = 0.5,
  #                                               dodge.width = 0.5)) +
  
  
  pdf(file="UnsampledCasesPlot.pdf",
      width=7,
      height=5)
  print(unsampledcases_plot)
  dev.off()
  
}


# Function to get transmission source plot -----
get_transm_source_plot <- function(output_dir,
                                   transphylo_output_dir,
                                   admit_date_path,
                                   patientID_path,
                                   tmrca_df_path,
                                   summary_path){
  # Set working directory
  setwd(dir = output_dir)
  # Read Summary RDS 
  cluster_summary <- readRDS(summary_path)[["clusters"]]
  
  tmrca_df <- readRDS(tmrca_df_path)
  
  tmrca_df <- tmrca_df %>%
    mutate(across(starts_with("tmrca"), convert_to_unix))
  
  
  # If the sum of TransPhylo exists, skip this step
  if(file.exists(file.path("~/lab/cluster/NICU/TransPhylo/plots/transphylo_sum.RDS"))){
    transphylo_sum <- readRDS("~/lab/cluster/NICU/TransPhylo/plots/transphylo_sum.RDS")
  }else{
    transphylo_sum <- analyze_transphylo(transphylo_plot_dir = "~/lab/cluster/NICU/TransPhylo/plots",
                                         transphylo_output_dir = transphylo_output_dir)
  }
  
  # Get the index case duration
  cluster_idx_df <- do.call(rbind, lapply(transphylo_sum, function(sum_entry) {
    sum_entry$index$cluster <- gsub("cluster", "", sum_entry$index$cluster)
    sum_entry$index
  }))
  
  cluster_idx_df <- cluster_idx_df %>%
    mutate(across(starts_with("transm"), convert_to_unix))
  
  # Read admit dates 
  admit_date <- readRDS(admit_date_path)
  # Patient ID and MRNs
  patientID <- readRDS(patientID_path)
  #make the data frame as the input of the plot
  plot_df <- do.call(rbind,
                     lapply(tmrca_df$cluster,
                            function(cluster_id){
                              # Find the cluster
                              cluster <- cluster_summary[[which(sapply(cluster_summary,function(cluster) cluster$cluster == cluster_id))]]
                              # Find the admission date of the index patient
                              patients_admit_dates  <- setNames(as.Date(admit_date$date[admit_date$MRN %in% patientID$MRN[patientID$patientID %in% cluster$patient]]),
                                                                patientID$MRN[patientID$patientID %in% cluster$patient])
                              
                              data.frame(
                                cluster = cluster$cluster,
                                admission_date = min(as.numeric(as.POSIXct(as.Date(patients_admit_dates)))),
                                detection_date = as.numeric(as.POSIXct(as.Date(paste0("20",cluster$first_seen)))),
                                tmrca_mean = tmrca_df$tmrca_mean[tmrca_df$cluster == cluster$cluster],
                                tmrca_upper =  tmrca_df$tmrca_upper[tmrca_df$cluster == cluster$cluster],
                                tmrca_lower =  tmrca_df$tmrca_lower[tmrca_df$cluster == cluster$cluster],
                                transm_idx_medoid = cluster_idx_df$transm_index_mean[cluster_idx_df$cluster == cluster$cluster],
                                transm_idx_upper = cluster_idx_df$transm_index_upper[cluster_idx_df$cluster == cluster$cluster],
                                transm_idx_lower = cluster_idx_df$transm_index_lower[cluster_idx_df$cluster == cluster$cluster]
                              )
                            }))
  
  
  
  plot_df$cluster <- factor(plot_df$cluster,
                            # order by the lower limit of transmission origin date
                            levels = plot_df$cluster[sort(plot_df$transm_idx_mean,index.return = TRUE)$ix])
  
  xaxis_limit <- c(min(unlist(plot_df[,2:ncol(plot_df)])),
                   max(unlist(plot_df[,2:ncol(plot_df)])))
  
  plot_df$plot_y <- sapply(plot_df$cluster,
                           function(cluster_idx){
                             length(plot_df$cluster) - which(levels(plot_df$cluster) == cluster_idx)
                           })
  
  
  transm_phylo_plot <- ggplot(plot_df) +
    geom_hline(aes(yintercept = plot_y), 
               linetype = "solid", 
               color = "lightgrey", 
               alpha = 0.75) +
    # TMRCA 95%HPD
    geom_rect(aes(xmin = tmrca_lower,
                  xmax = tmrca_upper, 
                  ymin = plot_y - 0.4,
                  ymax = plot_y),
              fill = "#6FB3E1",
              color = "transparent",
              alpha = 0.75) +
    # TMRCA mean 
    geom_segment(aes(x = tmrca_mean,
                     xend = tmrca_mean,
                     y = plot_y - 0.4,
                     yend = plot_y),
                 color = "#305483",
                 linewidth = 1.5) +
    # Transmission 95%HPD
    geom_rect(aes(xmin = transm_idx_lower,
                  xmax = transm_idx_upper, 
                  ymin = plot_y,
                  ymax = plot_y + 0.4),
              fill = "#CBD69B",
              color = "transparent",
              alpha = 0.75) +
    # Transmission Index mean
    geom_segment(aes(x = transm_idx_medoid,
                     xend = transm_idx_medoid,
                     y = plot_y,
                     yend = plot_y + 0.4),
                 color = "#829825",
                 linewidth = 1.5) +
    # Index Patient Admission Date
    geom_point(aes(x = admission_date,
                   y = plot_y),
               shape = 17,
               size = 2.5,
               color = "#C72374",
               alpha = 0.75) +
    # Strain WGS Detection Date
    geom_point(aes(x = detection_date,
                   y = plot_y),
               shape = 19,
               size = 2.5, 
               color = "#C72374",
               alpha = 0.75) +
    ggplot2::scale_x_continuous(name = "Month",
                                limits = xaxis_limit,
                                labels = function(x) format(as.Date(as.POSIXct(x, origin = "1970-01-01"), tz = "UTC"), "%b %Y"),
                                n.breaks = 6) + 
    # Manual adjustment here 
    scale_x_break(breaks = c(1585699200,
                             1604188800),
                  scales = 7.5,
                  space = 0.2,
                  expand = FALSE) +
    scale_y_continuous(name = "Cluster",
                       # Use plot_y positions as break points
                       breaks = plot_df$plot_y,
                       # But display cluster values as labels
                       labels = plot_df$cluster,      
                       expand = c(0,0)) + 
    theme(axis.text.x = element_text(size = 12.5),
          axis.line.x = element_line(),
          axis.text.y = element_text(size = 12.5),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 15,
                                      face = "bold"),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_line(),
          plot.background = element_blank(),
          panel.background = element_rect(fill = "transparent"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_blank(),
          legend.title = element_blank(),
          legend.key = element_blank(),
          axis.text.x.top = element_blank(),
          axis.ticks.x.top = element_blank(),
          axis.line.x.top = element_blank()
    )
  
  pdf(file="tmrca_transmission_origin_index_admission_v3.pdf",
      width=6,
      height=10)
  print(transm_phylo_plot)
  dev.off()
}

# Execute the functions
get_transm_source_plot()
