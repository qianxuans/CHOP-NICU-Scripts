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
#output_dir <- #
#tmrca_df_path <- #
#summary_path <- #
#patientID_path <- #
#admit_date_path <- #

# Function to get the tmrca admission plot ----
tmrca_admission_plot <- function(output_dir,
                                 admit_date_path,
                                 patientID_path,
                                 tmrca_df_path,
                                 summary_path){
  # Set working directory
  setwd(dir = output_dir)
  # Read RDS 
  # Only keep the clusters with no less than 4 genomes
  cluster_summary <- readRDS(summary_path)[["clusters"]]
  
  tmrca_df <- readRDS(tmrca_df_path)
  #convert tMRCA to Unix timestamp (seconds since January 1, 1970)
  convert_to_unix <- function(x) {
    date <- as.Date(paste0(floor(x), "-01-01")) + (x - floor(x)) * 365.25
    as.numeric(as.POSIXct(as.Date(date)))
  }
  
  tmrca_df <- tmrca_df %>%
    mutate(across(starts_with("tmrca"), convert_to_unix))
  
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
                                group = "tMRCA"
                              )
                            }))
  
  #swim plot of all cluster####
  
  tmrca_admission_plot <- swimmer_plot(df = plot_df,
                                       id = "cluster",
                                       id_order='tmrca_mean',
                                       start = "tmrca_lower",
                                       end = "tmrca_upper",
                                       name_fill ='group',
                                       alpha = 0.5) +
    ggplot2::scale_fill_manual(name=NA,
                               values = "#6FB3E1",
                               na.value=NA,
                               na.translate=FALSE) +
    guides(fill = "none") +
    swimmer_points(df_points = plot_df,
                   id = 'cluster',
                   time = 'tmrca_mean',
                   size= 3.5,
                   shape = 15,
                   color = "#305483",
                   alpha = 0.75) + 
    swimmer_points(df_points = plot_df,
                   id = 'cluster',
                   time = 'detection_date',
                   size= 3.5,
                   shape = 19,
                   color = "#C72374",
                   alpha = 0.75) + 
    swimmer_points(df_points = plot_df,
                   id = 'cluster',
                   time = 'admission_date',
                   size= 3.5,
                   shape = 17,
                   color = "#C72374",
                   alpha = 0.75) + 
    ggplot2::scale_x_discrete(name = "Cluster") + 
    ggplot2::scale_y_continuous(name = "Month",
                                labels = function(x) format(as.Date(as.POSIXct(x, origin = "1970-01-01"), tz = "UTC"), "%b %Y"),
                                n.breaks = 6)+ 
    coord_flip(ylim = c(min(plot_df$tmrca_lower),
                        max(plot_df$tmrca_upper))) + 
    theme(axis.text.x = element_text(face = "bold",
                                     angle = 0,
                                     vjust = 0.5,
                                     size = 20),
          axis.text.y = element_text(face = "bold",
                                     size = 15),
          axis.title.x = element_text(face = "bold", size = 20),
          axis.title.y = element_text(face = "bold", size = 20),
          plot.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(),
          axis.ticks.x = element_line(),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(1, 5, 1, 1), "cm"))
  
  pdf(file=paste0("clusters_origin_analysis_plot.pdf"),
      width=10,
      height=10)
  print(tmrca_admission_plot)
  dev.off()
  
}

tmrca_admission_plot(output_dir,
                     admit_date_path,
                     patientID_path,
                     tmrca_df_path,
                     summary_path)