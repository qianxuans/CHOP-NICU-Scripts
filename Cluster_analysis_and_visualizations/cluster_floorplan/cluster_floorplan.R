# Floor plan showing epidemiological links in the clusters
# Libraries ----
library(ggplot2)
library(gganimate)
library(ggmap)
library(colorspace)
library(dplyr)
library(ggrepel)
library(openxlsx2)
library(magick)
library(here)
# Input ----
output_dir <- "./"
#summary_path NOT PROVIDED
#bed_coor_path NOT PROVIDED
#floorplan_png_path NOT PROVIDED
#pt_timeline_path NOT PROVIDED
#environmental_metadata_path NOT PROVIDED

# Get the transmission floor plan of clusters -----
cluster_floorplan <- function(summary_path,
                              floorplan_png_path,
                              bed_coor_path,
                              pt_timeline_path,
                              environmental_metadata_path){
  setwd(dir = output_dir)
  ## Universal variables / Input ----
  # floorplan png
  map_plot <- image_read(floorplan_png_path)
  map_plot <- image_ggplot(map_plot)
  ### Bed X Y coordinates
  bed_coor <- read.table(bed_coor_path,
                         header = TRUE)
  ### Patient Timeline
  pt_timeline <- readRDS(pt_timeline_path)
  ### Cluster summary
  cluster_summary <- readRDS(summary_path)[["clusters"]]
  ### Environmental sample metadata if applicable
  
  for(i in seq_along(cluster_summary)){
    ### Get the data frame as the input of the plots
    pt_movement_df <- do.call(rbind,
                              lapply(cluster_summary[[i]]$patient,
                                     function(patient_entry){
                                       pt_timeline_tmp <- pt_timeline[pt_timeline$PatientID == patient_entry,]
                                       if(nrow(pt_timeline_tmp) > 0){
                                         # if the patient_timeline_tmp$BED_NAME[k] == patient_timeline_tmp$BED_NAME[k+1]
                                         # there was actually no changing the bed
                                         # fix this:
                                         pt_timeline_tmp <- pt_timeline_tmp %>%
                                           mutate(
                                             to_merge = lead(BED_NAME) == BED_NAME & 
                                               BED_EXIT_DATE == lead(BED_ENTER_DATE)
                                           ) %>%
                                           mutate(
                                             BED_EXIT_DATE = if_else(to_merge, lead(BED_EXIT_DATE), BED_EXIT_DATE)
                                           ) %>%
                                           filter(is.na(lag(to_merge)) | !lag(to_merge)) %>%
                                           select(-to_merge) %>%
                                           select(PatientID,
                                                  BED_NAME,
                                                  BED_ENTER_DATE,
                                                  BED_EXIT_DATE)
                                         
                                         pt_movement_df_tmp <- data.frame(
                                           X = sapply(pt_timeline_tmp$BED_NAME,function(bed) bed_coor$X[bed_coor$bed == bed]),
                                           Xend = lead(sapply(pt_timeline_tmp$BED_NAME,function(bed) bed_coor$X[bed_coor$bed == bed])),
                                           Y = sapply(pt_timeline_tmp$BED_NAME,function(bed) bed_coor$Y[bed_coor$bed == bed]),
                                           Yend = lead(sapply(pt_timeline_tmp$BED_NAME,function(bed) bed_coor$Y[bed_coor$bed == bed])),
                                           bed = pt_timeline_tmp$BED_NAME,
                                           patient = patient_entry,
                                           date = pt_timeline_tmp$BED_ENTER_DATE,
                                           month = format(as.Date(pt_timeline_tmp$BED_ENTER_DATE),"%b %Y")
                                         )
                                         
                                         # Make the Xend and Yend NA if the destination is the same although at different times
                                         pt_movement_df_tmp$Xend[pt_movement_df_tmp$Xend == pt_movement_df_tmp$X & pt_movement_df_tmp$Yend == pt_movement_df_tmp$Y] <- NA
                                         pt_movement_df_tmp$Yend[pt_movement_df_tmp$Xend == pt_movement_df_tmp$X & pt_movement_df_tmp$Yend == pt_movement_df_tmp$Y] <- NA
                                         
                                         return(pt_movement_df_tmp)
                                       }
                                       
                                     }))
    #re-order the patient
    #from small number to large number
    
    pt_movement_df$patient <- factor(pt_movement_df$patient,
                                              levels = sort(unique(as.integer(pt_movement_df$patient))))
    
    
    #If this is run for the first time 
    #the color scheme of the patients in the cluster is created and saved 
    #otherwise load the patient color scheme
    if(file.exists(paste0("./clusters_color_schemes/cluster",cluster_summary[[i]]$cluster))){
      patient_colors_df <- read.table(paste0("./clusters_color_schemes/cluster",
                                             cluster_summary[[i]]$cluster),
                                      header = FALSE,
                                      comment.char = "")
      
      patient_colors <- setNames(as.character(patient_colors_df$V2),
                                 as.character(patient_colors_df$V1))
    }else{
      #this color scheme generation is randomly everytime 
      #so for the 1st time, I save the color schemes in the directory "clusters_color_schemes"
      if(!dir.exists("./clusters_color_schemes")){
        system("mkdir clusters_color_schemes")
      }
      patient_colors <- setNames(qualitative_hcl(length(levels(pt_movement_df$patient)), palette = "Dark 3"),
                                 levels(pt_movement_df$patient))
      
      write.table(patient_colors,
                  paste0("./clusters_color_schemes/cluster",
                         cluster_summary[[i]]$cluster),
                  quote = FALSE,
                  col.names = FALSE)
    }
    
    # Add environmental sample if applicable
    if(cluster_summary[[i]]$environmental_sample){
      env_metadata <- read_xlsx(environmental_metadata_path)
      genomes_sum <- readRDS(summary_path)[["genomes"]]
      genomes_sum <- genomes_sum[genomes_sum$MRN == "Environmental",]
      env_entry <- sapply(cluster_summary[[i]]$genomes[cluster_summary[[i]]$genomes %in% genomes_sum$genome],
                            function(genome) strsplit(genome,split = "\\.")[[1]][2])
      env_df <- unique(do.call(rbind,
                        lapply(env_entry,
                               function(entry){
                                 
                                   beds <- strsplit(env_metadata$locationID[env_metadata$sampleID == entry],
                                                    split = ", ")[[1]]
                                   
                                   site <- env_metadata$site[env_metadata$sampleID == entry]
                                   site <- paste0(toupper(substr(site, 1, 1)), tolower(substr(site, 2, nchar(site))))
                                     
                                   do.call(rbind,
                                           lapply(beds,function(bed_entry){
                                             data.frame(
                                               X = bed_coor$X[bed_coor$bed == bed_entry],
                                               Y = bed_coor$Y[bed_coor$bed == bed_entry],
                                               label = paste0(site,
                                                              ": ",
                                                              env_metadata$collection_date[env_metadata$sampleID == entry])
                                             )}))
                               }
                        )))
    }
    
    #Static floorplan 
    
    static_floorplan <- map_plot +
      geom_point(data = pt_movement_df,
                 aes(x=X,
                     y=Y,
                     colour = patient),
                 alpha= 0.35,
                 size = 5,
                 shape=16)+
      guides(fill = guide_legend(order=1))+
      geom_curve(data= pt_movement_df[complete.cases(pt_movement_df),],
                 aes(x = X,
                     y = Y,
                     xend = Xend,
                     yend = Yend,
                     colour = patient),
                 arrow = arrow(length = unit(0.25, "cm")),
                 alpha = 1,
                 curvature = -0.5) +
      geom_text_repel(data= pt_movement_df,
                      segment.linetype = "dotted",
                      segment.size = unit(0.75, "cm"),
                      min.segment.length = 0.5,
                      aes(x=X,
                          y=Y,
                          label = month,
                          colour = patient),
                      size = 3.5,
                      box.padding = 0.35,
                      max.overlaps = 100,
                      show.legend=F,
                      direction = "both",
                      max.iter = 100000)+
      # Add the environmental sample info if applicable
      {
        if (cluster_summary[[i]]$environmental_sample) {
          list(
            geom_point(data = env_df,
                       aes(x = X,
                           y = Y),
                       color = "#9B938E",
                       alpha = 0.35,
                       size = 10,
                       shape = 18),
            geom_text_repel(data = env_df,
                            segment.linetype = "dotted",
                            segment.size = unit(0.75, "cm"),
                            min.segment.length = 0.5,
                            aes(x = X,
                                y = Y,
                                label = label),
                            colour = "#5C5F61",
                            size = 3.5,
                            box.padding = 0.35,
                            max.overlaps = 100,
                            show.legend = FALSE,
                            direction = "both",
                            max.iter = 100000)
          )
        } else {
          NULL
        }
      }+
      guides(color = guide_legend(order = 2))+
      scale_color_manual(values = patient_colors,
                         name = "Patient Movement") +
      theme(legend.text = element_text(size = 15),
            legend.title = element_text(size = 15),
            legend.key.size = unit(0.75,"cm"))
    
    pdf(file= paste0("Cluster",
                     cluster_summary[[i]]$cluster,
                     "_Static.pdf"),
        width=10,
        height=10)
    print(static_floorplan)
    dev.off()
    
    #Animated floorplan
    #A different version of data frame for the animated floorplan
    
    pt_movement_animated_df <- do.call(rbind,
                              lapply(cluster_summary[[i]]$patient,
                                     function(patient_entry){
                                       
                                       pt_timeline_tmp <- pt_timeline[pt_timeline$PatientID== patient_entry,]
                                       # if the patient_timeline_tmp$BED_NAME[k] == patient_timeline_tmp$BED_NAME[k+1]
                                       # there was actually no changing the bed
                                       # fix this:
                                       pt_timeline_tmp <- pt_timeline_tmp %>%
                                         mutate(
                                           to_merge = lead(BED_NAME) == BED_NAME & 
                                             BED_EXIT_DATE == lead(BED_ENTER_DATE)
                                         ) %>%
                                         mutate(
                                           BED_EXIT_DATE = if_else(to_merge, lead(BED_EXIT_DATE), BED_EXIT_DATE)
                                         ) %>%
                                         filter(is.na(lag(to_merge)) | !lag(to_merge)) %>%
                                         select(-to_merge) %>%
                                         select(PatientID,
                                                BED_NAME,
                                                BED_ENTER_DATE,
                                                BED_EXIT_DATE)
                                      
                                       
                                       return(do.call(rbind,
                                               lapply(seq_len(nrow(pt_timeline_tmp)),
                                                      function(row_entry){
                                                        
                                                        seq_dates <- if(!is.na(pt_timeline_tmp$BED_EXIT_DATE[row_entry])){
                                                          as.character(seq(as.Date(pt_timeline_tmp$BED_ENTER_DATE[row_entry]),
                                                                           as.Date(pt_timeline_tmp$BED_EXIT_DATE[row_entry]),
                                                                           by = "1 day"))
                                                        }else{
                                                          as.character(pt_timeline_tmp$BED_ENTER_DATE[row_entry])
                                                        }
                                                        
                                                        data.frame(
                                                          X = bed_coor$X[bed_coor$bed == pt_timeline_tmp$BED_NAME[row_entry]],
                                                          Y = bed_coor$Y[bed_coor$bed == pt_timeline_tmp$BED_NAME[row_entry]],
                                                          bed = pt_timeline_tmp$BED_NAME[row_entry],
                                                          patient = patient_entry,
                                                          date = seq_dates)
                                                         
                                                      })))
                                     }))
    
    pt_movement_animated_df$date <- as.Date(pt_movement_animated_df$date)
    
    
    animated_floorplan1 <- map_plot + 
      geom_point(data = pt_movement_animated_df,
                 aes(x = X,
                     y = Y,
                     colour = patient,
                     group = patient),
                 alpha = 0.4,
                 size = 17.5,
                 shape=16) + 
      geom_text(data = pt_movement_animated_df,
                aes(x=X,
                    y=Y,
                    label = patient,
                    colour = patient),
                size = 7.5,
                fontface='bold') +
      scale_color_manual(values = patient_colors) +
      guides(color = "none", shape = "none") + 
      scale_alpha_identity() + 
      transition_time(date) + 
      labs(title = paste0("Cluster",
                          cluster_summary[[i]]$cluster,
                          '\nDate: {frame_time}')) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(size = 25,
                                  hjust = 0.5)
      ) 
    
    animated_floorplan1 <- animate(animated_floorplan1,
                                              fps = 4,
                                              width=800,
                                              height=800,
                                              renderer = gifski_renderer(loop = FALSE))
    
    anim_save(animation = animated_floorplan1,
              filename = paste0("Cluster",
                                cluster_summary[[i]]$cluster,
                                "_animated_1.gif"))
    
    #Animated Version 2 
    
    pt_movement_df$date <- as.Date(pt_movement_df$date)
    
    curve_df <- pt_movement_df %>%
      group_by(patient) %>%
      arrange(date) %>%
      mutate(date = lead(date)) %>%
      filter(!is.na(date))
    
    animated_floorplan2 <- map_plot +
      geom_point(data = pt_movement_df,
                 aes(x=X,
                     y=Y,
                     colour = patient),
                 alpha= 0.4,
                 size = 17.5,
                 shape=16) +
      {if(nrow(curve_df) > 0){
        geom_curve(data= curve_df,
                   aes(x = X,
                       y = Y,
                       xend = Xend,
                       yend = Yend,
                       colour = patient),
                   arrow = arrow(length = unit(0.25, "cm")),
                   alpha = 1,
                   curvature = -0.5)
      }} +
      geom_text(data = pt_movement_df,
                aes(x=X,
                    y=Y,
                    label = patient,
                    colour = patient),
                size = 7.5,
                fontface='bold') +
      labs(title = 'Date: {closest_state}') + 
      guides(fill = guide_legend(order=1),
             color = "none")+
      scale_color_manual(values = patient_colors,
                         name = "Patient Movement") +
      theme(plot.title = element_text(size = 25,
                                      hjust = 0.5)) + 
      transition_states(date,
                        wrap = FALSE)+
      shadow_mark()
    
    animated_floorplan2 <- animate(animated_floorplan2,
                                                   fps = 4,
                                                   nframes = 100,
                                                   width=800,
                                                   height=800,
                                                   renderer = gifski_renderer(loop = FALSE))
    
    anim_save(animation = animated_floorplan2,
              filename = paste0("Cluster",
                                cluster_summary[[i]]$cluster,
                                "_animated_2.gif"))
    
  }
  
}


cluster_floorplan(summary_path,
                  floorplan_png_path,
                  bed_coor_path,
                  pt_timeline_path,
                  environmental_metadata_path)
