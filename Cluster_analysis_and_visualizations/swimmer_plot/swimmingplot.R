# Libraries ----
library(openxlsx2)
library(dplyr)
library(ggplot2)
library(swimplot)
library(colorspace)
library(here)
# Paths and Input ----
output_path <- "./"
setwd(output_path)
#bed_coor_path NOT PROVIDED
#summary_path NOT PROVIDED
#pt_timeline_path NOT PROVIDED
#pt_tx_team_path NOT PROVIDED
# Plots indicating the patient movement ----
movement_treatment_plot <- function(summary_path,
                                    bed_coor_path,
                                    pt_timeline_path,
                                    pt_tx_team){
  ## Universal variables / Input ----
  ### Patient Timeline 
  pt_timeline <- readRDS(pt_timeline_path)
  ### Patient Treatment Team
  pt_tx_team <- readRDS(pt_tx_team_path)
  ### Treatment team colors 
  treatment_team_colors <- setNames(
    c("#89CFF0", "#98FF98", "#FFB347", "#E0B0FF", "#FF6B6B", "#FFE066", "#F4A6B7"),
    c("Blue","Green","Orange","Purple","Red","Yellow","KOPH")
  )
  ### Section colors 
  section_color <- setNames(
    c("#FFD8A9","#BDBDDC", "#B6FFB6","#6FCCB2","#D9C78E","#BDBDBD","#87CEEB","#FFB6C1"),
    c("A","B","C","D","F","G","H","KOPH")
  )
  ### Bed X Y coordinates
  bed_coor <- read.table(bed_coor_path,
                         header = TRUE)
  
  ### Cluster summary
  cluster_summary <- readRDS(summary_path)[["clusters"]]
  
  ### Genome summary
  genome_summary <- readRDS(summary_path)[["genomes"]]
  
  for(i in seq_along(cluster_summary)){
    
    cluster_patients <- cluster_summary[[i]]$patient
    
    # Get the data frame for patient movement (Patient Movement Plot)
    pt_treatment_team <- do.call(rbind,
                              lapply(cluster_patients,
                                     function(patient_entry){
                                       
                                       pt_timeline_patient <- pt_timeline[pt_timeline$PatientID == patient_entry,] %>%
                                         arrange(as.Date(BED_ENTER_DATE))
                                       
                                       if(nrow(pt_timeline_patient)>0){
                                         return(
                                           data.frame(
                                             patient = patient_entry,
                                             start_date = as.numeric(as.POSIXct(as.Date(pt_timeline_patient$BED_ENTER_DATE))),
                                             end_date = as.numeric(as.POSIXct(as.Date(pt_timeline_patient$BED_EXIT_DATE))),
                                             # Find the treatment team for each bed
                                             treatment_team = sapply(seq_len(nrow(pt_timeline_patient)),
                                                                     function(bed_entry){
                                                                       as.character(ifelse(nrow(pt_tx_team %>%
                                                                                                  filter(PatientID == patient_entry,
                                                                                                         between(as.Date(specimen_taken_date),
                                                                                                                 as.Date(pt_timeline_patient$BED_ENTER_DATE[bed_entry]),
                                                                                                                 as.Date(pt_timeline_patient$BED_EXIT_DATE[bed_entry])))) > 0,
                                                                                           {pt_tx_team %>%
                                                                                               filter(PatientID == patient_entry,
                                                                                                      between(as.Date(specimen_taken_date),
                                                                                                              as.Date(pt_timeline_patient$BED_ENTER_DATE[bed_entry]),
                                                                                                              as.Date(pt_timeline_patient$BED_EXIT_DATE[bed_entry]))) %>%
                                                                                               pull(treatment_team) %>%
                                                                                               unique() %>%
                                                                                               gsub("NICU ", "", .)},
                                                                                           NA))
                                                                     })
                                           )
                                         )
                                       }
                                       
                                     }))
    
    # Get the data frame for treatment team (Treatment Team Plot) 
    pt_movement <- do.call(rbind,
                           lapply(cluster_patients,
                                  function(patient_entry){
                                    pt_timeline_patient <- pt_timeline[pt_timeline$PatientID == patient_entry,] %>%
                                      arrange(as.Date(BED_ENTER_DATE))
                                    if(nrow(pt_timeline_patient)>0){
                                      return(
                                        data.frame(
                                          patient = patient_entry,
                                          start_date = as.numeric(as.POSIXct(as.Date(pt_timeline_patient$BED_ENTER_DATE))),
                                          end_date = as.numeric(as.POSIXct(as.Date(pt_timeline_patient$BED_EXIT_DATE))),
                                          section = sapply(pt_timeline_patient$BED_NAME,function(bed) if(length(bed_coor$area[bed_coor$bed == bed]) == 1) bed_coor$area[bed_coor$bed == bed] else "Non-NICU")
                                        )
                                      )
                                    }
                                  }))
    
    # Get the data frame for cluster collection dates
    collection_dates <- do.call(rbind,
                                lapply(cluster_patients,
                                       function(patient_entry){
                                         
                                         genomes_patient <- genome_summary$genome[genome_summary$patientID == patient_entry]
                                         #the genomes must be within the cluster
                                         genomes_patient <- genomes_patient[genomes_patient %in% cluster_summary[[i]]$genomes]
                                         
                                         
                                         return(
                                           data.frame(
                                             patient = patient_entry,
                                             collection_date = as.numeric(as.POSIXct(as.Date(paste0("20",genome_summary$collection_date[genome_summary$genome %in% genomes_patient])))),
                                             status = sapply(genome_summary$genome[genome_summary$genome %in% genomes_patient],
                                                             function(genome)
                                                               if(grepl("marc.",genome)) "Invasive" else "Colonizing")
                                           )
                                         )
                                         
                                       }))
    # Get the data frame for patients' admission and discharge date
    admission_discharge <- do.call(rbind,
                                   lapply(cluster_patients,
                                          function(patient_entry){
                                            admission_discharge_patient <- unique(pt_timeline[pt_timeline$PatientID == patient_entry,] %>%
                                                                                    select(PatientID,
                                                                                           HOSPITAL_ADMIT_DATE,
                                                                                           HOSPITAL_DISCHARGE_DATE))
                                            
                                            return(data.frame(
                                              patient = patient_entry,
                                              date = c(as.numeric(as.POSIXct(as.Date(admission_discharge_patient$HOSPITAL_ADMIT_DATE))),
                                                       as.numeric(as.POSIXct(as.Date(admission_discharge_patient$HOSPITAL_DISCHARGE_DATE)))),
                                              event = c(rep("Hospital Admission",length(admission_discharge_patient$HOSPITAL_ADMIT_DATE)),
                                                        rep("Hospital Discharge",length(admission_discharge_patient$HOSPITAL_DISCHARGE_DATE)))
                                            )
                                            )
                                          }))
    
    #Use the 3 data frames above to make the patient-movement plot
    # Id Order
    
    id_order <- admission_discharge %>% 
      filter(event == "Hospital Admission") %>% 
      group_by(patient) %>%
      slice(1) %>%
      arrange(desc(date)) %>%
      pull(patient) %>%
      factor()
    
    
    ### Movement plot
    movement_plot <- swimmer_plot(df = pt_movement,
                                id='patient',
                                start = 'start_date',
                                id_order = id_order,
                                end='end_date',
                                name_fill ='section') +
      #show the legend for Department
      guides(fill = guide_legend(order=1))+
      #isolate collection dates and isolate source
      swimmer_points(df_points=collection_dates,
                     id='patient',
                     time='collection_date',
                     size=5,
                     alpha = 0.75,
                     fill = "white",
                     name_col = 'status') +
      ggplot2::scale_color_manual(name="Isolate Source",values=c('#011F5B',
                                                                 '#990000'))+
      ggplot2::scale_fill_manual(name="NICU Section",
                                 values=section_color,
                                 na.value=NA,
                                 na.translate=FALSE) +
      guides(color = guide_legend(override.aes = list(size = 15,
                                                      linetype = 0,
                                                      alpha = 1,
                                                      fill = NA),
                                  order = 2)) +
      #patient admission,discharge/death dates
      swimmer_points(df_points=admission_discharge[complete.cases(admission_discharge),],
                     id='patient',
                     time='date',
                     size=3,
                     alpha = 0.5,
                     stroke = 3.5,
                     fill = "white",
                     color = "#B22163",
                     name_shape = 'event') +
      ggplot2::scale_shape_manual(name="Patients Hospitalization",
                                  values=c("Hospital Admission"=17, "Hospital Discharge"=15)) +
      guides(shape = guide_legend(override.aes = list(size = 3,
                                                      stroke = 3.5,
                                                      alpha = 0.25),
                                  order = 3))+
      ggplot2::scale_y_continuous(name = "Date",
                                  labels = function(x) format(as.Date(as.POSIXct(x, origin = "1970-01-01"), tz = "UTC"), "%b %Y")) +
      ggplot2::scale_x_discrete(name = "Patient ID") +
      coord_flip(ylim = c(min(c(na.omit(admission_discharge$date),
                                na.omit(collection_dates$collection_date),
                                na.omit(pt_movement$start_date))),
                          max(c(na.omit(admission_discharge$date),
                                na.omit(collection_dates$collection_date),
                                na.omit(pt_movement$end_date))))) +
      theme(axis.text.x = element_text(face = "bold",
                                       angle = 45,
                                       hjust = 1,
                                       size = 20),
            axis.text.y = element_text(face = "bold",
                                       size = 20),
            axis.title.x = element_text(face = "bold", size = 35),
            axis.title.y = element_text(face = "bold", size = 35),
            plot.background = element_blank(),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line.x = element_line(),
            axis.ticks.x = element_line(),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 20))
    
    
    ### Treatment team plot
    treatment_plot <- swimmer_plot(df = pt_treatment_team,
                                  id='patient',
                                  start = 'start_date',
                                  id_order = id_order,
                                  end='end_date',
                                  name_fill ='treatment_team') +
      #show the legend for Department
      guides(fill = guide_legend(order=1))+
      #isolate collection dates and isolate source
      swimmer_points(df_points=collection_dates,
                     id='patient',
                     time='collection_date',
                     size=5,
                     alpha = 0.75,
                     fill = "white",
                     name_col = 'status') +
      ggplot2::scale_color_manual(name="Isolate Source",values=c('#011F5B',
                                                                 '#990000'))+
      ggplot2::scale_fill_manual(name="NICU Treatment Team",
                                 values=treatment_team_colors,
                                 na.value=NA,
                                 na.translate=FALSE) +
      guides(color = guide_legend(override.aes = list(size = 15,
                                                      linetype = 0,
                                                      alpha = 1,
                                                      fill = NA),
                                  order = 2)) +
      #patient admission,discharge/death dates
      swimmer_points(df_points=admission_discharge[complete.cases(admission_discharge),],
                     id='patient',
                     time='date',
                     size=3,
                     alpha = 0.5,
                     stroke = 3.5,
                     fill = "white",
                     color = "#B22163",
                     name_shape = 'event') +
      ggplot2::scale_shape_manual(name="Patients Hospitalization",
                                  values=c("Hospital Admission"=17, "Hospital Discharge"=15)) +
      guides(shape = guide_legend(override.aes = list(size = 3,
                                                      stroke = 3.5,
                                                      alpha = 0.25),
                                  order = 3))+
      ggplot2::scale_y_continuous(name = "Date",
                                  labels = function(x) format(as.Date(as.POSIXct(x, origin = "1970-01-01"), tz = "UTC"), "%b %Y")) +
      ggplot2::scale_x_discrete(name = "Patient ID") +
      coord_flip(ylim = c(min(c(na.omit(admission_discharge$date),
                                na.omit(collection_dates$collection_date),
                                na.omit(pt_movement$start_date))),
                          max(c(na.omit(admission_discharge$date),
                                na.omit(collection_dates$collection_date),
                                na.omit(pt_movement$end_date))))) +
      theme(axis.text.x = element_text(face = "bold",
                                       angle = 45,
                                       hjust = 1,
                                       size = 20),
            axis.text.y = element_text(face = "bold",
                                       size = 20),
            axis.title.x = element_text(face = "bold", size = 35),
            axis.title.y = element_text(face = "bold", size = 35),
            plot.background = element_blank(),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line.x = element_line(),
            axis.ticks.x = element_line(),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 20))
    
    
    
    pdf(file= paste0("Cluster",
                    cluster_summary[[i]]$cluster,
                    "_movement.pdf"),
        width=12.5,
        height=10)
    print(movement_plot)
    dev.off()
    
    
    pdf(file= paste0("Cluster",
                     cluster_summary[[i]]$cluster,
                     "_treatment.pdf"),
        width=12.5,
        height=10)
    print(treatment_plot)
    dev.off()
    
  }
}

movement_treatment_plot(summary_path,
                        bed_coor_path,
                        pt_timeline_path,
                        pt_tx_team)
