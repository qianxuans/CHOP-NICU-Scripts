# Libraries ----
library(ggplot2)
library(dplyr)
#for rounded rectangle
library(ggchicklet)
#geom_circle
library(ggforce)
library(openxlsx2)
library(shadowtext)
#other
library(here)
# Input ----
summary_path <- "./input/cluster_summary.RDS"
icon_path <- "./input/cluster_icon.RDS"

# Function to draw the cluster plot ----
clusterplot <- function(output_path,
                        summary_path,
                        icon_path){
  ## Universal parameters for plots ----
  icon <- readRDS(icon_path)
  
  #Size for export 
  pdf_export_size <- data.frame(num_circle = 2:100) %>%
    mutate(
      width = 5,
      height = case_when(
        num_circle <= 5 ~ 3.75,
        TRUE ~ ((ceiling(num_circle/5)*2.5*4 + 19 + 15)/44) * 3.75
      )
    )
  
  # Colors for invasive status
  status_color <- c("Colonizing" = "#011F5B",
                    "Invasive" = "#990000",
                    "Invasive/Colonizing" = "#F4A261")
  
  ## Read the summary and plot for each cluster ----
  all_summary <- readRDS(summary_path)
  ## iterate and make the plot for each cluster ----
  for(i in seq_along(all_summary$clusters)){
    
   
    if(all_summary$clusters[[i]]$environmental_sample){
      # If there is environmental sample,
      # The environmental section is also a member in the plot
      # The status is the name of the string in character plot_members
      plot_members <-  c(sort(gsub("NICU ","",unique(all_summary$clusters[[i]]$environmental_sample_section)),
                              decreasing = TRUE),
                         sort(as.integer(all_summary$clusters[[i]]$patient),
                              decreasing = TRUE))
    }else{
      plot_members <- sort(as.integer(all_summary$clusters[[i]]$patient),
                           decreasing = TRUE)
    }
    
    plot_members_num <- length(plot_members)
    
    # The data frame as ggplot input
    # There are 2 situations when making the plot:
    # 1: no more than 5 patients(1 line in the plot)
    # 2: more than 5 patients (more than 1 line in the plot)
    ### The plot_df for ggplot input ----
    if(plot_members_num <= 5){
      #situation1:
      #less than 5 patients in the cluster
      #there will always be only 1 line in the plot
      #so y of every patient is always 0.4 
      #the x of each patient will be set proportionally
      if(plot_members_num == 5){
        #If there are exactly 5 patients in the cluster
        #only 1 line
        #the x of each patient will be set proportionally(from n=1 to n=5, x = 0.1,0.3,0.5,0.7,and 0.9 respectively)
        #the y will be 0.5
        plot_df <- do.call(rbind,
                           lapply(seq_along(plot_members),
                                  function(entry){
                                    
                                    member_status <- if(grepl("^[A-Za-z]+$",plot_members[entry])){
                                      # If letter indicating NICU section, Environmental sample
                                      "Environmental"
                                    }else{
                                      # If number, the patient ID
                                      # All the genomes of this patient of this cluster
                                      
                                      member_genomes <- all_summary$genomes$genome[all_summary$genomes$patientID == plot_members[entry]] 
                                      #the member genomes have to be within the cluster
                                      member_genomes <- member_genomes[member_genomes %in% all_summary$clusters[[i]]$genomes]
                                      
                                      if(any(grepl("marc.",member_genomes))){
                                        
                                        if(all(grepl("marc.",member_genomes))){
                                          
                                          "Invasive"
                                        }else{
                                          
                                          "Invasive/Colonizing"
                                        }
                                        
                                      }else{
                                        "Colonizing"
                                      }
                                    }
                                    
                                    
                                    data.frame(
                                      ID = plot_members[entry],
                                      status = member_status,
                                      x = 1-((2*entry-1)/(2*5)),
                                      y = 0.5
                                    )
                                  }))
        
      }else{
        #less than 5 patients in the cluster
        #there will always be only 1 line in the plot
        #so y of every patient is always 0.4 
        #the x of each patient will be set proportionally
        
        plot_df <- do.call(rbind,
                           lapply(seq_along(plot_members),
                                  function(entry){
                                    
                                    member_status <- if(grepl("^[A-Za-z]+$",plot_members[entry])){
                                      # If letter indicating NICU section, Environmental sample
                                      "Environmental"
                                    }else{
                                      # If number, the patient ID
                                      # All the genomes of this patient of this cluster
                                      
                                      member_genomes <- all_summary$genomes$genome[all_summary$genomes$patientID == plot_members[entry]] 
                                      #the member genomes have to be within the cluster
                                      member_genomes <- member_genomes[member_genomes %in% all_summary$clusters[[i]]$genomes]
                                      
                                      if(any(grepl("marc.",member_genomes))){
                                        
                                        if(all(grepl("marc.",member_genomes))){
                                          
                                          "Invasive"
                                        }else{
                                          
                                          "Invasive/Colonizing"
                                        }
                                        
                                      }else{
                                        "Colonizing"
                                      }
                                    }
                                    
                                    data.frame(
                                      ID = plot_members[entry],
                                      status = member_status,
                                      x = 1-((2*(entry %% 5)-1)/(2*(plot_members_num %% 5))),
                                      y = 0.5
                                    )
                                  }))
        
      }
    }else{
      #situation2:
      #more 5 patients in the cluster
      #there will be more than 1 line in the plot
      #the x of each patient will be set proportionally
      #the y will vary depending on which line n is at
      #ceiling(nrow(plot_df)/5) is the total number of lines 
      #ceiling(n/5) indicates which line "n" is at
      #n %% 5 == 0 indicates this is the last one in this line
      plot_df <- do.call(rbind,
                         lapply(seq_along(plot_members),
                                function(entry){
                                  
                                  entry_id <- plot_members[entry]
                                  
                                  member_status <- if(grepl("^[A-Za-z]+$",plot_members[entry])){
                                    # If letter indicating NICU section, Environmental sample
                                    "Environmental"
                                  }else{
                                    # If number, the patient ID
                                    # All the genomes of this patient of this cluster
                                    
                                    member_genomes <- all_summary$genomes$genome[all_summary$genomes$patientID == plot_members[entry]] 
                                    #the member genomes have to be within the cluster
                                    member_genomes <- member_genomes[member_genomes %in% all_summary$clusters[[i]]$genomes]
                                    
                                    if(any(grepl("marc.",member_genomes))){
                                      
                                      if(all(grepl("marc.",member_genomes))){
                                        
                                        "Invasive"
                                      }else{
                                        
                                        "Invasive/Colonizing"
                                      }
                                      
                                    }else{
                                      "Colonizing"
                                    }
                                  }
                                    
                                  if(ceiling(entry/5) < ceiling(plot_members_num/5)){
                                    #ceiling(n/5) < ceiling(nrow(plot_df)/5
                                    #this line IS NOT the last line in the plot
                                    if(entry %% 5 == 0){
                                      data.frame(
                                        ID = entry_id,
                                        status = member_status,
                                        x = 0.1,
                                        y = (2*ceiling(entry/5)-1)/(ceiling(plot_members_num/5)*2)
                                      )
                                    }else{
                                      data.frame(
                                        ID = entry_id,
                                        status = member_status,
                                        x = 1-((2*(entry %% 5)-1)/10),
                                        y = (2*ceiling(entry/5)-1)/(ceiling(plot_members_num/5)*2)
                                      )
                                    }
                                  }else if(ceiling(entry/5) == ceiling(plot_members_num/5) & plot_members_num %% 5 == 0){
                                    
                                    #the LAST line in the plot
                                    #and also if nrow(plot_df) %% 5 == 0, which means the patient number IS 10,15,20,25, etc. 
                                    if(entry %% 5 == 0){
                                      data.frame(
                                        ID = entry_id,
                                        status = member_status,
                                        x = 0.1,
                                        y = (2*ceiling(entry/5)-1)/(ceiling(plot_members_num/5)*2)
                                      )
                                      
                                    }else{
                                      data.frame(
                                        ID = entry_id,
                                        status = member_status,
                                        x = 1-((2*(entry %% 5)-1)/(2*5)),
                                        y = (2*ceiling(entry/5)-1)/(ceiling(plot_members_num/5)*2)
                                      )
                                    }
                                  }else if(ceiling(entry/5) == ceiling(plot_members_num/5) & plot_members_num %% 5 != 0){
                                    #the LAST line in the plot
                                    #and also if nrow(plot_df) %% 5 != 0, which means the patient number IS NOT 10,15,20,25, etc. 
                                    if(entry %% 5 == 0){
                                      
                                      data.frame(
                                        ID = entry_id,
                                        status = member_status,
                                        x = 0.1,
                                        y = (2*ceiling(entry/5)-1)/(ceiling(plot_members_num/5)*2)
                                      )
                                      
                                    }else{
                                      data.frame(
                                        ID = entry_id,
                                        status = member_status,
                                        x = 1-((2*(entry %% 5)-1)/(2*(plot_members_num %% 5))),
                                        y = (2*ceiling(entry/5)-1)/(ceiling(plot_members_num/5)*2)
                                      )
                                    }
                                  }
                                }))
    }
    
    
    #now use plot_df to make the plot using ggplot2
    #In the rectangle of the plot, I want each line to have at most 5 dots
    #otherwise it will be too crowded
    #Assume the length of the rectangle is L
    #the number of the circles is N, the r of the circle is R = 2.5
    #the position of the circle center at X place is:  L(2X-1/2N)
    #I set L(the length of the rectangle) based on the number of the circles: L = 2*N*R (N < 5)
    #if N >= 5, L =  10R = 25
    #H(height of the rectangle) = ceiling(N/5)*R*4
    H <- ceiling(plot_members_num/5)*2.5*4
    
    L <- 2*10*2.5
  
    ### generate the plot ----
    #icon of sections and treatment team
    if(any(grepl("Invasive",plot_df$status))){
      treatment_team_icon <- icon$treatment_team_invasive
      section_icon <- icon$section_invasive
    }else{
      treatment_team_icon <- icon$treatment_team_colonzing
      section_icon <- icon$section_colonizing
    }
    
    #ST of the strain 
    cluster_st <- if(any(grepl("Unassigned",all_summary$clusters[[i]]$ST))) "-" else unique(all_summary$clusters[i][[1]]$ST)[1]
    #CC of the cluster 
    cluster_cc <- if(any(grepl("Unassigned",all_summary$clusters[[i]]$CC))) "-" else unique(all_summary$clusters[i][[1]]$CC)[1]
    
    #data frame for treatment team and NICU section 
    
    #nicu section 
    
    section_colors <- setNames(
      c("#FFD8A9", "#BDBDDC", "#B6FFB6", "#6FCCB2", "#D9C78E", "#BDBDBD", "#87CEEB", "#FFB6C1"),
      c("A", "B", "C", "D", "F", "G", "H", "K")
    )
    
    sections <- sort(unique(c(all_summary$clusters[[i]]$patient_section,
                              all_summary$clusters[[i]]$environmental_sample_section)))
    
    sections_df <- data.frame(
      x = 4.5 + (50-4.5-2.5) * seq_along(sections)/(length(sections) + 1),
      y = ((7.5+3)+(7.5+3+7))/2,
      chr = gsub("NICU |OPH 5 NICU", "", sections),
      color = section_colors[gsub("NICU |OPH 5 NICU", "", sections)]
    )
    # treatment team 
    treatment_team_colors <- setNames(
      c("#89CFF0", "#98FF98", "#FFB347", "#E0B0FF", "#FF6B6B", "#FFE066", "#F4A6B7"),
      c("B","G","O","P","R","Y","K")
    )
    
    treatment_teams <- sapply(gsub("NICU ","",sort(all_summary$clusters[[i]]$treatment_team)),
                              function(team) strsplit(team,split = "")[[1]][1])
    
    treatment_team_df <- data.frame(
      x = 4.5 + (50-4.5-2.5) * seq_along(treatment_teams)/(length(treatment_teams) + 1),
      y = (2.5+9.5)/2,
      chr = treatment_teams,
      color = treatment_team_colors[treatment_teams]
    )
    
    cluster_plot <- ggplot() +
      ggchicklet:::geom_rrect(aes(xmin = 0,
                                  xmax = L,
                                  ymin = 0,
                                  ymax = H+19+15),
                              fill = ifelse(any(grepl("Invasive",plot_df$status)),
                                            "#FFAFAF",
                                            "#DAEAF8"),
                              color = ifelse(any(grepl("Invasive",plot_df$status)),
                                             "#CD6868",
                                             "#25537D"),
                              size = 3,
                              alpha = 1,
                              r = unit(0.1, 'npc')) +
      # upper bar showing clusterID, strainID, ST, CC and MRSA
      ggchicklet:::geom_rrect(aes(xmin = 2.5,
                                  xmax = L-2.5,
                                  ymin = H+19+2.5,
                                  ymax = H+19+15-2.5),
                              fill = ifelse(any(grepl("Invasive",plot_df$status)),
                                            "#FFEFEF",
                                            "#F2FAFF"),
                              color = "transparent",
                              size = 0,
                              alpha = 1,
                              r = unit((5/15), 'npc')) + 
      # 1st lower bar showing the treatment team
      ggchicklet:::geom_rrect(aes(xmin = 2.5,
                                  xmax = L-2.5,
                                  ymin = 2.5,
                                  ymax = 9.5),
                              fill = ifelse(any(grepl("Invasive",plot_df$status)),
                                            "#FFEFEF",
                                            "#F2FAFF"),
                              color =  "transparent",
                              size = 0,
                              alpha = 1,
                              r = unit(0.5, 'npc'))  + 
      # 2nd lower bar showing the NICU sections
      ggchicklet:::geom_rrect(aes(xmin = 2.5,
                                  xmax = L-2.5,
                                  ymin = 7.5+3,
                                  ymax = 7.5+3+7),
                              fill = ifelse(any(grepl("Invasive",plot_df$status)),
                                            "#FFEFEF",
                                            "#F2FAFF"),
                              color =  "transparent",
                              size = 0,
                              alpha = 1,
                              r = unit(0.5, 'npc'))  + 
      #Only patient are presented as dots
      geom_circle(data = plot_df[plot_df$status != "Environmental",],
                  aes(x0 = x*L,
                      y0 = 19+y*H,
                      r = 4,
                      fill = status),
                  color = "transparent",
                  linewidth = 1) + 
      #Environmental samples presented as diamonds (if applicable)
      {if(any(plot_df$status == "Environmental")) {
        geom_polygon(data = do.call(rbind,
                                    lapply(seq_len(nrow(plot_df[plot_df$status == "Environmental",])),
                                           function(entry){
                                             center_x <- plot_df[plot_df$status == "Environmental",][entry,3] * L
                                             center_y <- (plot_df[plot_df$status == "Environmental",][entry,4] * H) + 19
                                             return(data.frame(
                                               x = c(center_x,
                                                     center_x+5,
                                                     center_x,
                                                     center_x-5),
                                               y = c(center_y-5,
                                                     center_y,
                                                     center_y+5,
                                                     center_y),
                                               section = rep(plot_df[plot_df$status == "Environmental",][entry,1],4)
                                             ))
                                           })),
                     mapping=aes(x=x,
                                 y=y,
                                 group=section),
                     fill = "#C7C0BB",
                     color = "transparent")
      }} +
      scale_fill_manual(values = status_color) + 
      guides(fill = "none") + 
      #The patient ID inside the circles and NICU section of environmental sample within the diamonds (if applicable)
      geom_text(data = plot_df,
                aes(x = x*L,
                    y = 19+y*H,
                    label = ID),
                size = 9,
                color = "white",
                fontface = "bold") + 
      #Cluster ID
      geom_text(aes(x = L/2,
                    y = H+19+10,
                    label = paste0("Cluster",
                                   all_summary$clusters[i][[1]]$cluster)),
                size = 10.5,
                color = "black",
                fontface = "bold")+
      #ST/CC
      geom_text(aes(x = L/2,
                    y = H+19+2.5+3,
                    label = paste0(cluster_st,
                                   "/",
                                   cluster_cc)),
                size = 7.5,
                color = "black") +
      ggplot2::coord_fixed(ratio = 1,
                           xlim = c(0, L),
                           ylim = c(0, H+19+15)) +
      #MRSA or MSSA
      geom_shadowtext(aes(x = 0.925*L-2.5,
                          y = H+19+(15/2),
                          label = ifelse(all_summary$clusters[i][[1]]$MRSA == "MRSA",
                                         "R",
                                         "S")),
                      color = ifelse(all_summary$clusters[i][[1]]$MRSA == "MRSA",
                                     "#FF69B4",
                                     "#FFD700"),
                      size = 21.5,
                      fontface = "bold",
                      bg.colour = NA) + 
      # NICU sections 
      geom_shadowtext(data = sections_df,
                      aes(x = x,
                          y = y,
                          label = chr),
                      color = sections_df$color,
                      size = 15,
                      fontface = "bold",
                      bg.colour = NA) + 
      # Treatment team
      geom_shadowtext(data = treatment_team_df,
                      aes(x = x,
                          y = y,
                          label = chr),
                      color = treatment_team_df$color,
                      size = 15,
                      fontface = "bold",
                      bg.colour = NA) + 
      # NICU section svg icon
      annotation_custom(section_icon,
                        xmin = 2.75,
                        xmax = 8,
                        ymin = 11.5,
                        ymax = Inf) +
      # Treatment team svg icon
      annotation_custom(treatment_team_icon,
                        xmin = 3.5,
                        xmax = 8,
                        ymin = 3.5,
                        ymax = Inf) +
      
      theme_void()
    
    
    pdf(file = paste0("Cluster",all_summary$clusters[[i]]$cluster,".pdf"),
        width = as.numeric(pdf_export_size$width[pdf_export_size$num_circle == plot_members_num]),
        height = as.numeric(pdf_export_size$height[pdf_export_size$num_circle == plot_members_num]))
    print(cluster_plot)
    dev.off()
  }
}

  


clusterplot(output_path,summary_path,icon_path)
