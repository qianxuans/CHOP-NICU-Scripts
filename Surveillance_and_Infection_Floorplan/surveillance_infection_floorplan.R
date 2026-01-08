# Libraries ----
library(dplyr)
library(openxlsx2)
library(magick)
library(ggplot2)
library(gganimate)

# Input ----
# The input data is not provided as it contains sensitive patient information that cannot be shared to ensure privacy and confidentiality
# Function to visualize the carriage and infections ----

visualize_carriage_infection <- function(output_path,
                                         input_df,
                                         bed_coor_path,
                                         map_png_path){
  # Working directory
  setwd(dir = output_path)
  # Import the floorplan
  map_plot <- image_read(map_png_path)
  map_plot <- image_ggplot(map_plot)
  # Import the bed coordinates
  bed_coor <- read.table(bed_coor_path,
                         header = TRUE)
  
  input_df <- merge(input_df,
                    bed_coor,
                    by = "bed",
                    all.x = TRUE,
                    all.y = FALSE) %>%
    filter(!is.na(X))
  
  # Generate the animated plot
  animated_plot <- map_plot +
    geom_point(data = input_df,
               aes(x=X,
                   y=Y,
                   size=ifelse(status == "Invasive", 4, frequency),
                   colour=status,
                   shape=status,
                   stroke = ifelse(status == "Invasive", 4, 0)),
               alpha= 0.75) +
    scale_color_manual(values = c("Colonizing" = "#011F5B", "Invasive" = "#990000")) +
    scale_shape_manual(values = c("Colonizing" = 16, "Invasive" = 4),
                       name = "Invasive Status") +
    scale_size_continuous(breaks = c(1,5,10,15,20,25,30),
                          name = "Frequency")+
    transition_states(collection_date,
                      wrap = FALSE) +
    shadow_mark() +
    labs(title = "Date: {closest_state}") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(), 
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      plot.title = element_text(size = 25,
                                hjust = 0.5),
      legend.text = element_text(size = 17.5),
      legend.title = element_text(size = 22.5)
    ) + 
    guides(
      shape = "none",
      color = "none"
    )
  
  animated_plot <- animate(animated_plot,
                               width=2700,
                               height=3000,
                               units = "px",
                               res = 377,
                               renderer = gifski_renderer(loop = FALSE))
  
  anim_save(animation = animated_plot,
            filename = "animated_surveillance_infection_floorplan.gif")
  
  # Static Plot
  static_plot <- map_plot +
    geom_point(data = input_df[input_df$status == "Colonizing",],
               aes(x=X,
                   y=Y,
                   size=frequency),
               colour="#011F5B",
               alpha= 1) +
    scale_size_continuous(breaks = c(1,5,10,15,20,25,30),
                          name = "Frequency")+
    geom_point(data = input_df[input_df$status == "Invasive",],
               aes(x=X,
                   y=Y),
               colour="#990000",
               alpha= 0.75,
               shape = 4,
               size = 4,
               stroke = 4)+ 
    theme(legend.text = element_text(size = 25),
          legend.title = element_text(size = 25))
  pdf(file="static_surveillance_infection_floorplan.pdf",
      width=10,
      height=10)
  print(static_plot)
  dev.off()
   
}

visualize_carriage_infection(output_path,
                             input_df,
                             bed_coor_path,
                             map_png_path)

