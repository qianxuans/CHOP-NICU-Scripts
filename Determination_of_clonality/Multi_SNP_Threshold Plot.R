#Visualize the clones/singletons/discrepancy/strain bootstrap support at each SNP threshold

# Libraries ----
library(ggplot2)
library(scales)
library(reshape2)
library(patchwork)
library(here)
# Input ----
output_path <- "./"
thresher_input <- readRDS("./input/thresher_input.RDS")
final_strains <- readRDS("./input/final_strains.RDS")
setwd(dir = output_path)
# Function to generate the plots ----
get_plots <- function(output_path,
                      determine_strains_input,
                      final_strains){
  #Universal variables/setting
  line_color <- c("#C72373",
                  "#86C2E8",
                  "#8F8578")
  names(line_color) <- c("Discrepancy",
                         "Clones",
                         "Singletons")
  # Iterate groups with no less than 2 genomes 
  for(i in seq_along(determine_strains_input)){
    group_id <- determine_strains_input[[i]][[1]]$HC_group
    group_input <- determine_strains_input[[i]]
    ## The main plot showing the numbers of clones/singletons/discrepancy at SNP thresholds ----
    ### Data frame for main plot 
    main_plot_df <- melt(do.call(rbind, 
                                 lapply(group_input, function(row) 
                                   data.frame(threshold = row$cutoff, 
                                              Discrepancy = row$discrepancy, 
                                              Singletons = row$after_correction_singletons, 
                                              Clones = row$after_correction_clones))), 
                         id.vars = "threshold", 
                         measure.vars = c("Discrepancy", "Singletons", "Clones"), 
                         variable.name = "category", 
                         value.name = "value")
    # The rectangle showing the plateau
    plateau_pos <- final_strains$plateaus$plateau[final_strains$plateaus$group == group_id]
    plateau_length <- unique(final_strains$plateaus$plateau_length)
    
    ### Main plot
    main_plot <- ggplot(main_plot_df,
                        aes(x = threshold,
                            y = value,
                            color = category)) +
      scale_color_manual(values = line_color,
                         name = "Number") +
      geom_line(alpha = 0.75,
                linewidth = 2) +
      # The plateau
      annotate("rect",
               xmin=plateau_pos,
               xmax=plateau_pos + plateau_length,
               ymin=-Inf,
               ymax=Inf,
               alpha=0.35,
               fill="#B3C16D") + 
      scale_x_continuous(name = "SNP Threshold",
                         breaks = seq(0, max(main_plot_df$threshold),50),
                         limits = c(min(main_plot_df$threshold),
                                    max(main_plot_df$threshold))) + 
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
        legend.key.size = unit(1.5, 'cm'),
        legend.position = "right",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 2)
      )
    
    
    ## The secondary plot showing the mean/median bootstrap supports for strains ----
    
    secondary_plot_df <- melt(do.call(rbind, 
                                      lapply(group_input, function(row) 
                                        data.frame(threshold = row$cutoff, 
                                                   Mean = row$mean_strain_bootstrap_support,
                                                   Median = row$median_strain_bootstrap_support))), 
                              id.vars = "threshold", 
                              measure.vars = c("Mean", "Median"), 
                              variable.name = "Bootstrap", 
                              value.name = "value")
    
    
    secondary_plot <- ggplot(secondary_plot_df,
                             aes(x = threshold,
                                 y = Bootstrap,
                                 fill = value)) + 
      geom_tile() + 
      scale_fill_viridis_c(name = "Value") + 
      scale_x_continuous(name = "SNP Threshold",
                         breaks = seq(0, max(secondary_plot_df$threshold),50),
                         limits = c(min(secondary_plot_df$threshold),
                                    max(secondary_plot_df$threshold))) + 
      theme(
        axis.title.y.left = element_text(colour = "black",
                                         size = 20,
                                         face = "bold"),
        axis.text.y.left = element_text(colour = "black",
                                        size = 20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        legend.title = element_text(size=20),
        legend.text = element_text(size=15),
        legend.key.size = unit(0.75, 'cm'),
        legend.position = "right",
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    linewidth = 2)
      ) + 
      coord_fixed(ratio=50)
    
    # Combine the main and secondary plots
    combined_plot <- secondary_plot + main_plot +
      plot_layout(ncol=1)
    
    # Export the combined plot
    pdf(file= paste0("Group",
                     group_id,
                     "_ThresholdVisualization.pdf"),
        width=15,
        height=10)
    print(combined_plot)
    dev.off()
  }
}

get_plots(output_path,
          determine_strains_input,
          final_strains)
