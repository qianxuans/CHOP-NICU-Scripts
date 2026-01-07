# Use msm package for HMM analysis to estimate the probability of false negative
# Compare colonizing and invasive clusters
# Libraries ----
library(msm)
library(dplyr)
library(scales)
library(openxlsx2)
library(ggplot2)
library(ggsignif)
library(parallel)
library(ggforce)
# Input ----
cluster_msm_input <- get_cluster_msm_input(surveillance,
                                           cluster_sum,
                                           patient_id_mrn)
cluster_msm_input <- readRDS("input/cluster_msm_input.RDS")
# Perform msm ----
exec_cluster_msm <- function(cluster_msm_input,
                             definition,
                             ncores = 4,
                             msm_seed = 1205){
  
  # set the seed for reproducibility
  set.seed(msm_seed)
  
  if(definition == "culture"){
    
    msm_input <- do.call(rbind,
                         lapply(cluster_msm_input,
                                function(input_entry){
                                  input_entry$culture_input
                                }))
    
  }else if(definition == "phylogenomics"){
    
    msm_input <- do.call(rbind,
                         lapply(cluster_msm_input,
                                function(input_entry){
                                  input_entry$phylogenomics_input
                                }))
    
  }
  
  msm_input$invasive_status <- factor(msm_input$invasive_status,
                                      levels = unique(msm_input$invasive_status))
  
  msm <- msm(result ~ month,
             data = msm_input,
             subject = clusterID,
             cl = 0.95,
             misccovariates = ~ invasive_status, 
             qmatrix = rbind(
               c(0.7,0.3),
               c(0.3,0.7)
             ),
             # underlying states as rows, and observed states as columns
             ematrix = rbind(
               c(0.95,0.05),
               c(0.05,0.95)
             ),
             control = list(maxit=2000000000))
  
  # Misclassification probability
  # Estimated misclassification probability matrix.
  # The rows correspond to true states, and columns observed states.
  colonizing_msm_ematrix <- ematrix.msm(
    msm,
    covariates = list(invasive_status = "Colonizing"),
    ci = "bootstrap",
    cl = 0.95,
    B = 100,
    cores = ncores
  )
  
  invasive_msm_ematrix <- ematrix.msm(
    msm,
    covariates = list(invasive_status = "Invasive"),
    ci = "bootstrap",
    cl = 0.95,
    B = 100,
    cores = ncores
  )
  
  nocovar_msm_ematrix <- ematrix.msm(
    msm,
    covariates = 0,
    ci = "bootstrap",
    cl = 0.95,
    B = 100,
    cores = ncores
  )
  
  
  # helper function to get the data frame for ematrix
  get_plot_df <- function(msm_ematrix){
   return(
     data.frame(observed = c("Positive","Positive","Negative","Negative"),
                true = c("Positive","Negative","Positive","Negative"),
                # The rows correspond to true states, and columns observed states.
                probability = c(msm_ematrix$estimates[4]*100,
                                msm_ematrix$estimates[3]*100,
                                msm_ematrix$estimates[2]*100,
                                msm_ematrix$estimates[1]*100),
                low = c(msm_ematrix$L[4]*100,
                        msm_ematrix$L[3]*100,
                        msm_ematrix$L[2]*100,
                        msm_ematrix$L[1]*100),
                up = c(msm_ematrix$U[4]*100,
                       msm_ematrix$U[3]*100,
                       msm_ematrix$U[2]*100,
                       msm_ematrix$U[1]*100))
   ) 
  }
  
  # return the results
  list(
    full = list(
      msm = msm,
      colonizing_ematrix = colonizing_msm_ematrix,
      invasive_ematrix = invasive_msm_ematrix,
      nocovar_ematrix = nocovar_msm_ematrix
    ),
    colonizing_misclassification_df = get_plot_df(colonizing_msm_ematrix),
    invasive_misclassification_df = get_plot_df(invasive_msm_ematrix),
    nocor_misclassification_df = get_plot_df(nocovar_msm_ematrix)
  )
}

# Get msm result
# culture definition 
culture_msm_results <- exec_cluster_msm(cluster_msm_input = cluster_msm_input,
                                        definition = "culture",
                                        ncores = 4,
                                        msm_seed = 1205)
saveRDS(culture_msm_results,
        file = "culture_msm_result.RDS")
# phylogenomics definition 
phylogenomics_msm_results <- exec_cluster_msm(cluster_msm_input = cluster_msm_input,
                                              definition = "phylogenomics",
                                              ncores = 4,
                                              msm_seed = 1205)
saveRDS(phylogenomics_msm_results,
        file = "phylogenomics_msm_result.RDS")


# Visualize misclassification Probability ----
visualize_misclassification <- function(cluster_msm_results,
                                        definition){
  
  
  # Sort the data frame for ggplot2
  
  cluster_emission_df <- rbind(
    cbind(
      invasive_status = "Colonizing",
      result = paste0("Observed ",cluster_msm_results$colonizing_misclassification_df$observed,
                      " True ",cluster_msm_results$colonizing_misclassification_df$true),
      cluster_msm_results$colonizing_misclassification_df
    ),
    cbind(
      invasive_status = "Invasive",
      result = paste0("Observed ",cluster_msm_results$invasive_misclassification_df$observed,
                      " True ",cluster_msm_results$invasive_misclassification_df$true),
      cluster_msm_results$invasive_misclassification_df
    )
  )
  
  
  
  
  
  # misclassification
  misclassification_df <- cluster_emission_df[cluster_emission_df$result %in% c("Observed Negative True Positive"),]
  
  misclassification_df$invasive_status <- factor(misclassification_df$invasive_status,
                                       levels = c("Colonizing","Invasive"))
  
  # color schemes
  
  invasive_fill <- setNames(c("#990000",
                              "#011F5B"),
                            c("Invasive",
                              "Colonizing"))
  
  invasive_color <- setNames(c("#420D09",
                               "#29465B"),
                             c("Invasive",
                               "Colonizing"))
    
  # misclassification
  misclassification_p_plot <- ggplot(misclassification_df,
                              aes(x = invasive_status,
                                  y = probability,
                                  fill = invasive_status,
                                  color = invasive_status)) +
    
    geom_bar(stat = "identity",
             position = position_dodge(),
             alpha = 0.5) +
    geom_text(aes(label = round(probability,
                                digits = 5),
                  color = invasive_status),
              position = position_stack(vjust= 0.5),
              size = 7.5) + 
    {
      if (if_scientific){
        list(
          scale_y_continuous(labels = scales::scientific)
        )
      }
    } +
    scale_color_manual(values = invasive_color) +
    scale_fill_manual(values = invasive_fill) +
    geom_errorbar(aes(ymin = low, ymax = up),
                  width = 0.25,
                  position = position_dodge(0.9)) +
    theme_minimal() +
    labs(title = paste0(definition," HMM False Negative"),
         x = "Invasive Status", 
         y = "%Probability") +
    theme(plot.title = element_text(hjust = 0.5,
                                    face = "bold",
                                    size = 25),
          strip.text = element_blank(),
          axis.text.x = element_text(face = "bold",
                                     size = 20),
          axis.text.y = element_text(face = "bold",
                                     size = 20),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face = "bold",
                                      size = 25),
          axis.ticks.y = element_line(colour = "black",
                                      linewidth = 1),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  
  
  
  pdf(file = paste0(definition,"_Misclassification_Probability.pdf"),
      width=5,
      height=5)
  print(misclassification_p_plot)
  dev.off()
  
}

visualize_misclassification(cluster_msm_results = culture_msm_results,
                            definition = "Culture")

visualize_misclassification(cluster_msm_results = phylogenomics_msm_results,
                            definition = "Phylogenomics")
