#use the output from checkM, mash, and mlst to do the QC of the NICU genome assemblies
library(openxlsx2)
library(dplyr)
library(Biostrings)
library(here)
# Paths to the output ----
checkm_output_path <- "./input/checkm_results.txt"
mlst_output_path <- "./input/nicu_mlst.csv"
mash_output_path <- "./input/Mash_contamination_report.txt"
assembly_scan_path <- "./input/assembly_scan/"
# CheckM ----
# Read and combine all input files
checkm_df <- read.delim(checkm_output_path,
                        header = TRUE) %>%
  select(Bin.Id,
         Marker.lineage,
         Completeness,
         Contamination,
         Strain.heterogeneity)

colnames(checkm_df) <- c("Genome",
                         "CheckM_Linage",
                         "CheckM_Completeness",
                         "CheckM_Contamination",
                         "CheckM_Heterogeneity")

checkm_df[c("CheckM_Completeness",
            "CheckM_Contamination",
            "CheckM_Heterogeneity")] <- lapply(checkm_df[c("CheckM_Completeness",
                                                           "CheckM_Contamination",
                                                           "CheckM_Heterogeneity")],
                                               as.numeric)


# MLST ----
mlst_df <- read.csv(mlst_output_path,
                    header = FALSE)
              
mlst_df <- mlst_df %>% 
  select(V1,V2,V3) %>%
  filter(!grepl("\\(|\\)",V2)) %>%
  filter(!grepl("\\(|\\)",V1))

colnames(mlst_df) <- c("Genome",
                       "MLST_Species",
                       "MLST_ST")

mlst_df$Genome <- sapply(X = 1:nrow(mlst_df),function(X){
  gsub(".fasta|.fa",
       "",
       strsplit(mlst_df$Genome[X],split = "/")[[1]][length(strsplit(mlst_df$Genome[X],split = "/")[[1]])])
})
# Mash ----
#use Ahmed's script: 
#Mash_contamination_checker.py to check contamination
mash_df <- read.delim(mash_output_path,
                      header = TRUE)

mash_df <- mash_df %>% select(File,
                              Conclusion,
                              Contamination)

colnames(mash_df) <- c("Genome",
                       "Mash_Conclusion",
                       "Mash_Contamination")

mash_df$Genome <- gsub("_sorted_mash.tab|_sorted_rebinned_mash.tab",
                       "",
                       mash_df$Genome)

#if there is Staphylococcus sp. Scaffold, this would not be considered the contamination
mash_df$Mash_Conclusion <- ifelse(
  mash_df$Mash_Conclusion == "Failed" &
    grepl("Staphylococcus sp.", mash_df$Mash_Contamination) &
    grepl("Scaffold", mash_df$Mash_Contamination),
  "Passed",
  mash_df$Mash_Conclusion
)

#combine all data frames into one
nicu_qc <- merge(checkm_df,
                 mlst_df,
                 by = "Genome",
                 all.x = TRUE,
                 all.y = FALSE)

nicu_qc <- merge(nicu_qc,
                 mash_df,
                 by = "Genome",
                 all.x = TRUE,
                 all.y = FALSE)

remove(checkm_df)
remove(mlst_df)
remove(mash_df)

# Assembly-scan ----
nicu_qc <- nicu_qc %>% mutate(genome_size = NA,
                              gc_content = NA,
                              average_coverage = NA,
                              num_contigs = NA,
                              n50 = NA,
                              average_contig_length = NA
                              )

assembly_scan_list <- as.data.frame(list.files(path = assembly_scan_path,
                                               pattern = "*assembly_scan.txt",
                                               all.files = TRUE,
                                               full.names = TRUE,
                                               recursive = FALSE))
assembly_scan_list <- assembly_scan_list %>% mutate(genome = NA)
colnames(assembly_scan_list) <- c("path",
                                  "Genome")

assembly_scan_list$Genome <- sapply(X = 1:nrow(assembly_scan_list),
                                    function(X){
                                      gsub("_assembly_scan.txt",
                                           "",
                                           strsplit(assembly_scan_list$path[X],
                                               split = "/")[[1]][length(strsplit(assembly_scan_list$path[X],
                                                                                 split = "/")[[1]])])
                                    })

for(i in 1:nrow(nicu_qc)){
  assembly_scan_df <- read.table(assembly_scan_list$path[assembly_scan_list$Genome == nicu_qc$Genome[i]],
                                 header = FALSE)
  nicu_qc$num_contigs[i] <- as.numeric(assembly_scan_df$V3[assembly_scan_df$V2 == "total_contig"])
  nicu_qc$average_contig_length[i] <- as.numeric(assembly_scan_df$V3[assembly_scan_df$V2 == "mean_contig_length"])
  nicu_qc$n50[i] <- as.numeric(assembly_scan_df$V3[assembly_scan_df$V2 == "n50_contig_length"])
  nicu_qc$genome_size[i] <- as.numeric(assembly_scan_df$V3[assembly_scan_df$V2 == "total_contig_length"])
  nicu_qc$gc_content[i] <- (as.numeric(assembly_scan_df$V3[assembly_scan_df$V2 == "contig_percent_c"]) + as.numeric(assembly_scan_df$V3[assembly_scan_df$V2 == "contig_percent_g"])) / (as.numeric(assembly_scan_df$V3[assembly_scan_df$V2 == "contig_percent_c"]) + as.numeric(assembly_scan_df$V3[assembly_scan_df$V2 == "contig_percent_g"]) + as.numeric(assembly_scan_df$V3[assembly_scan_df$V2 == "contig_percent_a"]) + as.numeric(assembly_scan_df$V3[assembly_scan_df$V2 == "contig_percent_t"]))
  remove(assembly_scan_df)
  
  if(grepl("marc.bacteremia",nicu_qc$Genome[i]) | 
     grepl("Path_Reg",nicu_qc$Genome[i]) | 
     grepl("Sta.aur.bead",nicu_qc$Genome[i]) | 
     grepl("Staph.bead",nicu_qc$Genome[i]) |
     grepl("NIICU.Sau.",nicu_qc$Genome[i])){
    contig_coverage <- names(Biostrings::readDNAStringSet(paste0("~/lab/cluster/NICU/genome_assemblies/",
                                                                 nicu_qc$Genome[i],
                                                                 ".fasta")))
    contig_coverage <- sapply(X = seq_len(length(contig_coverage)),
                              function(X){
                                as.numeric(strsplit(contig_coverage[X],
                                         split = "\\_")[[1]][length(strsplit(contig_coverage[X],
                                                                             split = "\\_")[[1]])])
                              })
    nicu_qc$average_coverage[i] <- mean(contig_coverage)
    remove(contig_coverage)
  }else{
    tmp <- read.delim("~/lab/cluster/NICU/microcart_reads/metadata/microcart_summary.tsv",
             header = TRUE)
    nicu_qc$average_coverage[i] <- tmp$Avg.Coverage[tmp$Sample == nicu_qc$Genome[i]]
    remove(tmp)
  }
}
remove(assembly_scan_list)

nicu_qc$genome_size <- as.numeric(nicu_qc$genome_size)
nicu_qc$gc_content <- as.numeric(nicu_qc$gc_content)
nicu_qc$CheckM_Completeness <- as.numeric(nicu_qc$CheckM_Completeness)
nicu_qc$CheckM_Contamination <- as.numeric(nicu_qc$CheckM_Contamination)
nicu_qc$CheckM_Heterogeneity <- as.numeric(nicu_qc$CheckM_Heterogeneity)
nicu_qc$average_coverage <- as.numeric(nicu_qc$average_coverage)
openxlsx2::write_xlsx(nicu_qc,
                      file = "all_nicu_genomes.xlsx")
saveRDS(nicu_qc,
        "all_nicu_genomes.RDS")
# Filtering ----

nicu_qc_passed <- nicu_qc %>% filter(CheckM_Linage == "g__Staphylococcus (UID301)" & 
                                       CheckM_Contamination <= 5 & 
                                       CheckM_Completeness >= 95 & 
                                       MLST_ST != "2250" &
                                       MLST_ST != "1223" &
                                       genome_size > 2550000 & 
                                       genome_size < 3150000 &
                                       Mash_Conclusion == "Passed")

openxlsx2::write_xlsx(nicu_qc_passed,
                      file = "passed_nicu_genomes.xlsx")

saveRDS(nicu_qc_passed,
        "passed_nicu_genomes.RDS")
