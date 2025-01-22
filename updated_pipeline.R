# From Ludo: 
# this is a test modification note to see whether my fork is working correctly

# load arguments
# use bash: Rscript pipeline.R params.yaml

library(yaml)
library(protti)
library(dplyr)
library(missForest)
library(doParallel)

registerDoParallel(cores = 10) 

# Read the YAML file
args <- commandArgs(trailingOnly = TRUE)
# Expect the first argument to be the path to the YAML file
yaml_file <- args[1]
params <- yaml::read_yaml(yaml_file)

# Access the parameters
group_id <- params$group_id
input_file <- params$input_file
input_file_tryptic_control <- params$input_file_tryptic_control
experiment_ids <- params$experiment_id
treatment <- params$treatment
ref_condition <- params$ref_condition
comparisons <- params$comparison
output_dir <- params$output_dir

if (!is.null(params$input_file_tryptic_control)) {
  input_file_tryptic_control <- params$input_file_tryptic_control
  tryptic_control_data <- protti::read_protti(input_file_tryptic_control)
} else {
  input_file_tryptic_control <- NULL
  tryptic_control_data <- NULL
}

group_folder_path <- file.path(output_dir, group_id)

if (!dir.exists("./preprocessed")) {
  dir.create("./preprocessed", recursive = TRUE)
}

if (!dir.exists(group_folder_path)) {
  dir.create(group_folder_path, recursive = TRUE)
  cat("Created group folder:", group_folder_path, "\n")
}

# Log outputs to check for errors 
# Redirect output to a log file
logfile_dir <- file.path(output_dir, group_id, "processing_log.txt")
sink(logfile_dir, append = TRUE)

# Print or use the parameters
print(paste("Input File:", input_file))
print(paste("Treatment:", treatment))
print(paste("Ref Condtion:", ref_condition))
print(paste("Output Directory:", output_dir))

# Log R session information, including loaded packages
cat("Logging R session info:\n")
sessionInfo() 

# Create a list to store the plots
plot_list <- list()

# load file
DIA_raw <- read_protti(input_file)

# ------------------------------------------------------------------------------
# LIP QUANT
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Preprocessing 
# ------------------------------------------------------------------------------

DIA_raw$intensity_log2 <- log2(DIA_raw$fg_ms2raw_quantity)
DIA_raw$condrep <- paste(DIA_raw$r_condition, DIA_raw$r_replicate, sep = "_")

DIA_raw_norm <- protti::normalise(
  DIA_raw,
  sample = r_file_name,
  intensity_log2 = intensity_log2,
  method = "median"
)

DIA_clean <- DIA_raw_norm %>%
  dplyr::filter(intensity_log2 > 10) 

DIA_clean$fg_id <- paste0(DIA_clean$fg_labeled_sequence, DIA_clean$fg_charge)
unis <- unique(DIA_clean$pg_protein_accessions) # make vector for fetch_uniprot

## Load data from uniprot and join with DIA dataframe
uniprot <-
  protti::fetch_uniprot(
    unis,
    columns =  c(
      "protein_name",
      "gene_names",
      "length",
      "sequence",
      "xref_pdb",
      "go_f",
      "go_p",
      "go_c"
    )
  ) 

DIA_clean_uniprot <- DIA_clean %>%
  left_join(uniprot, by = c("pg_protein_accessions" = "accession")) %>% 
  find_peptide(sequence, pep_stripped_sequence) %>%
  assign_peptide_type(aa_before, last_aa, aa_after) %>%
  distinct()

DIA_clean_uniprot$condrep <- paste(DIA_clean_uniprot$r_condition, DIA_clean_uniprot$r_replicate, sep = "_")
DIA_raw$condrep <- paste(DIA_raw$r_condition, DIA_raw$r_replicate, sep = "_")

proteins_identified <- uniprot %>%
  distinct(accession)

# Imputation on precursor lvl 
imputed <- impute_randomforest(
  DIA_clean_uniprot,
  sample = r_file_name,
  grouping = fg_id,
  intensity_log2 = normalised_intensity_log2,
  retain_columns = c("pep_stripped_sequence", "pg_protein_accessions", 
                     "gene_names", "go_f", "r_condition", "start", "end"),
  parallelize = "variables"
)

# sum up precursors to peptide lvl and keep only one entry per pep_stripped_sequence
DIA_clean_uniprot_summed_protti <- protti::calculate_protein_abundance(
  imputed,
  sample = r_file_name,
  protein_id = pep_stripped_sequence,
  precursor = fg_id,
  peptide = pep_stripped_sequence,
  intensity_log2 = normalised_intensity_log2,
  min_n_peptides = 1,
  method = "sum",
  for_plot = FALSE,
  retain_columns = c("pep_stripped_sequence", "pg_protein_accessions", 
                     "gene_names", "go_f", "r_condition", "start", "end")
)

dia_clean_file <- file.path(group_folder_path, paste0("dia_clean_uniprot.csv"))
write.csv(DIA_clean_uniprot_summed_protti, dia_clean_file)

# ------------------------------------------------------------------------------
# QC
# ------------------------------------------------------------------------------

plot_list[[1]] <- protti::qc_median_intensities(
  data = DIA_raw_norm,
  sample = r_file_name,
  grouping = pep_grouping_key,
  intensity = normalised_intensity_log2,
  plot = TRUE,
  interactive = FALSE
)

plot_list[[2]] <- protti::qc_cvs(
  data = DIA_clean_uniprot,
  grouping = fg_id,
  condition = r_condition,
  intensity = fg_quantity,
  plot = FALSE
)

plot_list[[3]] <- protti::qc_cvs(
  data = DIA_clean_uniprot,
  grouping = fg_id,
  condition = r_condition,
  intensity = fg_quantity,
  plot_style = "density",
  plot = TRUE
)

plot_list[[4]] <- protti::qc_cvs(
  data = DIA_clean_uniprot,
  grouping = fg_id,
  condition = r_condition,
  intensity = fg_quantity,
  plot_style = "violin",
  plot = TRUE
)

plot_list[[5]] <- protti::qc_intensity_distribution(
  data = DIA_clean_uniprot,
  sample = condrep,
  grouping = pep_grouping_key,
  intensity_log2 = intensity_log2,
  plot_style = "histogram"
)

## Missed cleavages
plot_list[[6]] <- protti::qc_missed_cleavages(
  data = DIA_clean_uniprot,
  sample = condrep,
  grouping = fg_id,
  missed_cleavages = pep_nr_of_missed_cleavages,
  intensity = fg_quantity,
  method = "intensity",
  plot = TRUE,
  interactive = FALSE
)

plot_list[[7]] <- protti::qc_missed_cleavages(
  data = DIA_clean_uniprot,
  sample = condrep,
  grouping = fg_id,
  missed_cleavages = pep_nr_of_missed_cleavages,
  intensity = fg_quantity,
  method = "count",
  plot = TRUE,
  interactive = FALSE
)

plot_list[[8]] <- protti::qc_peptide_type(
  DIA_clean_uniprot,
  condrep,
  fg_id,
  pep_type,
  intensity = fg_quantity,
  method = "count",
  plot = TRUE,
  interactive = FALSE
)

plot_list[[9]] <- protti::qc_peptide_type(
  DIA_clean_uniprot,
  condrep,
  fg_id,
  pep_type,
  intensity = fg_quantity,
  method = "intensity",
  plot = TRUE,
  interactive = FALSE
)

plot_list[[10]] <- protti::qc_ids(
  data = DIA_raw, 
  sample = condrep, 
  grouping = pep_grouping_key, 
  condition = r_condition, 
  intensity = fg_quantity
)

## Principal component analysis (PCA)
plot_list[[11]] <-DIA_clean_uniprot %>%
  protti::qc_pca(
    sample = condrep, 
    grouping = pep_grouping_key, 
    intensity = intensity_log2, 
    condition = r_condition
  )

## corelation_map
plot_list[[12]] <- protti::qc_sample_correlation(
  data = DIA_clean_uniprot,
  sample = r_file_name,
  grouping = fg_id,
  intensity_log2 = intensity_log2,
  condition = r_condition
)

# Save QC plots
output_qc_pdf <- file.path(group_folder_path, "qc_plots.pdf")
pdf(output_qc_pdf, width = 7, height = 5)  # Open the PDF device
# Loop through the list and print each ggplot object
lapply(plot_list, print)
dev.off()  

# ------------------------------------------------------------------------------
# Analysis ...
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# differential abundance and go term
# ------------------------------------------------------------------------------

# Data analysis
## Volcano plots peptide lvl
# Iterate over each comparison and associated experiment ID
for (i in seq_along(comparisons)) {
  comparison_filter <- comparisons[[i]]
  experiment_id <- experiment_ids[[i]]
  
  comparison_parts <- strsplit(comparison_filter, "_vs_")[[1]] 
  print(comparison_parts)
  
  filtered_data <- DIA_clean_uniprot_summed_protti %>%
    dplyr::filter(r_condition %in% comparison_parts)
  
  Volcano_input <- filtered_data %>%
    unique() %>%
    protti::assign_missingness(
      sample = r_file_name,
      condition = r_condition,
      grouping = pep_stripped_sequence,
      intensity = normalised_intensity_log2,
      ref_condition = ref_condition,
      retain_columns = all_of(c("pg_protein_accessions", "r_file_name", "r_condition", 
                                "normalised_intensity_log2", 
                                "pep_stripped_sequence", "go_f", "start", "end"))
    )
  
  df_diff_abundance <- protti::calculate_diff_abundance(
    data = Volcano_input,
    sample = r_file_name,
    condition = r_condition,
    grouping = pep_stripped_sequence,
    intensity_log2 = normalised_intensity_log2,
    missingness = missingness,
    comparison = comparison,
    method = "moderated_t-test",
    retain_columns = all_of(c("pg_protein_accessions", "pep_stripped_sequence", 
                              "comparison", "go_f", "start", "end"))
  )
  
  # Save Differential Abundance results
  diff_abundance_file <- file.path(
    group_folder_path, 
    paste0("differential_abundance_", experiment_id, "_", comparison_filter, ".csv")
  )
  write.csv(df_diff_abundance, diff_abundance_file)
  
  # Calculate GO term enrichment with error handling
  tryCatch({
    df_diff_abundance_significant <- df_diff_abundance %>%
      dplyr::mutate(significant = ifelse(!is.na(adj_pval) & adj_pval < 0.05, TRUE, FALSE)) 
    
    df_go_term <- protti::calculate_go_enrichment(
      data = df_diff_abundance_significant,
      protein_id = pg_protein_accessions,
      go_annotations_uniprot = go_f,
      is_significant = significant,
      min_n_detected_proteins_in_process = 3,
      plot = FALSE
    )
    
    # Save GO Term enrichment results
    go_term_file <- file.path(group_folder_path, paste0("go_term_", experiment_id, "_", comparison_filter, ".csv"))
    write.csv(df_go_term, go_term_file)
  }, error = function(e) {
    message(paste("Error in GO term enrichment for comparison", comparison_filter, ":", e))
  })
}

# ------------------------------------------------------------------------------
# LiP-Quant Analysis
# ------------------------------------------------------------------------------
# Adding dose-response analysis using `parallel_fit_drc_4p`

# Ensure a numeric concentration column exists for dose-response fitting
DIA_clean_uniprot_summed_protti$conc_frag <- as.numeric(gsub("_mM", "", DIA_clean_uniprot_summed_protti$r_condition))

# Perform parallel dose-response fitting
lipquant_results <- protti::parallel_fit_drc_4p(
  data = DIA_clean_uniprot_summed_protti,
  sample = r_file_name,
  grouping = pep_stripped_sequence,
  intensity_log2 = normalised_intensity_log2,
  conc_frag = conc_frag,
  filter = "post",
  replicate_completeness = 0.7,
  condition_completeness = 0.5,
  correlation_cutoff = 0.8,
  log_logarithmic = TRUE,
  retain_columns = c("pg_protein_accessions", "start", "end", "pep_type", "gene_names", "go_f")
)

# Save LiP-Quant results
lipquant_output_file <- file.path(group_folder_path, "lipquant_results.csv")
write.csv(lipquant_results, lipquant_output_file)


lipquant_plots <- protti::plot_drc_results(
  lipquant_results,
  grouping = pep_stripped_sequence,
  conc_frag = conc_frag,
  intensity_log2 = normalised_intensity_log2
)

# Save diagnostic plots for dose-response fitting
lipquant_plot_pdf <- file.path(group_folder_path, "lipquant_dose_response_plots.pdf")
pdf(lipquant_plot_pdf, width = 7, height = 5)
lapply(lipquant_plots, print)
dev.off()


# ------------------------------------------------------------------------------
# TRYPTIC CONTROL
# ------------------------------------------------------------------------------

# if tryptic control present
# Tryptic Control (if provided)
if (!is.null(input_file_tryptic_control)) {
  tryptic_control_data <- protti::read_protti(input_file_tryptic_control)
  tryptic_control_data$intensity_log2 <- log2(tryptic_control_data$pg_quantity)
  
  # Preprocess tryptic control data
  tryptic_clean <- tryptic_control_data %>%
   dplyr::filter(intensity_log2 > 10) %>%
    protti::normalise(
      sample = r_file_name,
      intensity_log2 = intensity_log2,
      method = "median"
    )
  
  # Fetch UniProt annotations for tryptic control
  unis_tryptic <- unique(tryptic_clean$pg_protein_accessions)
  uniprot_tryptic <- protti::fetch_uniprot(unis_tryptic, columns = c("protein_name", "gene_names", "length"))
  
  tryptic_clean_uniprot <- tryptic_clean %>%
    left_join(uniprot_tryptic, by = c("pg_protein_accessions" = "accession"))
  
  # Save preprocessed tryptic control data
  write.csv(tryptic_clean_uniprot, file.path(group_folder_path, "tryptic_control_clean.csv"))
}


# copy yaml file into the output as well
yaml_file_path <- file.path(group_folder_path, "params.yaml")
file.copy(yaml_file, yaml_file_path)

# Stop redirecting output to the log file
sink()
