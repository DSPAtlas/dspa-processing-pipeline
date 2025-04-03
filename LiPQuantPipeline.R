
# LiPQuant pipeline
closeAllConnections()
library(yaml)
library(protti)
library(tidyr)
library(dplyr)
library(missForest)
library(doParallel)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(magrittr)
library(purrr)

registerDoParallel(cores = 10) 

# Read the YAML file
args <- commandArgs(trailingOnly = TRUE)
# Expect the first argument to be the path to the YAML file
yaml_file <- args[1]
params <- yaml::read_yaml(file=yaml_file)

group_id <- params$group_id
input_file <- params$input_file
input_file_tryptic_control <- params$input_file_tryptic_control
experiment_ids <- params$dpx_comparison
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

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
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
plot_list2 <- list()

# load file2
read_protti2 <-
  function(filename, ...) {
    data.table::fread(file=filename, ...) %>%
      janitor::clean_names() %>%
      tibble::as_tibble()
  }
df <- read_protti2(filename=input_file)

# ------------------------------------------------------------------------------
# Preprocessing 
# ------------------------------------------------------------------------------

df %<>%
  dplyr::mutate(pg_protein_accessions_split = ifelse(base::grepl(";", pg_protein_accessions, fixed = FALSE), 
                                                     base::sort(base::strsplit(pg_protein_accessions, ";", fixed = TRUE)[[1]])[1], pg_protein_accessions)) %>%
  dplyr::mutate(intensity_log2 = log2(fg_ms2raw_quantity)) %>%
  dplyr::mutate(condrep = paste(r_condition, r_replicate, sep = "_"))

ids <- df %>% dplyr::pull(pg_protein_accessions_split) %>% base::unique()
uniprot <-
  protti::fetch_uniprot(
    ids,
    columns =  c("length", "sequence") ) 

# ------------------------------------------------------------------------------
# Normalise
# ------------------------------------------------------------------------------

df %<>%
  protti::normalise(
    sample = r_file_name,
    intensity_log2 = intensity_log2,
    method = "median"
  ) %>%
  dplyr::left_join(
    uniprot, 
    by = c("pg_protein_accessions_split" = "accession")
  ) %>%
  protti::find_peptide(sequence, pep_stripped_sequence) %>%
  protti::assign_peptide_type(aa_before, last_aa, aa_after) %>%
  dplyr::distinct() %>% 
  protti::calculate_sequence_coverage(
    protein_sequence = sequence, 
    peptides = pep_stripped_sequence
  ) %>% 
  dplyr::mutate(normalised_intensity = 2^normalised_intensity_log2) %>%
  dplyr::mutate(imputed = if_else(is.na(intensity_log2), TRUE, FALSE)) %>%
  dplyr::filter(intensity_log2 > 10) %>%
  dplyr::mutate(fg_id =paste0(fg_labeled_sequence, fg_charge)) %>%
  dplyr::mutate(normalised_intensity = 2^normalised_intensity_log2)

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------

plot_list[[1]] <- protti::qc_ids(
  data = df, 
  sample = condrep, 
  grouping = pep_grouping_key, 
  condition = r_condition, 
  intensity = fg_quantity
)+ ggtitle('Precursor ID count per sample')

plot_list[[2]] <- protti::qc_ids(
  data = df, 
  sample = condrep, 
  grouping = pg_protein_accessions, 
  condition = r_condition, 
  intensity = pg_quantity
)+ ggtitle('Protein ID count per sample')

plot_list[[3]] <- protti::qc_intensity_distribution(
  data = df,
  sample = condrep,
  grouping = pep_grouping_key,
  intensity_log2 = intensity_log2,
  plot_style = "histogram"
) + ggtitle('Overall log2 Intensity distribution before normalisation')

plot_list[[4]] <- protti::qc_intensity_distribution(
  data = df,
  sample = condrep,
  grouping = pep_grouping_key,
  intensity_log2 = normalised_intensity_log2,
  plot_style = "boxplot"
) + ggtitle('Run intensities after normalisation')

plot_list[[5]] <- protti::qc_cvs(
  data = df,
  grouping = fg_id,
  condition = r_condition,
  intensity = normalised_intensity,
  plot_style = "violin",
  plot = TRUE
)

## Missed cleavages
plot_list[[6]] <- protti::qc_missed_cleavages(
  data = df,
  sample = condrep,
  grouping = fg_id,
  missed_cleavages = pep_nr_of_missed_cleavages,
  intensity = normalised_intensity,
  method = "count",
  plot = TRUE,
  interactive = FALSE
)

plot_list[[7]] <- protti::qc_missed_cleavages(
  data = df,
  sample = condrep,
  grouping = fg_id,
  missed_cleavages = pep_nr_of_missed_cleavages,
  intensity = normalised_intensity,
  method = "intensity",
  plot = TRUE,
  interactive = FALSE
)

plot_list[[8]] <- protti::qc_peptide_type(
  df,
  condrep,
  fg_id,
  pep_type,
  intensity = normalised_intensity,
  method = "count",
  plot = TRUE,
  interactive = FALSE
)

plot_list[[9]] <- protti::qc_peptide_type(
  df,
  condrep,
  fg_id,
  pep_type,
  intensity = normalised_intensity,
  method = "intensity",
  plot = TRUE,
  interactive = FALSE
)

## Principal component analysis (PCA)
plot_list[[10]] <- df %>%
  protti::qc_pca(
    sample = condrep, 
    grouping = pep_grouping_key, 
    intensity = intensity_log2, 
    condition = r_condition
  )

## corelation_map
plot_list[[11]] <- protti::qc_sample_correlation(
  data = df,
  sample = condrep,
  grouping = fg_id,
  intensity_log2 = intensity_log2,
  condition = r_condition
)[[4]]


# ------------------------------------------------------------------------------
# Impute
# ------------------------------------------------------------------------------

df %<>% 
  distinct(r_file_name, fg_id, normalised_intensity_log2, eg_modified_peptide, pep_stripped_sequence, pg_protein_accessions,  r_condition) %>% 
  tidyr::complete(nesting(r_file_name, r_condition), nesting(pg_protein_accessions, fg_id, eg_modified_peptide, pep_stripped_sequence))

df %<>% impute_randomforest(
  sample = r_file_name,
  grouping = fg_id,
  intensity_log2 = normalised_intensity_log2,
  retain_columns = c("eg_modified_peptide", "pep_stripped_sequence", "pg_protein_accessions",
                     "r_condition", "start", "end", "imputed",
                     "coverage", "length"),
  parallelize = "variables"
)

imputed_file <- file.path(group_folder_path, paste0("imputed.tsv"))
write.table(df, imputed_file, sep = "\t", row.names= FALSE, quote = FALSE)

df %<>% protti::calculate_protein_abundance(
  sample = r_file_name,
  protein_id = eg_modified_peptide,
  precursor = fg_id,
  peptide = eg_modified_peptide,
  intensity_log2 = normalised_intensity_log2,
  min_n_peptides = 1,
  method = "sum",
  for_plot = FALSE,
  retain_columns = c("pg_protein_accessions", "r_condition")#  "start", "end", "coverage"
)


plot_list[[12]] <- df %>%
  dplyr::mutate(imputed = factor(imputed, levels = c(TRUE, FALSE), labels = c("Imputed", "Observed"))) %>%
  ggplot(aes(x = normalised_intensity_log2, fill = imputed)) +
  labs(title = "Histogram of Intensities Before and After Imputation (Log2)",
       x = "Log2 Intensity",
       y = "Frequency",
       fill = "Type") +
  geom_histogram(
    binwidth = 0.5,
    color = "black",
    position = "identity"
  ) +
  scale_fill_manual(values = protti_colours[c(2, 1)]) +
  theme_bw() + 
  coord_cartesian(xlim = c(5, 30))

# ------------------------------------------------------------------------------
# Save QC plots
# ------------------------------------------------------------------------------

output_qc_pdf <- file.path(group_folder_path, "qc_plots.pdf")

ggsave(
  filename = output_qc_pdf, 
  plot = marrangeGrob(plot_list, nrow=1, ncol=1), 
  width = 8, height = 6
)
rm(plot_list, uniprot)
gc()

# ------------------------------------------------------------------------------
# Differential analysis
# ------------------------------------------------------------------------------

for (i in seq_along(comparisons)) {
  comparison_filter <- comparisons[[i]]
  experiment_id <- experiment_ids[[i]]
  
  comparison_parts <- strsplit(comparison_filter, "_vs_")[[1]] 
  
  df_filtered <- df %>%
    dplyr::filter(r_condition %in% comparison_parts)
  
  df_diff <- df_filtered  %>%
    unique() %>%
    protti::assign_missingness(
      sample = r_file_name,
      condition = r_condition,
      grouping = eg_modified_peptide,
      intensity = normalised_intensity_log2,
      ref_condition = "CTR_LiP",
      retain_columns = all_of(c("pg_protein_accessions","r_file_name", "r_condition", 
                                "normalised_intensity_log2")))%>%
    protti::calculate_diff_abundance(
      sample = r_file_name,
      condition = r_condition,
      grouping = eg_modified_peptide,
      intensity_log2 = normalised_intensity_log2,
      missingness = missingness,
      comparison = comparison,
      method = "t-test",
      retain_columns = all_of(c("pg_protein_accessions","eg_modified_peptide", 
                                "comparison"))
    )
  
  diff_abundance_file <- file.path(
    group_folder_path, 
    paste0("differential_abundance_", experiment_id, "_", comparison_filter, ".tsv")
  )
  write.table(df_diff, diff_abundance_file, sep = "\t", row.names= FALSE, quote = FALSE)
  
  unis <- df_diff %>%
    dplyr::mutate(pg_protein_accessions_split = ifelse(base::grepl(";", pg_protein_accessions, fixed = FALSE), 
                                                       base::sort(base::strsplit(pg_protein_accessions, ";", fixed = TRUE)[[1]])[1], pg_protein_accessions)) %>% 
    pull(pg_protein_accessions_split)  %>%# make vector for fetch_uniprot
    unique()
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
  
  joined_df <- df_diff %>%
    dplyr::mutate(pg_protein_accessions_split = ifelse(
      base::grepl(";", pg_protein_accessions, fixed = FALSE), 
      base::sort(base::strsplit(pg_protein_accessions, ";", fixed = TRUE)[[1]])[1], 
      pg_protein_accessions
    )) %>%
    dplyr::left_join(uniprot, by = c("pg_protein_accessions_split" = "accession"))
  
  
  tryCatch({
    go_terms <- c("go_f", "go_p", "go_c")
    
    # Loop over the GO terms and combine results into a single data frame
    df_go_term <- purr::map_dfr(go_terms, function(go_col) {
      joined_df %>%
        dplyr::mutate(significant = ifelse(!is.na(adj_pval) & adj_pval < 0.05, TRUE, FALSE)) %>%
        protti::calculate_go_enrichment(
          protein_id = pg_protein_accessions,   
          go_annotations_uniprot = !!sym(go_col),  
          is_significant = significant,          
          min_n_detected_proteins_in_process = 3,
          plot=FALSE
        ) %>%
        dplyr::mutate(go_type = go_col)  # Add a column to indicate the GO type
    })
    
    
    # Save GO Term enrichment results
    go_term_file <- file.path(group_folder_path, paste0("go_term_", experiment_id, "_", comparison_filter, ".tsv"))
    write.table(df_go_term, go_term_file, sep = "\t", row.names= FALSE, quote = FALSE)
  }, error = function(e) {
    message(paste("Error in GO term enrichment for comparison", comparison_filter, ":", e))
  })
}


# ------------------------------------------------------------------------------
# LiP-Quant Analysis
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Dose Response
# ------------------------------------------------------------------------------
# Adding dose-response analysis using `parallel_fit_drc_4p`

# Ensure a numeric concentration column exists for dose-response fitting
df %<>%
  drop_na(r_condition) %>%
  dplyr::mutate(
    dose = case_when(
      r_condition == "CTR_LiP" ~ 0,  # CTR is assigned a dose of 0
      grepl(".*_\\d+\\.?\\d*uM", r_condition) ~ as.numeric(sub(".*_(\\d+\\.?\\d*)uM", "\\1", r_condition)),  # Extract number after any prefix
      TRUE ~ NA_real_  # For other cases, set as NA
    )
  ) 

# Perform parallel dose-response fitting
lipquant_results <- protti::parallel_fit_drc_4p(
  data = df,
  sample = r_file_name,
  grouping = eg_modified_peptide,
  response = normalised_intensity_log2,
  dose = dose,
  filter = "post",
  replicate_completeness = 0.7,
  condition_completeness = 0.5,
  correlation_cutoff = 0.8,
  log_logarithmic = TRUE,
  retain_columns = c("pg_protein_accessions") # "start", "end", "pep_type"
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


# copy yaml file into the output as well
yaml_file_path <- file.path(group_folder_path, "params.yaml")
file.copy(yaml_file, yaml_file_path)

# Stop redirecting output to the log file
if (sink.number() > 0) sink(NULL)
closeAllConnections()