# load arguments
# use bash: Rscript pipeline.R params.yaml

# LiP MS pipeline
# with optional Tryptic control
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

# Access the parameters
group_id <- params$group_id
input_file <- params$input_file
input_file_tryptic_control <- params$input_file_tryptic_control
experiment_ids <- params$dpx_comparison
treatment <- params$treatment
ref_condition <- params$ref_condition
comparisons <- params$comparison
output_dir <- params$output_dir

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

# Log R session information, including loaded packages
cat("Logging R session info:\n")
sessionInfo() 

# Create a list to store the plots
plot_list <- list()
plot_list2 <- list()

# load file
read_file <-
  function(filename, ...) {
    data.table::fread(file = filename, ...) %>%
      janitor::clean_names() %>%
      tibble::as_tibble()
  }
df <- read_file(input_file)

# ------------------------------------------------------------------------------
# LiP 
# ------------------------------------------------------------------------------

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
  data =df,
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

## correlation_map
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
  distinct(r_file_name, fg_id, normalised_intensity_log2, eg_modified_peptide, pep_stripped_sequence, pg_protein_accessions,  r_condition, start, end, coverage, length) %>% 
  tidyr::complete(nesting(r_file_name, r_condition), nesting(pg_protein_accessions, fg_id, eg_modified_peptide, pep_stripped_sequence)) %>%
  dplyr::mutate(imputed = is.na(normalised_intensity_log2)
  )

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

# sum up precursors to peptide level and keep only one entry per pep_stripped_sequence
df %<>% protti::calculate_protein_abundance(
  sample = r_file_name,
  protein_id = eg_modified_peptide,
  precursor = fg_id,
  peptide = eg_modified_peptide,
  intensity_log2 = normalised_intensity_log2,
  min_n_peptides = 1,
  method = "sum",
  for_plot = FALSE,
  retain_columns = c("pg_protein_accessions", "r_condition", "start", "end", "coverage")
)

dia_clean_file <- file.path(group_folder_path, paste0("dia_clean_uniprot.tsv"))
write.table(df, dia_clean_file, sep = "\t", row.names= FALSE, quote = FALSE)



# ------------------------------------------------------------------------------
# TRYPTIC CONTROL
# ------------------------------------------------------------------------------

# Tryptic Control (if provided)
if (!is.null(input_file_tryptic_control)) {
  
  input_file_tryptic_control <- params$input_file_tryptic_control
  
  df_tryptic <- read_file(input_file_tryptic_control) %>%
    dplyr::mutate(
      pg_protein_accessions_split = ifelse(base::grepl(";", pg_protein_accessions, fixed = FALSE), 
      base::sort(base::strsplit(pg_protein_accessions, ";", fixed = TRUE)[[1]])[1], pg_protein_accessions)) %>%
      dplyr::mutate(intensity_log2 = log2(fg_ms2raw_quantity)) %>%
      dplyr::mutate(condrep = paste(r_condition, r_replicate, sep = "_"))
  
  ids <- df_tryptic %>% pull(pg_protein_accessions_split) %>% unique()
  uniprot <-
    protti::fetch_uniprot(
      ids,
      columns = c("length", "sequence")) 
  
  df_tryptic %<>%
    dplyr::left_join(uniprot, by = c("pg_protein_accessions_split" = "accession")) %>%
    protti::find_peptide(sequence, pep_stripped_sequence) %>%
    protti::assign_peptide_type(aa_before, last_aa, aa_after) %>%
    distinct() %>% 
    protti::calculate_sequence_coverage(
      protein_sequence = sequence, 
      peptides = pep_stripped_sequence
    ) %>% 
    dplyr::mutate(imputed = if_else(is.na(intensity_log2), TRUE, FALSE)) %>%
    dplyr::filter(intensity_log2 > 10) %>% 
    protti::normalise(
      sample = r_file_name,
      intensity_log2 = intensity_log2,
      method = "median"
    ) %>% 
    dplyr::mutate(normalised_intensity = 2^normalised_intensity_log2) %>%
    dplyr::mutate(fg_id = paste0(fg_labeled_sequence,fg_charge))

  # make the QC plots as for LiP
  plot_list2[[1]] <- protti::qc_ids(
    data = df_tryptic, 
    sample = condrep, 
    grouping = pep_grouping_key, 
    condition = r_condition, 
    intensity = fg_quantity
  ) + ggtitle('Tryptic control: Precursor ID count per sample')
  
  plot_list2[[2]] <- protti::qc_ids(
    data = df_tryptic, 
    sample = condrep, 
    grouping = pg_protein_accessions, 
    condition = r_condition, 
    intensity = pg_quantity
  ) + ggtitle('Tryptic control: Protein ID count per sample')
  
  plot_list2[[3]] <- protti::qc_intensity_distribution(
    data = df_tryptic,
    sample = condrep,
    grouping = fg_id,
    intensity_log2 = intensity_log2,
    plot_style = "histogram"
  ) + ggtitle('Tryptic control: Overall log2 Intensity distribution before normalisation')
  
  plot_list2[[4]] <- protti::qc_intensity_distribution(
    data = df_tryptic,
    sample = condrep,
    grouping = pep_grouping_key,
    intensity_log2 = normalised_intensity_log2,
    plot_style = "boxplot"
  ) + ggtitle('Tryptic control: Run intensities after normalisation')
  
  plot_list2[[5]] <- protti::qc_cvs(
    data = df_tryptic,
    grouping = fg_id,
    condition = r_condition,
    intensity = normalised_intensity,
    plot_style = "violin",
    plot = TRUE
  ) + ggtitle('Tryptic control: Violon plot')
  
  ## Missed cleavages
  plot_list2[[6]] <- protti::qc_missed_cleavages(
    data = df_tryptic,
    sample = condrep,
    grouping = fg_id,
    missed_cleavages = pep_nr_of_missed_cleavages,
    intensity = normalised_intensity,
    method = "count",
    plot = TRUE,
    interactive = FALSE
  ) + ggtitle('Tryptic control: number of missed cleavages')
  
  plot_list2[[7]] <- protti::qc_missed_cleavages(
    data = df_tryptic,
    sample = condrep,
    grouping = fg_id,
    missed_cleavages = pep_nr_of_missed_cleavages,
    intensity = normalised_intensity,
    method = "intensity",
    plot = TRUE,
    interactive = FALSE
  ) + ggtitle('Tryptic control: intensity of missed cleavages')
  
  plot_list2[[8]] <- protti::qc_peptide_type(
    df_tryptic,
    condrep,
    fg_id,
    pep_type,
    intensity = normalised_intensity,
    method = "count",
    plot = TRUE,
    interactive = FALSE
  ) + ggtitle('Tryptic control: number of half tryptics')
  
  plot_list2[[9]] <- protti::qc_peptide_type(
    df_tryptic,
    condrep,
    fg_id,
    pep_type,
    intensity = normalised_intensity,
    method = "intensity",
    plot = TRUE,
    interactive = FALSE
  ) + ggtitle('Tryptic control: intensity of half tryptics')
  
  df_tryptic %<>%
    distinct(r_file_name, fg_id, normalised_intensity_log2, eg_modified_peptide, pep_stripped_sequence, pg_protein_accessions, r_condition, condrep) %>% 
    tidyr::complete(
      nesting(r_file_name, r_condition, condrep), 
      nesting(pg_protein_accessions, fg_id, eg_modified_peptide, pep_stripped_sequence)
    )

  # calculate protein abundance
  df_tryptic %<>% calculate_protein_abundance(
    sample = r_file_name,
    protein_id = pg_protein_accessions,
    precursor = fg_id,
    peptide = eg_modified_peptide,
    intensity_log2 = normalised_intensity_log2,
    method = "iq",
    for_plot = FALSE,
    retain_columns = c("pg_protein_accessions", "r_condition", "condrep", "eg_modified_peptide")
  )
  
  write.table(df_tryptic, file.path(group_folder_path, "tryptic_control_clean.tsv"), sep = "\t", row.names= FALSE, quote = FALSE)
  
  ## Principal component analysis (PCA)
  plot_list2[[10]] <- df_tryptic %>%
    protti::qc_pca(
      sample = condrep, 
      grouping = pg_protein_accessions, 
      intensity = normalised_intensity_log2, 
      condition = r_condition
    ) + ggtitle('Tryptic control: PCA at protein level')
  
  ## correlation_map
  plot_list2[[11]] <- protti::qc_sample_correlation(
    data = df_tryptic,
    sample = condrep,
    grouping = pg_protein_accessions,
    intensity_log2 = normalised_intensity_log2,
    condition = r_condition
  )[[4]] 
  
  # Save QC plots
  output_tryptic_qc_pdf <- file.path(group_folder_path, "tryptic_controls_qc_plots.pdf")
  ggsave(
    filename = output_tryptic_qc_pdf, 
    plot = marrangeGrob(plot_list2, nrow=1, ncol=1), 
    width = 8, height = 6
  )
  rm(plot_list2, uniprot)
  gc()
 
}


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
      ref_condition = ref_condition,
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
  

  if (!is.null(input_file_tryptic_control)) {
    ref_condition_tryptic <- params$ref_condition_trp
    comparison_parts_tryptic <- strsplit(comparison_filter, "_vs_")[[1]] 
    comparison_parts_tryptic <- gsub("LiP", "TrP", comparison_parts_tryptic)
    
    df_trp_filtered <- df_tryptic %>%
      dplyr::mutate(
        pg_protein_accessions_split = ifelse(base::grepl(";", pg_protein_accessions, fixed = FALSE), 
         base::sort(base::strsplit(pg_protein_accessions, ";", fixed = TRUE)[[1]])[1], pg_protein_accessions))  %>% 
      dplyr::filter(r_condition %in% comparison_parts_tryptic) %>%
      dplyr::distinct(pg_protein_accessions_split, r_file_name, r_condition, .keep_all = TRUE)
    
    df_trp_filtered_diff <- df_trp_filtered  %>%
      assign_missingness(
        sample = r_file_name,
        condition = r_condition,
        intensity = normalised_intensity_log2,
        grouping =  pg_protein_accessions,
        ref_condition = ref_condition_tryptic,
        retain_columns = all_of(c("pg_protein_accessions","r_file_name", "r_condition", 
                                  "normalised_intensity_log2", "pg_protein_accessions_split",
                                   "eg_modified_peptide"))
      ) %>%
      calculate_diff_abundance(
        sample = r_file_name,
        condition = r_condition,
        grouping =  pg_protein_accessions,
        intensity_log2 = normalised_intensity_log2,
        comparison = comparison,
        method = "t-test",
        retain_columns = all_of(c("pg_protein_accessions","r_file_name", "r_condition", 
                                  "normalised_intensity_log2", 
                                  "pg_protein_accessions_split", "eg_modified_peptide"))
      )
    
    # Save Differential Abundance results
    diff_trp_file_path <- file.path(
      group_folder_path, 
      paste0("trp_differential_abundance_", experiment_ids, "_", comparisons, ".tsv")
    )
    write.table( df_trp_filtered_diff, diff_trp_file_path, sep = "\t", row.names= FALSE, quote = FALSE)
    
    # perform TrP protein correction on LiP:
    df_trp_filtered_diff %<>%
      dplyr::mutate(
        pg_protein_accessions_split = ifelse(base::grepl(";", pg_protein_accessions, fixed = FALSE), 
        base::sort(base::strsplit(pg_protein_accessions, ";", fixed = TRUE)[[1]])[1], pg_protein_accessions)) %>%
      dplyr::mutate(Trp_candidates = adj_pval < 0.05 & abs(diff) > 1) %>%
      dplyr::mutate(comparison = gsub("_TrP", "_LiP", x = comparison, fixed = TRUE)) 
  
    df_diff <- protti::correct_lip_for_abundance(
      lip_data = df_diff,
      trp_data =  df_trp_filtered_diff,
      protein_id = pg_protein_accessions,
      grouping = eg_modified_peptide,
      comparison = comparison, 
      diff = diff,
      n_obs = n_obs,
      std_error = std_error,
      p_adj_method = "BH",
      retain_columns = all_of(c("missingness")),
      method = "satterthwaite"
    )
  }
  df_diff <- df_diff %>%
    left_join(df %>%
                select(eg_modified_peptide, start, end) %>%
                distinct(eg_modified_peptide, .keep_all = TRUE),
              by = "eg_modified_peptide")
  
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
    df_go_term <- map_dfr(go_terms, function(go_col) {
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


# copy yaml file into the output as well
yaml_file_path <- file.path(group_folder_path, "params.yaml")
file.copy(yaml_file, yaml_file_path)

if (sink.number() > 0) sink(NULL)
closeAllConnections()

#joined <-df %>%
#  dplyr::select(coverage, end, start, eg_modified_peptide) %>%
#  dplyr::left_join(
#  df_diff, 
#  by = "eg_modified_peptide"
#)
