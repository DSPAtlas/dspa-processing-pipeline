# ---- CONFIG ----
yaml_dir <- "../yaml_files_sanitized"
script_path <- "../dspa-processing-pipeline/LiPTrPpipeline.R"
max_jobs <- 10

# ---- SETUP ----
files <- list.files(yaml_dir, pattern = "\\.yaml$", full.names = TRUE)

if (length(files) == 0) {
  stop("No YAML files found.")
}

library(callr)

cat("Found", length(files), "files\n")

# ---- PROCESS MANAGEMENT ----
procs <- list()

for (file in files) {

  # wait until a slot is free
  while (length(procs) >= max_jobs) {
    Sys.sleep(2)
    procs <- procs[sapply(procs, function(p) p$is_alive())]
  }

  # start new process
  p <- r_bg(function(script, input) {
    system2("Rscript", c(script, input))
  }, args = list(script_path, file))

  procs[[length(procs) + 1]] <- p

  cat("Started:", basename(file), "\n")
}

# ---- WAIT FOR COMPLETION ----
while (any(sapply(procs, function(p) p$is_alive()))) {
  Sys.sleep(2)
}

cat("All jobs finished\n")
