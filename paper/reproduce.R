args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)

if (length(script_arg) == 1L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg))
  repository_dir <- dirname(dirname(script_path))
} else {
  repository_dir <- normalizePath(".")
}

required_packages <- c("rmarkdown", "rjtools", "nnet", "pdp", "statnnet")
available <- vapply(
  required_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)

if (!all(available)) {
  stop(
    "Install the required packages before reproducing the paper: ",
    paste(required_packages[!available], collapse = ", ")
  )
}

source_file <- file.path(repository_dir, "paper", "new-draft.Rmd")
message("Rendering ", source_file)

rmarkdown::render(
  source_file,
  output_format = "all",
  clean = TRUE,
  envir = new.env(parent = globalenv())
)

expected_outputs <- file.path(
  repository_dir,
  "paper",
  c("new-draft.html", "new-draft.pdf", "new-draft.tex", "new-draft.R")
)

if (!all(file.exists(expected_outputs))) {
  stop(
    "Paper rendering did not create: ",
    paste(basename(expected_outputs[!file.exists(expected_outputs)]), collapse = ", ")
  )
}

message("Paper reproduction completed successfully.")
