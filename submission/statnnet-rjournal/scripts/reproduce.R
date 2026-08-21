args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)

if (length(script_arg) == 1L) {
  submission_dir <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
} else {
  submission_dir <- normalizePath(".")
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
    "Install the required packages before reproducing the article: ",
    paste(required_packages[!available], collapse = ", ")
  )
}

source_file <- file.path(submission_dir, "mcinerney-burke.Rmd")
message("Rendering ", source_file)

rmarkdown::render(
  source_file,
  output_format = "all",
  clean = TRUE,
  envir = new.env(parent = globalenv())
)

expected_outputs <- file.path(
  submission_dir,
  c(
    "mcinerney-burke.html",
    "mcinerney-burke.pdf",
    "mcinerney-burke.tex",
    "mcinerney-burke.R",
    "RJwrapper.tex",
    "RJournal.sty"
  )
)

if (!all(file.exists(expected_outputs))) {
  stop(
    "Rendering did not create: ",
    paste(basename(expected_outputs[!file.exists(expected_outputs)]), collapse = ", ")
  )
}

message("Submission article reproduced successfully.")
