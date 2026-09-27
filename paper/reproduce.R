args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)

if (length(script_arg) == 1L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg))
  repository_dir <- dirname(dirname(script_path))
} else {
  repository_dir <- normalizePath(".")
}

required_packages <- c("rmarkdown", "rjtools", "nnet", "statnnet")
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

if (!rmarkdown::pandoc_available() && identical(.Platform$OS.type, "windows")) {
  pandoc_candidates <- c(
    "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools",
    "C:/Program Files/Positron/resources/app/quarto/bin/tools",
    "C:/Program Files/Quarto/bin/tools"
  )
  pandoc_candidates <- pandoc_candidates[
    file.exists(file.path(pandoc_candidates, "pandoc.exe"))
  ]
  if (length(pandoc_candidates) > 0L) {
    Sys.setenv(RSTUDIO_PANDOC = pandoc_candidates[[1L]])
  }
}

if (!rmarkdown::pandoc_available()) {
  stop("Pandoc is required to render the paper but was not found.")
}

source_file <- file.path(repository_dir, "paper", "new-draft.Rmd")
message("Rendering ", source_file)

rmarkdown::render(
  source_file,
  output_format = "rjtools::rjournal_pdf_article",
  clean = TRUE,
  envir = new.env(parent = globalenv())
)

intermediate_files <- c(
  file.path(
    repository_dir,
    "paper",
    c(
      "new-draft.R", "new-draft.html", "new-draft.tex",
      "new-draft.knit.md", "new-draft.utf8.md", "session-info.txt"
    )
  ),
  list.files(
    file.path(repository_dir, "paper"),
    pattern = "\\.(aux|bbl|blg|brf|log|out|toc)$",
    full.names = TRUE
  )
)
unlink(intermediate_files[file.exists(intermediate_files)])

intermediate_directories <- file.path(
  repository_dir,
  "paper",
  c("new-draft_files", "new_draft_files")
)
unlink(
  intermediate_directories[dir.exists(intermediate_directories)],
  recursive = TRUE,
  force = TRUE
)

expected_outputs <- file.path(repository_dir, "paper", "new-draft.pdf")

if (!all(file.exists(expected_outputs))) {
  stop(
    "Paper rendering did not create: ",
    paste(basename(expected_outputs[!file.exists(expected_outputs)]), collapse = ", ")
  )
}

message("Paper reproduction completed successfully.")
