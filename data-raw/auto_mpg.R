# Rebuild data/auto_mpg.rda from the UCI Machine Learning Repository.
# Source: https://doi.org/10.24432/C5859H (CC BY 4.0)

local({
  source_url <- "https://archive.ics.uci.edu/static/public/9/data.csv"
  source_file <- tempfile(fileext = ".csv")
  on.exit(unlink(source_file), add = TRUE)

  utils::download.file(source_url, source_file, mode = "wb")
  raw_auto_mpg <- utils::read.csv(
    source_file,
    na.strings = c("NaN", "?"),
    check.names = FALSE
  )

  expected_names <- c(
    "car_name", "cylinders", "displacement", "horsepower", "weight",
    "acceleration", "model_year", "origin", "mpg"
  )
  stopifnot(
    identical(names(raw_auto_mpg), expected_names),
    nrow(raw_auto_mpg) == 398L
  )

  auto_mpg <- raw_auto_mpg[
    stats::complete.cases(raw_auto_mpg$horsepower),
    c(
      "mpg", "cylinders", "displacement", "horsepower", "weight",
      "acceleration", "model_year", "origin"
    )
  ]
  auto_mpg$cylinders <- factor(auto_mpg$cylinders)
  auto_mpg$origin <- factor(
    auto_mpg$origin,
    levels = 1:3,
    labels = c("USA", "Europe", "Japan")
  )
  rownames(auto_mpg) <- NULL

  stopifnot(
    nrow(auto_mpg) == 392L,
    !anyNA(auto_mpg),
    identical(levels(auto_mpg$origin), c("USA", "Europe", "Japan"))
  )

  save(
    auto_mpg,
    file = file.path("data", "auto_mpg.rda"),
    compress = "xz",
    version = 2
  )
})
