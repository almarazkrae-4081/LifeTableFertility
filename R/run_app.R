#' Run the LifeTableFertility Shiny App
#'
#' Launches the bundled Shiny application included in this package.
#'
#' @return Runs a Shiny app (interactive).
#' @export
LifeTableFertility <- function() {
  app_dir <- system.file("app", package = "LifeTableFertility")
  if (app_dir == "") {
    stop("Could not find app directory. Try reinstalling the package.", call. = FALSE)
  }
  shiny::runApp(app_dir, display.mode = "normal")
}
