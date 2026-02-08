req_cols <- function(df, cols) {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) {
    stop("Missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }
}

as_numeric_safe <- function(x) suppressWarnings(as.numeric(x))

clean_names <- function(nm) {
  nm <- gsub("^\ufeff", "", nm) # remove UTF-8 BOM if present
  trimws(nm)
}

read_table_smart <- function(path) {
  first <- readLines(path, n = 1, warn = FALSE, encoding = "UTF-8")
  
  if (grepl("\t", first, fixed = TRUE)) {
    df <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (grepl(";", first, fixed = TRUE)) {
    df <- utils::read.csv2(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
  
  names(df) <- clean_names(names(df))
  df
}

clean_day1 <- function(df) {
  names(df) <- clean_names(names(df))
  req_cols(df, c("female_id", "age_day", "eggs"))
  
  df$female_id <- trimws(as.character(df$female_id))
  df$female_id[df$female_id == ""] <- NA
  
  df$age_day <- as_numeric_safe(df$age_day)
  df$eggs    <- as_numeric_safe(df$eggs)
  
  df <- df[stats::complete.cases(df[, c("female_id", "age_day", "eggs")]), , drop = FALSE]
  df <- df[df$age_day >= 1, , drop = FALSE]
  df
}