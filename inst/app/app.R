# LifeTableFertility - Shiny app
# Location in package: inst/app/app.R
# Runs via: LifeTableFertility()

# ---------- Helpers ----------
req_cols <- function(df, cols) {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) {
    stop("Missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }
}

as_numeric_safe <- function(x) suppressWarnings(as.numeric(x))

clean_names <- function(nm) {
  nm <- gsub("^\ufeff", "", nm) # remove BOM if present
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

# ---------- Core ----------
build_life <- function(df_clean, development_days = 0, female_ratio = 0.5) {
  if (nrow(df_clean) < 5) stop("Too few valid rows after cleaning.", call. = FALSE)
  if (any(df_clean$age_day < 1, na.rm = TRUE)) stop("age_day must start at 1.", call. = FALSE)
  if (any(df_clean$eggs < 0, na.rm = TRUE)) stop("eggs cannot be negative.", call. = FALSE)
  if (!is.finite(development_days) || development_days < 0) stop("Immature development must be >= 0.", call. = FALSE)
  if (!is.finite(female_ratio) || female_ratio <= 0 || female_ratio > 1) stop("Female proportion must be in (0, 1].", call. = FALSE)
  
  ids <- unique(df_clean$female_id)
  n0 <- length(ids)
  if (n0 < 2) stop(paste0("At least 2 females are required. After cleaning, n_females = ", n0, "."), call. = FALSE)
  
  df_clean$age <- df_clean$age_day + development_days
  
  alive <- stats::aggregate(female_id ~ age, data = df_clean, FUN = function(x) length(unique(x)))
  names(alive)[2] <- "n_alive"
  alive$lx <- alive$n_alive / n0
  
  mx <- stats::aggregate(eggs ~ age, data = df_clean, FUN = function(x) mean(x, na.rm = TRUE))
  names(mx)[2] <- "mx_raw"
  mx$mx <- mx$mx_raw * female_ratio
  
  tab <- merge(alive[, c("age", "n_alive", "lx")], mx[, c("age", "mx_raw", "mx")], by = "age", all = TRUE)
  tab <- tab[order(tab$age), ]
  tab$lxmx  <- tab$lx * tab$mx
  tab$xlxmx <- tab$age * tab$lxmx
  
  R0 <- sum(tab$lxmx, na.rm = TRUE)
  if (!is.finite(R0) || R0 <= 0) {
    T <- NA_real_; r <- NA_real_; lambda <- NA_real_; DT <- NA_real_
  } else {
    T <- sum(tab$xlxmx, na.rm = TRUE) / R0
    r <- log(R0) / T
    lambda <- exp(r)
    DT <- log(2) / r
  }
  
  list(
    life = tab,
    params = c(R0 = R0, T = T, r = r, lambda = lambda, DT = DT),
    n_females = n0
  )
}

bootstrap_percentile <- function(df_clean, development_days = 0, female_ratio = 0.5, B = 2000, seed = 123) {
  set.seed(seed)
  
  ids <- unique(df_clean$female_id)
  n0 <- length(ids)
  if (n0 < 3) stop("Bootstrap requires at least 3 females.", call. = FALSE)
  if (!is.finite(B) || B < 200) stop("B must be >= 200.", call. = FALSE)
  
  mat <- matrix(NA_real_, nrow = B, ncol = 5)
  colnames(mat) <- c("R0", "T", "r", "lambda", "DT")
  
  b <- 1
  tries <- 0
  max_tries <- B * 30
  
  while (b <= B && tries < max_tries) {
    tries <- tries + 1
    ids_b <- sample(ids, size = n0, replace = TRUE)
    if (length(unique(ids_b)) < 2) next
    
    db <- do.call(rbind, lapply(ids_b, function(id) df_clean[df_clean$female_id == id, , drop = FALSE]))
    
    mat[b, ] <- tryCatch(
      as.numeric(build_life(db, development_days, female_ratio)$params),
      error = function(e) rep(NA_real_, 5)
    )
    b <- b + 1
  }
  
  if (sum(is.finite(mat[, 1])) < 50) stop("Too few valid bootstrap resamples. Increase the number of females.", call. = FALSE)
  
  ci <- apply(mat, 2, function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 50) return(c(NA_real_, NA_real_))
    stats::quantile(x, c(0.025, 0.975), na.rm = TRUE)
  })
  
  list(method = "bootstrap_percentile", ci = ci, n_females = n0)
}

jackknife_delete1 <- function(df_clean, development_days = 0, female_ratio = 0.5) {
  ids <- unique(df_clean$female_id)
  n <- length(ids)
  if (n < 3) stop("Jackknife requires at least 3 females.", call. = FALSE)
  
  full <- build_life(df_clean, development_days, female_ratio)$params
  
  loo <- matrix(NA_real_, nrow = n, ncol = 5)
  colnames(loo) <- names(full)
  
  for (i in seq_along(ids)) {
    db <- df_clean[df_clean$female_id != ids[i], , drop = FALSE]
    loo[i, ] <- tryCatch(
      as.numeric(build_life(db, development_days, female_ratio)$params),
      error = function(e) rep(NA_real_, 5)
    )
  }
  
  pseudo <- matrix(NA_real_, nrow = n, ncol = 5)
  colnames(pseudo) <- names(full)
  for (k in seq_len(ncol(pseudo))) {
    theta <- full[k]
    pseudo[, k] <- n * theta - (n - 1) * loo[, k]
  }
  
  theta_j <- colMeans(pseudo, na.rm = TRUE)
  
  se_j <- apply(pseudo, 2, function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 3) return(NA_real_)
    m <- mean(x)
    sqrt(sum((x - m)^2) / (length(x) * (length(x) - 1)))
  })
  
  ci <- rbind(theta_j - 1.96 * se_j, theta_j + 1.96 * se_j)
  colnames(ci) <- names(full)
  
  list(method = "jackknife", ci = ci, n_females = n)
}

# ---------- Example data (Plutella xylostella) ----------
# Preferred: packaged CSV at inst/extdata/plutella_xylostella_example.csv
load_plutella_example <- function() {
  p <- system.file("extdata", "plutella_xylostella_example.csv", package = "LifeTableFertility")
  if (!is.null(p) && nzchar(p) && file.exists(p)) {
    df <- utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
    names(df) <- clean_names(names(df))
    return(df)
  }
  
  # Fallback minimal example (only used if the packaged file is missing)
  data.frame(
    female_id = c(rep(1, 5), rep(2, 5), rep(3, 5)),
    age_day   = rep(1:5, 3),
    eggs      = c(12, 25, 31, 18, 5,
                  10, 22, 28, 16, 4,
                  11, 24, 30, 17, 5),
    stringsAsFactors = FALSE
  )
}

# ---------- UI ----------
ui <- shiny::fluidPage(
  shiny::titlePanel("Life Table & Fertility"),
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::fileInput("file", "Upload file", accept = c(".csv", ".txt", ".tsv")),
      shiny::actionButton("load_example", "Load example (Plutella xylostella)"),
      shiny::tags$hr(),
      shiny::numericInput("dev_days", "Immature development (days):", value = 0, min = 0, step = 1),
      shiny::numericInput("sex_ratio", "Female proportion (0–1):", value = 0.5, min = 0.01, max = 1, step = 0.01),
      shiny::tags$hr(),
      shiny::selectInput(
        "ci_method", "95% CI method:",
        choices = c(
          "None" = "none",
          "Bootstrap (percentile)" = "bootstrap",
          "Jackknife (delete-1)" = "jackknife"
        ),
        selected = "jackknife"
      ),
      shiny::conditionalPanel(
        condition = "input.ci_method == 'bootstrap'",
        shiny::numericInput("B", "Replicates (B):", value = 2000, min = 200, step = 200),
        shiny::numericInput("seed", "Seed:", value = 123, min = 1)
      ),
      shiny::tags$hr(),
      shiny::actionButton("calc", "Calculate", class = "btn-primary"),
      shiny::tags$hr(),
      shiny::downloadButton("download_table", "Download life table (CSV)"),
      shiny::downloadButton("download_params", "Download parameters (CSV)"),
      shiny::tags$div(
        style = "margin-top: 22px; border-top: 1px solid #ddd; padding-top: 12px;",
        shiny::tags$strong("How to cite"),
        shiny::tags$p(
          style = "font-size: 12.5px; line-height: 1.35;",
          "Please use: citation('LifeTableFertility')"
        )
      )
    ),
    shiny::mainPanel(
      shiny::tabsetPanel(
        shiny::tabPanel("Raw data", DT::DTOutput("raw")),
        shiny::tabPanel("Life table", DT::DTOutput("life_tbl")),
        shiny::tabPanel("Parameters", DT::DTOutput("params_tbl")),
        shiny::tabPanel(
          "Plots",
          shiny::plotOutput("p_lx", height = 260),
          shiny::downloadButton("dl_lx", "Download lx figure (PNG)"),
          shiny::tags$hr(),
          shiny::plotOutput("p_mx", height = 260),
          shiny::downloadButton("dl_mx", "Download mx figure (PNG)"),
          shiny::tags$hr(),
          shiny::plotOutput("p_lxmx", height = 260),
          shiny::downloadButton("dl_lxmx", "Download lxmx figure (PNG)")
        )
      )
    )
  )
)

# ---------- Server ----------
server <- function(input, output, session) {
  
  rv <- shiny::reactiveValues(raw = load_plutella_example(), out = NULL, ci = NULL)
  
  shiny::observeEvent(input$load_example, {
    rv$raw <- load_plutella_example()
    rv$out <- NULL
    rv$ci  <- NULL
  })
  
  shiny::observeEvent(input$file, {
    rv$raw <- read_table_smart(input$file$datapath)
    rv$out <- NULL
    rv$ci  <- NULL
  })
  
  output$raw <- DT::renderDT({
    shiny::req(rv$raw)
    DT::datatable(rv$raw, options = list(pageLength = 12))
  })
  
  shiny::observeEvent(input$calc, {
    shiny::req(rv$raw)
    tryCatch({
      req_cols(rv$raw, c("female_id", "age_day", "eggs"))
      df_clean <- clean_day1(rv$raw)
      
      rv$out <- build_life(df_clean, input$dev_days, input$sex_ratio)
      rv$ci  <- NULL
      
      if (input$ci_method == "bootstrap") {
        rv$ci <- bootstrap_percentile(df_clean, input$dev_days, input$sex_ratio, input$B, input$seed)
      } else if (input$ci_method == "jackknife") {
        rv$ci <- jackknife_delete1(df_clean, input$dev_days, input$sex_ratio)
      }
      
      shiny::showNotification("Done.", type = "message", duration = 4)
    }, error = function(e) {
      rv$out <- NULL
      rv$ci  <- NULL
      shiny::showNotification(paste("Error:", e$message), type = "error", duration = 10)
    })
  })
  
  output$life_tbl <- DT::renderDT({
    shiny::req(rv$out)
    DT::datatable(rv$out$life, options = list(pageLength = 15))
  })
  
  output$params_tbl <- DT::renderDT({
    shiny::req(rv$out)
    p <- rv$out$params
    
    ci2.5  <- rep(NA_real_, length(p))
    ci97.5 <- rep(NA_real_, length(p))
    method <- "none"
    
    if (!is.null(rv$ci)) {
      method <- rv$ci$method
      ci <- rv$ci$ci
      
      # supports both matrix orientations from bootstrap/jackknife
      if (is.matrix(ci) && nrow(ci) == 2 && all(names(p) %in% colnames(ci))) {
        ci2.5  <- as.numeric(ci[1, names(p)])
        ci97.5 <- as.numeric(ci[2, names(p)])
      } else if (is.matrix(ci) && ncol(ci) == 2 && all(names(p) %in% rownames(ci))) {
        ci2.5  <- as.numeric(ci[names(p), 1])
        ci97.5 <- as.numeric(ci[names(p), 2])
      }
    }
    
    DT::datatable(
      data.frame(
        n_females = rep(rv$out$n_females, length(p)),
        parameter = names(p),
        estimate  = as.numeric(p),
        CI2.5     = ci2.5,
        CI97.5    = ci97.5,
        CI_method = rep(method, length(p)),
        stringsAsFactors = FALSE
      ),
      options = list(pageLength = 10)
    )
  })
  
  # ---- Plots (reactive) ----
  plot_lx <- shiny::reactive({
    shiny::req(rv$out)
    d <- rv$out$life
    ggplot2::ggplot(d, ggplot2::aes(age, lx)) +
      ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::labs(title = "Survival (lx)", x = "Age (days)", y = "lx")
  })
  
  plot_mx <- shiny::reactive({
    shiny::req(rv$out)
    d <- rv$out$life
    ggplot2::ggplot(d, ggplot2::aes(age, mx)) +
      ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::labs(title = "Fertility (mx)", x = "Age (days)", y = "mx")
  })
  
  plot_lxmx <- shiny::reactive({
    shiny::req(rv$out)
    d <- rv$out$life
    ggplot2::ggplot(d, ggplot2::aes(age, lxmx)) +
      ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::labs(title = "Reproduction (lx·mx)", x = "Age (days)", y = "lx·mx")
  })
  
  output$p_lx   <- shiny::renderPlot({ plot_lx() })
  output$p_mx   <- shiny::renderPlot({ plot_mx() })
  output$p_lxmx <- shiny::renderPlot({ plot_lxmx() })
  
  # ---- Downloads (direct PNG) ----
  output$dl_lx <- shiny::downloadHandler(
    filename = function() paste0("figure_lx_", Sys.Date(), ".png"),
    contentType = "image/png",
    content = function(file) {
      shiny::req(rv$out)
      ggplot2::ggsave(filename = file, plot = plot_lx(), width = 7, height = 4.5, dpi = 300, device = "png")
    }
  )
  
  output$dl_mx <- shiny::downloadHandler(
    filename = function() paste0("figure_mx_", Sys.Date(), ".png"),
    contentType = "image/png",
    content = function(file) {
      shiny::req(rv$out)
      ggplot2::ggsave(filename = file, plot = plot_mx(), width = 7, height = 4.5, dpi = 300, device = "png")
    }
  )
  
  output$dl_lxmx <- shiny::downloadHandler(
    filename = function() paste0("figure_lxmx_", Sys.Date(), ".png"),
    contentType = "image/png",
    content = function(file) {
      shiny::req(rv$out)
      ggplot2::ggsave(filename = file, plot = plot_lxmx(), width = 7, height = 4.5, dpi = 300, device = "png")
    }
  )
  
  # ---- CSV downloads ----
  output$download_table <- shiny::downloadHandler(
    filename = function() paste0("life_table_", Sys.Date(), ".csv"),
    content = function(file) {
      shiny::req(rv$out)
      utils::write.csv(rv$out$life, file, row.names = FALSE)
    }
  )
  
  output$download_params <- shiny::downloadHandler(
    filename = function() paste0("life_params_", Sys.Date(), ".csv"),
    content = function(file) {
      shiny::req(rv$out)
      p <- rv$out$params
      
      ci2.5  <- rep(NA_real_, length(p))
      ci97.5 <- rep(NA_real_, length(p))
      method <- "none"
      
      if (!is.null(rv$ci)) {
        method <- rv$ci$method
        ci <- rv$ci$ci
        
        if (is.matrix(ci) && nrow(ci) == 2 && all(names(p) %in% colnames(ci))) {
          ci2.5  <- as.numeric(ci[1, names(p)])
          ci97.5 <- as.numeric(ci[2, names(p)])
        } else if (is.matrix(ci) && ncol(ci) == 2 && all(names(p) %in% rownames(ci))) {
          ci2.5  <- as.numeric(ci[names(p), 1])
          ci97.5 <- as.numeric(ci[names(p), 2])
        }
      }
      
      utils::write.csv(
        data.frame(
          n_females = rep(rv$out$n_females, length(p)),
          parameter = names(p),
          estimate  = as.numeric(p),
          CI2.5     = ci2.5,
          CI97.5    = ci97.5,
          CI_method = rep(method, length(p)),
          stringsAsFactors = FALSE
        ),
        file, row.names = FALSE
      )
    }
  )
}

shiny::shinyApp(ui, server)
