build_life <- function(df_clean, development_days = 0, female_ratio = 0.5) {
  if (nrow(df_clean) < 5) stop("Too few valid rows after cleaning.", call. = FALSE)
  if (any(df_clean$age_day < 1, na.rm = TRUE)) stop("age_day must start at 1.", call. = FALSE)
  if (any(df_clean$eggs < 0, na.rm = TRUE)) stop("eggs cannot be negative.", call. = FALSE)
  if (!is.finite(development_days) || development_days < 0) stop("Immature development must be >= 0.", call. = FALSE)
  if (!is.finite(female_ratio) || female_ratio <= 0 || female_ratio > 1) stop("Female proportion must be in (0, 1].", call. = FALSE)
  
  ids <- unique(df_clean$female_id)
  n0 <- length(ids)
  if (n0 < 2) {
    stop(paste0("At least 2 females are required. After cleaning, n_females = ", n0, "."), call. = FALSE)
  }
  
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
  
  b <- 1; tries <- 0; max_tries <- B * 30
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
  
  if (sum(is.finite(mat[, 1])) < 50) {
    stop("Too few valid bootstrap resamples. Increase the number of females.", call. = FALSE)
  }
  
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
