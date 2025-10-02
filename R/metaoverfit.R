#' metaoverfit: Optimism-Corrected Heterogeneity in Meta-Regression
#'
#' Tools for detecting and correcting overfitting in meta-regression analyses
#' through cross-validation and bootstrap methods.
#'
#' @keywords internal
"_PACKAGE"

# Imports needed by this package (generate proper NAMESPACE)
#' @importFrom stats coef cor model.matrix quantile sd
#' @importFrom utils txtProgressBar setTxtProgressBar
NULL

# Silence NOTES about non-standard evaluation column names used in ggplot2
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("Type", "R2het", "Predicted", "Observed"))
}

# ==============================
# Core functions
# ==============================

#' Calculate Apparent R-squared for Heterogeneity
#'
#' @param yi Vector of effect sizes.
#' @param vi Vector of sampling variances.
#' @param mods Model matrix or formula for moderators. If a matrix, include intercept if desired.
#' @param data Optional data frame (required if `mods` is a formula).
#' @param method Heterogeneity estimator (default: "REML").
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{r2het} Apparent R^2_het.
#'   \item \code{r2het_adj} Adjusted R^2_het with small-sample penalty.
#'   \item \code{tau2_null}, \code{tau2_full} Estimated tau^2 from null and full models.
#'   \item \code{I2} I^2 from the null model.
#'   \item \code{k}, \code{p} Number of studies and parameters.
#'   \item \code{model_null}, \code{model_full} \code{metafor::rma} fits.
#' }
#' @export
#'
#' @examples
#' k <- 30
#' yi <- rnorm(k, 0, sqrt(0.1))
#' vi <- runif(k, 0.01, 0.1)
#' mods <- cbind(1, rnorm(k))
#' r2het(yi, vi, mods)
r2het <- function(yi, vi, mods = NULL, data = NULL, method = "REML") {
  if (length(yi) != length(vi)) {
    stop("yi and vi must have the same length")
  }

  # Fit null model
  fit_null <- metafor::rma(yi = yi, vi = vi, data = data, method = method)
  tau2_null <- fit_null$tau2
  I2 <- fit_null$I2

  # If no moderators, return trivial
  if (is.null(mods)) {
    return(list(
      r2het = 0,
      r2het_adj = 0,
      tau2_null = tau2_null,
      tau2_full = tau2_null,
      I2 = I2,
      k = length(yi),
      p = 1,
      model_null = fit_null,
      model_full = fit_null
    ))
  }

  # Fit full model (formula vs matrix)
  if (inherits(mods, "formula")) {
    fit_full <- metafor::rma(yi = yi, vi = vi, mods = mods, data = data, method = method)
    p <- length(coef(fit_full))
  } else {
    mods <- as.matrix(mods)
    fit_full <- metafor::rma(yi = yi, vi = vi, mods = ~ mods - 1, data = data, method = method)
    p <- ncol(mods)
  }

  tau2_full <- fit_full$tau2

  # Apparent R^2_het
  r2het_val <- if (tau2_null > 0) max(0, 1 - tau2_full / tau2_null) else 0

  # Adjusted R^2_het
  k <- length(yi)
  r2het_adj <- if (tau2_null > 0 && k > p) max(0, 1 - (tau2_full / tau2_null) * ((k - 1) / (k - p))) else 0

  list(
    r2het = r2het_val,
    r2het_adj = r2het_adj,
    tau2_null = tau2_null,
    tau2_full = tau2_full,
    I2 = I2,
    k = k,
    p = p,
    model_null = fit_null,
    model_full = fit_full
  )
}

#' Cross-Validated R-squared for Heterogeneity
#'
#' @param yi Vector of effect sizes.
#' @param vi Vector of sampling variances.
#' @param mods Model matrix or formula for moderators.
#' @param data Optional data frame.
#' @param method Heterogeneity estimator (default: "REML").
#' @param cv_method "loo" for leave-one-out or "kfold" for k-fold CV.
#' @param k_folds Number of folds for k-fold CV (default: 5).
#' @param verbose Show progress bar.
#'
#' @return List containing apparent/corrected R^2_het, optimism, diagnostics, and predictions.
#' @export
#'
#' @examples
#' k <- 30
#' yi <- rnorm(k, 0, sqrt(0.1))
#' vi <- runif(k, 0.01, 0.1)
#' mods <- cbind(1, rnorm(k))
#' result <- r2het_cv(yi, vi, mods)
#' round(result$optimism * 100, 1)
r2het_cv <- function(yi, vi, mods = NULL, data = NULL, method = "REML",
                     cv_method = "loo", k_folds = 5, verbose = TRUE) {
  k <- length(yi)

  # Apparent
  apparent <- r2het(yi, vi, mods, data, method)

  if (is.null(mods)) {
    return(list(
      r2het_apparent = 0,
      r2het_corrected = 0,
      optimism = 0,
      tau2_null = apparent$tau2_null,
      tau2_cv = apparent$tau2_null,
      convergence_rate = 100,
      rmse = NA_real_,
      mae = NA_real_,
      correlation = NA_real_,
      predictions = rep(NA_real_, k),
      yi = yi  # include observed for plotting completeness
    ))
  }

  # Build moderator matrix
  if (inherits(mods, "formula")) {
    if (is.null(data)) stop("Data must be provided when using a formula for 'mods'.")
    mods_matrix <- model.matrix(mods, data)
  } else {
    mods_matrix <- as.matrix(mods)
  }
  p <- ncol(mods_matrix)

  if (k <= p + 1) {
    warning("Too few studies for stable cross-validation (k <= p+1)")
    return(list(
      r2het_apparent = apparent$r2het,
      r2het_corrected = 0,
      optimism = apparent$r2het,
      tau2_null = apparent$tau2_null,
      tau2_cv = NA_real_,
      convergence_rate = 0,
      rmse = NA_real_,
      mae = NA_real_,
      correlation = NA_real_,
      predictions = rep(NA_real_, k),
      yi = yi
    ))
  }

  # CV predictions
  yhat_cv <- rep(NA_real_, k)
  converged <- rep(FALSE, k)

  if (cv_method == "loo") {
    if (verbose) pb <- txtProgressBar(min = 0, max = k, style = 3)
    for (i in seq_len(k)) {
      tryCatch({
        fit_cv <- metafor::rma(
          yi = yi[-i], vi = vi[-i],
          mods = ~ mods_matrix[-i, , drop = FALSE] - 1,
          method = method
        )
        yhat_cv[i] <- sum(mods_matrix[i, ] * coef(fit_cv))
        converged[i] <- TRUE
      }, error = function(e) {
        yhat_cv[i] <- NA_real_
      })
      if (verbose) setTxtProgressBar(pb, i)
    }
    if (verbose) close(pb)
  } else if (cv_method == "kfold") {
    folds <- sample(rep(seq_len(k_folds), length.out = k))
    for (fold in seq_len(k_folds)) {
      test_idx <- which(folds == fold)
      train_idx <- which(folds != fold)
      tryCatch({
        fit_cv <- metafor::rma(
          yi = yi[train_idx], vi = vi[train_idx],
          mods = ~ mods_matrix[train_idx, , drop = FALSE] - 1,
          method = method
        )
        for (j in test_idx) {
          yhat_cv[j] <- sum(mods_matrix[j, ] * coef(fit_cv))
          converged[j] <- TRUE
        }
      }, error = function(e) {})
    }
  } else {
    stop("cv_method must be 'loo' or 'kfold'.")
  }

  valid <- converged & !is.na(yhat_cv)

  if (sum(valid) > k / 2) {
    tau2_cv <- mean(pmax(0, (yi[valid] - yhat_cv[valid])^2 - vi[valid]))
    r2het_cv_val <- if (apparent$tau2_null > 0) max(0, 1 - tau2_cv / apparent$tau2_null) else 0
    rmse <- sqrt(mean((yi[valid] - yhat_cv[valid])^2))
    mae <- mean(abs(yi[valid] - yhat_cv[valid]))
    correlation <- cor(yi[valid], yhat_cv[valid])
  } else {
    tau2_cv <- NA_real_
    r2het_cv_val <- NA_real_
    rmse <- NA_real_
    mae <- NA_real_
    correlation <- NA_real_
  }

  optimism <- if (!is.na(r2het_cv_val)) apparent$r2het - r2het_cv_val else NA_real_

  list(
    r2het_apparent = apparent$r2het,
    r2het_corrected = r2het_cv_val,
    optimism = optimism,
    tau2_null = apparent$tau2_null,
    tau2_cv = tau2_cv,
    convergence_rate = sum(valid) / k * 100,
    rmse = rmse,
    mae = mae,
    correlation = correlation,
    predictions = if (sum(valid) > 0) yhat_cv else NULL,
    yi = yi
  )
}

#' Bootstrap Confidence Intervals for R-squared
#'
#' @param yi Vector of effect sizes.
#' @param vi Vector of sampling variances.
#' @param mods Model matrix or formula for moderators.
#' @param data Optional data frame.
#' @param method Heterogeneity estimator.
#' @param B Number of bootstrap samples (default: 1000).
#' @param conf_level Confidence level (default: 0.95).
#' @param verbose Show progress.
#'
#' @return List with confidence intervals and diagnostics.
#' @export
r2het_boot <- function(yi, vi, mods = NULL, data = NULL, method = "REML",
                       B = 1000, conf_level = 0.95, verbose = TRUE) {
  k <- length(yi)
  alpha <- 1 - conf_level

  r2het_ap_vec <- rep(NA_real_, B)
  r2het_cv_vec <- rep(NA_real_, B)

  if (verbose) pb <- txtProgressBar(min = 0, max = B, style = 3)

  for (b in seq_len(B)) {
    idx <- sample(seq_len(k), k, replace = TRUE)
    yi_b <- yi[idx]
    vi_b <- vi[idx]

    # Prepare bootstrapped mods/data
    mods_b <- NULL
    data_b <- NULL
    has_matrix_mods <- !is.null(mods) && inherits(mods, "matrix")
    has_formula_mods <- !is.null(mods) && inherits(mods, "formula")

    if (has_matrix_mods) {
      mods_b <- mods[idx, , drop = FALSE]
    } else if (has_formula_mods) {
      if (is.null(data)) stop("When 'mods' is a formula, 'data' must be provided.")
      data_b <- data[idx, , drop = FALSE]
    }

    # Compute CV in bootstrap sample
    tryCatch({
      if (is.null(mods)) {
        res <- r2het_cv(yi_b, vi_b, mods = NULL, data = NULL, method = method,
                        cv_method = "loo", verbose = FALSE)
      } else if (has_formula_mods) {
        res <- r2het_cv(yi_b, vi_b, mods, data_b, method, cv_method = "loo", verbose = FALSE)
      } else {
        res <- r2het_cv(yi_b, vi_b, mods_b, NULL, method, cv_method = "loo", verbose = FALSE)
      }
      r2het_ap_vec[b] <- res$r2het_apparent
      r2het_cv_vec[b] <- res$r2het_corrected
    }, error = function(e) {
      # leave NA on failure
    })

    if (verbose) setTxtProgressBar(pb, b)
  }

  if (verbose) close(pb)

  r2het_ap_vec <- r2het_ap_vec[!is.na(r2het_ap_vec)]
  r2het_cv_vec <- r2het_cv_vec[!is.na(r2het_cv_vec)]

  if (length(r2het_ap_vec) > B / 2) {
    ci_apparent <- quantile(r2het_ap_vec, c(alpha / 2, 1 - alpha / 2))
    mean_apparent <- mean(r2het_ap_vec)
    se_apparent <- sd(r2het_ap_vec)
  } else {
    ci_apparent <- c(NA_real_, NA_real_)
    mean_apparent <- NA_real_
    se_apparent <- NA_real_
  }

  if (length(r2het_cv_vec) > B / 2) {
    ci_corrected <- quantile(r2het_cv_vec, c(alpha / 2, 1 - alpha / 2))
    mean_corrected <- mean(r2het_cv_vec)
    se_corrected <- sd(r2het_cv_vec)
    prop_at_zero <- mean(r2het_cv_vec < 0.01)
  } else {
    ci_corrected <- c(NA_real_, NA_real_)
    mean_corrected <- NA_real_
    se_corrected <- NA_real_
    prop_at_zero <- NA_real_
  }

  list(
    ci_apparent = ci_apparent,
    ci_corrected = ci_corrected,
    mean_apparent = mean_apparent,
    mean_corrected = mean_corrected,
    se_apparent = se_apparent,
    se_corrected = se_corrected,
    prop_at_zero = prop_at_zero,
    convergence = (length(r2het_ap_vec) / B) * 100,
    bootstrap_values = list(apparent = r2het_ap_vec, corrected = r2het_cv_vec)
  )
}

#' Check for Overfitting in Meta-Regression
#'
#' @param yi Vector of effect sizes.
#' @param vi Vector of sampling variances.
#' @param mods Model matrix or formula for moderators.
#' @param data Optional data frame.
#' @param method Heterogeneity estimator.
#' @param B Number of bootstrap samples.
#'
#' @return A list (class "metaoverfit") summarizing overfitting risk and metrics.
#' @export
check_overfitting <- function(yi, vi, mods = NULL, data = NULL,
                              method = "REML", B = 500) {
  k <- length(yi)

  if (is.null(mods)) {
    rpt <- list(
      k = k,
      p = 1,
      k_per_p = k,
      risk_category = "None",
      expected_optimism = "<10%",
      actual_optimism = 0,
      r2het_apparent = 0,
      r2het_corrected = 0,
      ci_corrected = c(NA_real_, NA_real_),
      recommendation = "No moderators - no overfitting risk",
      convergence_rate = 100
    )
    class(rpt) <- c("metaoverfit", "list")
    return(rpt)
  }

  # Determine p
  if (inherits(mods, "matrix")) {
    p <- ncol(mods)
  } else if (inherits(mods, "formula")) {
    tmp <- metafor::rma(yi = yi, vi = vi, mods = mods, data = data, method = method)
    p <- length(coef(tmp))
  } else {
    stop("'mods' must be a matrix or a formula, or NULL.")
  }

  k_per_p <- k / p

  # Heuristic risk bands
  if (k < 20 || k_per_p < 5) {
    risk_category <- "Extreme"
    expected_optimism <- ">40%"
    recommendation <- "DO NOT conduct meta-regression - sample size too small"
  } else if (k_per_p < 10) {
    risk_category <- "Severe"
    expected_optimism <- "20-40%"
    recommendation <- "Results highly unreliable - consider exploratory only"
  } else if (k_per_p < 15) {
    risk_category <- "Moderate"
    expected_optimism <- "10-20%"
    recommendation <- "Interpret with caution - report optimism correction"
  } else {
    risk_category <- "Low"
    expected_optimism <- "<10%"
    recommendation <- "Acceptable sample size - still report optimism correction"
  }

  cv_result <- r2het_cv(yi, vi, mods, data, method, verbose = FALSE)
  boot_result <- r2het_boot(yi, vi, mods, data, method, B = B, verbose = FALSE)

  rpt <- list(
    k = k,
    p = p,
    k_per_p = round(k_per_p, 1),
    risk_category = risk_category,
    expected_optimism = expected_optimism,
    actual_optimism = round(cv_result$optimism * 100, 1),
    r2het_apparent = round(cv_result$r2het_apparent * 100, 1),
    r2het_corrected = round(cv_result$r2het_corrected * 100, 1),
    ci_corrected = round(boot_result$ci_corrected * 100, 1),
    recommendation = recommendation,
    convergence_rate = cv_result$convergence_rate
  )
  class(rpt) <- c("metaoverfit", "list")
  rpt
}

#' Sample Size Recommendation for Meta-Regression
#'
#' @param p Number of parameters (including intercept).
#' @param target_optimism Maximum acceptable optimism (default: 0.10).
#'
#' @return Invisibly returns the minimum required k; prints a small report.
#' @export
sample_size_recommendation <- function(p, target_optimism = 0.10) {
  if (target_optimism <= 0.05) {
    min_ratio <- 20
  } else if (target_optimism <= 0.10) {
    min_ratio <- 15
  } else if (target_optimism <= 0.20) {
    min_ratio <- 10
  } else {
    min_ratio <- 5
  }

  min_k <- max(20, ceiling(p * min_ratio))

  cat("Sample Size Recommendation\n")
  cat(strrep("-", 40), "\n")
  cat("Parameters (p):", p, "\n")
  cat("Target optimism:", target_optimism * 100, "%\n")
  cat("Minimum k/p ratio:", min_ratio, "\n")
  cat("MINIMUM k required:", min_k, "\n")
  cat("\nNote: This ensures optimism <", target_optimism * 100, "%\n")

  invisible(min_k)
}

#' Print Method for metaoverfit Objects
#'
#' @param x metaoverfit object.
#' @param ... Additional arguments (unused).
#'
#' @export
print.metaoverfit <- function(x, ...) {
  cat("\n", strrep("=", 60), "\n")
  cat("META-REGRESSION OVERFITTING ASSESSMENT\n")
  cat(strrep("=", 60), "\n\n")

  cat("Sample Size:\n")
  cat("  k =", x$k, "studies\n")
  cat("  p =", x$p, "parameters\n")
  cat("  k/p ratio =", x$k_per_p, "\n\n")

  cat("Risk Assessment:\n")
  cat("  Category:", x$risk_category, "\n")
  cat("  Expected optimism:", x$expected_optimism, "\n")
  cat("  Actual optimism:", x$actual_optimism, "%\n\n")

  cat("R^2 for Heterogeneity:\n")
  cat("  Apparent R^2_het:", x$r2het_apparent, "%\n")
  cat("  Corrected R^2_het:", x$r2het_corrected, "%\n")
  cat("  95% CI:", paste0("[", x$ci_corrected[1], "%, ",
                          x$ci_corrected[2], "%]"), "\n\n")

  cat("RECOMMENDATION:\n")
  cat("  ", x$recommendation, "\n")
  cat("\n", strrep("=", 60), "\n")
}

#' Plot Overfitting Diagnostics
#'
#' @param result Output from \code{check_overfitting()} or \code{r2het_cv()}.
#' @param type "bar" for comparison of apparent vs corrected R^2_het,
#'   or "scatter" for observed vs predicted (requires predictions).
#'
#' @return A ggplot object.
#' @export
plot_overfitting <- function(result, type = "bar") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for plotting")
  }

  if (type == "bar") {
    df <- data.frame(
      Type  = c("Apparent", "Corrected"),
      R2het = c(result$r2het_apparent * 100,
                result$r2het_corrected * 100)
    )

    p <- ggplot2::ggplot(df, ggplot2::aes(x = Type, y = R2het, fill = Type)) +
      ggplot2::geom_bar(stat = "identity", width = 0.6) +
      ggplot2::scale_fill_manual(values = c("Apparent" = "#E74C3C",
                                            "Corrected" = "#27AE60")) +
      ggplot2::labs(
        x = "", y = "R^2_het (%)",
        title = "Overfitting in Meta-Regression",
        subtitle = paste(
          "Optimism:",
          round((result$r2het_apparent - result$r2het_corrected) * 100, 1), "%"
        )
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none")

  } else if (type == "scatter") {
    if (is.null(result$predictions) || all(is.na(result$predictions))) {
      stop("No predictions available in 'result' to draw a scatter plot.")
    }
    if (is.null(result$yi)) {
      stop("Observed values 'yi' not found in 'result' (call r2het_cv which returns 'yi').")
    }

    df <- data.frame(
      Observed  = result$yi,
      Predicted = result$predictions
    )

    p <- ggplot2::ggplot(df, ggplot2::aes(x = Predicted, y = Observed)) +
      ggplot2::geom_point(size = 3, alpha = 0.7) +
      ggplot2::geom_abline(slope = 1, intercept = 0,
                           linetype = "dashed", color = "red") +
      ggplot2::labs(
        x = "Predicted Effect Size",
        y = "Observed Effect Size",
        title = "Cross-Validation Performance",
        subtitle = paste("Correlation:", round(result$correlation, 3))
      ) +
      ggplot2::theme_minimal()
  } else {
    stop("type must be 'bar' or 'scatter'.")
  }

  p
}
