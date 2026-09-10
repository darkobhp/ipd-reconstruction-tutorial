# Scenario 1: no numbers-at-risk table
# Generated from the Quarto tutorial. Change only values in USER INPUTS.

# ---- packages ----

required_packages <- c(
  "dplyr",
  "purrr",
  "ggplot2",
  "survival",
  "survminer",
  "broom",
  "knitr",
  "gtsummary"
)

packages_to_install <- setdiff(
  required_packages,
  rownames(installed.packages())
)

if (length(packages_to_install) > 0) {
  install.packages(packages_to_install)
}

invisible(
  lapply(required_packages, library, character.only = TRUE)
)

# ---- user-inputs ----

# ====================== USER INPUTS: CHANGE THESE ======================

# STEP 1: Keep TRUE while learning with the dummy data.
# Change this to FALSE when you are ready to use your own two CSV files.
use_demo_data <- TRUE

# STEP 2: Enter the information for Arm 1.
# Arm 1 becomes the reference group in the hazard-ratio calculation.
arm1_file <- "~/Downloads/arm1_digitized_curve.csv" # CHANGE: path to Arm 1 CSV
arm1_N <- 62L                                        # CHANGE: number at risk at time 0
arm1_label <- "Control"                             # CHANGE: arm name from the paper

# Leave as NA_real_ when the CSV time column is already in months/years.
# If the digitized time is scaled from 0 to 1, replace NA_real_ with the
# maximum time shown on the paper's x-axis (for example, 36).
arm1_time_max <- NA_real_                            # USUALLY LEAVE UNCHANGED

# STEP 3: Enter the information for Arm 2.
arm2_file <- "~/Downloads/arm2_digitized_curve.csv" # CHANGE: path to Arm 2 CSV
arm2_N <- 80L                                        # CHANGE: number at risk at time 0
arm2_label <- "Treatment"                           # CHANGE: arm name from the paper
arm2_time_max <- NA_real_                            # USUALLY LEAVE UNCHANGED

# STEP 4: Enter the results reported in the paper for the sanity checks.
# The HR must be Arm 2 versus Arm 1. If the paper gives Arm 1 versus Arm 2,
# enter its reciprocal (1 / the reported HR) instead.
reported_hr_arm2_vs_arm1 <- 0.66       # CHANGE: published HR, Arm 2 vs Arm 1
reported_median_os_arm1 <- 16.75       # CHANGE: published Arm 1 median OS, months
reported_median_os_arm2 <- 22.93       # CHANGE: published Arm 2 median OS, months

# STEP 5: Review the plot settings.
time_unit <- "Months"
x_axis_max <- 36                      # CHANGE to the largest x-axis time to display
x_axis_break_by <- 6                  # CHANGE to the spacing used in the paper (spacing between ticks on the x-axis)
arm_colours <- c("#D55E00", "#0072B2") # OPTIONAL: Arm 1 and Arm 2 colours

# ==================== END OF VALUES TO CHANGE =========================

# ---- demo-data ----

if (use_demo_data) {
  arm1_curve <- data.frame(
    time = c(0, 3, 6, 9, 12, 15, 18, 21, 24, 30, 36),
    survival = c(1.00, 0.94, 0.84, 0.74, 0.64, 0.55, 0.47, 0.39, 0.32, 0.23, 0.16)
  )

  arm2_curve <- data.frame(
    time = c(0, 3, 6, 9, 12, 15, 18, 21, 24, 30, 36),
    survival = c(1.00, 0.97, 0.91, 0.85, 0.78, 0.70, 0.63, 0.56, 0.49, 0.38, 0.29)
  )
}

# ---- clean-km-function ----

clamp01 <- function(x) {
  pmin(pmax(x, 0), 1)
}

clean_km <- function(t_raw, S_raw) {
  if (length(t_raw) != length(S_raw)) {
    stop("Time and survival vectors must have the same length.")
  }

  curve <- data.frame(
    time = as.numeric(t_raw),
    survival = as.numeric(S_raw)
  ) |>
    dplyr::filter(!is.na(time), !is.na(survival)) |>
    dplyr::arrange(time) |>
    # At duplicate times, keep the lowest survival value (the full step drop).
    dplyr::group_by(time) |>
    dplyr::summarise(survival = min(survival), .groups = "drop")

  if (nrow(curve) < 2) {
    stop("The digitized curve must contain at least two valid points.")
  }

  if (curve$time[1] < 0) {
    stop("Time values cannot be negative.")
  }

  if (curve$time[1] > 0) {
    curve <- dplyr::bind_rows(
      data.frame(time = 0, survival = 1),
      curve
    )
  } else {
    curve$survival[1] <- 1
  }

  curve$survival <- clamp01(curve$survival)
  curve$survival <- cummin(curve$survival)
  curve
}

# ---- estimate-events-function ----

events_no_censor <- function(S, N) {
  if (length(S) < 2) {
    stop("At least two survival values are required.")
  }
  if (length(N) != 1 || is.na(N) || N <= 0 || N != as.integer(N)) {
    stop("Baseline N must be one positive integer.")
  }

  n_points <- length(S)
  at_risk <- integer(n_points)
  events <- integer(n_points - 1L)
  at_risk[1] <- as.integer(N)

  for (j in seq_len(n_points - 1L)) {
    if (S[j] <= 0 || at_risk[j] == 0) {
      estimated_events <- at_risk[j]
    } else {
      estimated_events <- round(
        at_risk[j] * (1 - S[j + 1L] / S[j])
      )
    }

    estimated_events <- max(
      0L,
      min(as.integer(estimated_events), at_risk[j])
    )

    events[j] <- estimated_events
    at_risk[j + 1L] <- at_risk[j] - events[j]
  }

  # Reconcile the total with the survival probability at the last point.
  target_survivors <- as.integer(round(N * S[n_points]))
  target_total_events <- as.integer(N - target_survivors)
  event_difference <- target_total_events - sum(events)

  if (event_difference > 0L) {
    for (j in rev(seq_along(events))) {
      available <- N - sum(events)
      add <- min(event_difference, available)
      events[j] <- events[j] + add
      event_difference <- event_difference - add
      if (event_difference == 0L) break
    }
  } else if (event_difference < 0L) {
    events_to_remove <- -event_difference
    for (j in rev(seq_along(events))) {
      remove <- min(events_to_remove, events[j])
      events[j] <- events[j] - remove
      events_to_remove <- events_to_remove - remove
      if (events_to_remove == 0L) break
    }
  }

  events
}

# ---- build-ipd-function ----

build_ipd_uniform <- function(curve, events, N, arm_label) {
  if (length(events) != nrow(curve) - 1L) {
    stop("There must be one event count per curve interval.")
  }

  event_times <- purrr::map2(
    seq_len(nrow(curve) - 1L),
    events,
    function(j, number_events) {
      if (number_events == 0L) return(numeric(0))

      seq(
        from = curve$time[j],
        to = curve$time[j + 1L],
        length.out = number_events + 2L
      )[-c(1L, number_events + 2L)]
    }
  ) |>
    unlist(use.names = FALSE)

  number_censored <- N - length(event_times)
  if (number_censored < 0) {
    stop("Estimated events exceed baseline N.")
  }

  data.frame(
    time = c(event_times, rep(max(curve$time), number_censored)),
    status = c(rep(1L, length(event_times)), rep(0L, number_censored)),
    arm = arm_label
  )
}

# ---- reconstruct-arm-function ----

read_digitized_curve <- function(file) {
  first_attempt <- tryCatch(
    read.csv(file, header = FALSE),
    error = function(e) NULL
  )

  if (!is.null(first_attempt) && ncol(first_attempt) >= 2 &&
      all(!is.na(suppressWarnings(as.numeric(first_attempt[[1]]))))) {
    return(first_attempt[, 1:2])
  }

  second_attempt <- read.csv(file, header = TRUE)
  if (ncol(second_attempt) < 2) {
    stop("The digitized file must contain at least two columns.")
  }
  second_attempt[, 1:2]
}

reconstruct_arm <- function(file = NULL,
                            curve_data = NULL,
                            baseline_N,
                            arm_label,
                            time_max = NA_real_) {
  if (is.null(curve_data)) {
    if (is.null(file) || !file.exists(path.expand(file))) {
      stop("Digitized file not found for arm: ", arm_label)
    }
    raw_curve <- read_digitized_curve(path.expand(file))
  } else {
    raw_curve <- curve_data[, 1:2]
  }

  raw_time <- suppressWarnings(as.numeric(raw_curve[[1]]))
  raw_survival <- suppressWarnings(as.numeric(raw_curve[[2]]))

  if (!is.na(time_max)) {
    raw_time <- raw_time * time_max
  }

  cleaned_curve <- clean_km(raw_time, raw_survival)
  estimated_events <- events_no_censor(
    cleaned_curve$survival,
    baseline_N
  )
  reconstructed_ipd <- build_ipd_uniform(
    cleaned_curve,
    estimated_events,
    baseline_N,
    arm_label
  )

  list(
    curve = cleaned_curve,
    interval_events = estimated_events,
    ipd = reconstructed_ipd
  )
}

# ---- reconstruct-two-arms ----

arm1_result <- reconstruct_arm(
  file = if (!use_demo_data) arm1_file else NULL,
  curve_data = if (use_demo_data) arm1_curve else NULL,
  baseline_N = arm1_N,
  arm_label = arm1_label,
  time_max = arm1_time_max
)

arm2_result <- reconstruct_arm(
  file = if (!use_demo_data) arm2_file else NULL,
  curve_data = if (use_demo_data) arm2_curve else NULL,
  baseline_N = arm2_N,
  arm_label = arm2_label,
  time_max = arm2_time_max
)

combined_ipd <- dplyr::bind_rows(
  arm1_result$ipd,
  arm2_result$ipd
) |>
  dplyr::mutate(
    status = as.integer(status),
    arm = factor(arm, levels = c(arm1_label, arm2_label))
  )

head(combined_ipd)

# ---- validation-tests ----

stopifnot(
  nrow(arm1_result$ipd) == arm1_N,
  nrow(arm2_result$ipd) == arm2_N,
  nrow(combined_ipd) == arm1_N + arm2_N,
  all(combined_ipd$status %in% c(0L, 1L)),
  all(is.finite(combined_ipd$time)),
  all(combined_ipd$time >= 0),
  !anyNA(combined_ipd)
)

validation_summary <- combined_ipd |>
  dplyr::group_by(arm) |>
  dplyr::summarise(
    participants = dplyr::n(),
    deaths = sum(status),
    censored = sum(status == 0L),
    event_rate_percent = 100 * mean(status),
    .groups = "drop"
  )

validation_summary |>
  dplyr::mutate(
    event_rate_percent = sprintf("%.2f%%", event_rate_percent)
  )

# ---- arm-diagnostic-function ----

plot_arm_check <- function(result, arm_label, colour) {
  reconstructed_fit <- survival::survfit(
    survival::Surv(time, status) ~ 1,
    data = result$ipd
  )

  reconstructed_curve <- data.frame(
    time = c(0, reconstructed_fit$time),
    survival = c(1, reconstructed_fit$surv),
    source = "Reconstructed IPD"
  )

  ggplot2::ggplot() +
    ggplot2::geom_step(
      data = reconstructed_curve,
      ggplot2::aes(time, survival, colour = source),
      linewidth = 1
    ) +
    ggplot2::geom_point(
      data = result$curve,
      ggplot2::aes(time, survival, shape = "Digitized points"),
      colour = colour,
      size = 2
    ) +
    ggplot2::scale_colour_manual(values = c("Reconstructed IPD" = colour)) +
    ggplot2::scale_shape_manual(values = c("Digitized points" = 16)) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      labels = scales::label_percent(accuracy = 1)
    ) +
    ggplot2::labs(
      title = paste("Reconstruction check:", arm_label),
      x = paste("Time (", time_unit, ")", sep = ""),
      y = "Survival probability",
      colour = NULL,
      shape = NULL
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}

# ---- plot-arm-1 ----

plot_arm_check(arm1_result, arm1_label, arm_colours[1])

# ---- plot-arm-2 ----

plot_arm_check(arm2_result, arm2_label, arm_colours[2])

# ---- survival-models ----

survival_object <- survival::Surv(
  time = combined_ipd$time,
  event = combined_ipd$status
)

km_fit <- survival::survfit(
  survival_object ~ arm,
  data = combined_ipd
)

logrank_fit <- survival::survdiff(
  survival_object ~ arm,
  data = combined_ipd
)

logrank_p <- stats::pchisq(
  logrank_fit$chisq,
  df = length(logrank_fit$n) - 1L,
  lower.tail = FALSE
)

cox_fit <- survival::coxph(
  survival_object ~ arm,
  data = combined_ipd
)

cox_result <- broom::tidy(
  cox_fit,
  exponentiate = TRUE,
  conf.int = TRUE
)

hr <- cox_result$estimate[1]
hr_ci_low <- cox_result$conf.low[1]
hr_ci_high <- cox_result$conf.high[1]
hr_p <- cox_result$p.value[1]

# ---- summary-statistics ----

format_p <- function(p) {
  ifelse(
    is.na(p),
    "NA",
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

median_table <- survminer::surv_median(km_fit) |>
  dplyr::transmute(
    arm = sub("^arm=", "", strata),
    median_os = median
  )

event_table <- combined_ipd |>
  dplyr::group_by(arm) |>
  dplyr::summarise(
    participants = dplyr::n(),
    deaths = sum(status),
    event_rate = 100 * mean(status),
    .groups = "drop"
  ) |>
  dplyr::mutate(arm = as.character(arm))

arm_summary <- event_table |>
  dplyr::left_join(median_table, by = "arm") |>
  dplyr::mutate(
    median_os_display = ifelse(
      is.na(median_os),
      "NR",
      sprintf("%.2f", median_os)
    ),
    event_rate_display = sprintf("%.2f%%", event_rate)
  )

model_summary <- data.frame(
  comparison = paste(arm2_label, "vs", arm1_label),
  hazard_ratio = sprintf("%.2f", hr),
  confidence_interval_95 = sprintf("%.2f to %.2f", hr_ci_low, hr_ci_high),
  cox_p_value = format_p(hr_p),
  logrank_p_value = format_p(logrank_p)
)

arm_summary |>
  dplyr::select(
    arm,
    participants,
    deaths,
    median_os = median_os_display,
    mortality_event_rate = event_rate_display
  )

model_summary

# ---- final-plot ----

arm_annotation_lines <- paste0(
  arm_summary$arm,
  ": median OS = ",
  arm_summary$median_os_display,
  " ",
  time_unit,
  "; event rate = ",
  arm_summary$event_rate_display
)

annotation_text <- paste(
  sprintf(
    "HR, %s vs %s: %.2f (95%% CI %.2f–%.2f); Cox p%s",
    arm2_label,
    arm1_label,
    hr,
    hr_ci_low,
    hr_ci_high,
    ifelse(hr_p < 0.001, "<0.001", paste0("=", sprintf("%.3f", hr_p)))
  ),
  paste0(
    "Log-rank p",
    ifelse(
      logrank_p < 0.001,
      "<0.001",
      paste0("=", sprintf("%.3f", logrank_p))
    )
  ),
  arm_annotation_lines[1],
  arm_annotation_lines[2],
  sep = "\n"
)

final_km_plot <- survminer::ggsurvplot(
  fit = km_fit,
  data = combined_ipd,
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.25,
  conf.int = FALSE,
  censor = TRUE,
  palette = arm_colours,
  legend.title = "Arm",
  legend.labs = c(arm1_label, arm2_label),
  xlab = paste0("Time (", time_unit, ")"),
  ylab = "Overall survival probability",
  xlim = c(0, x_axis_max),
  break.time.by = x_axis_break_by,
  surv.scale = "percent",
  ggtheme = ggplot2::theme_classic(base_size = 13),
  tables.theme = survminer::theme_cleantable(base_size = 11)
)

final_km_plot$plot <- final_km_plot$plot +
  ggplot2::annotate(
    "label",
    x = 0.98 * x_axis_max,
    y = 0.98,
    label = annotation_text,
    hjust = 1,
    vjust = 1,
    size = 3.5,
    linewidth = 0.25,
    fill = scales::alpha("white", 0.88)
  ) +
  ggplot2::theme(
    legend.position = "top",
    legend.text = ggplot2::element_text(size = 11),
    legend.title = ggplot2::element_text(size = 11, face = "bold"),
    axis.title = ggplot2::element_text(face = "bold")
  )

print(final_km_plot)

# ---- sanity-check-function ----

run_reconstruction_sanity_checks <- function(
    reconstructed_hr,
    reported_hr,
    reconstructed_medians,
    reported_medians,
    arm_labels,
    median_tolerance_months = 0.20,
    log_hr_tolerance = 0.10) {

  if (length(reconstructed_hr) != 1L || !is.finite(reconstructed_hr) ||
      reconstructed_hr <= 0 || length(reported_hr) != 1L ||
      !is.finite(reported_hr) || reported_hr <= 0) {
    stop("Reported and reconstructed HRs must each be one positive number.")
  }

  if (length(reconstructed_medians) != 2L ||
      length(reported_medians) != 2L || length(arm_labels) != 2L) {
    stop("Supply exactly two reconstructed medians, reported medians, and arm labels.")
  }

  median_difference <- abs(reconstructed_medians - reported_medians)
  median_pass <- is.finite(median_difference) &
    median_difference <= median_tolerance_months + sqrt(.Machine$double.eps)

  log_hr_difference <- abs(log(reported_hr / reconstructed_hr))
  hr_pass <- is.finite(log_hr_difference) &&
    log_hr_difference <= log_hr_tolerance + sqrt(.Machine$double.eps)

  results <- data.frame(
    check = c(
      paste0("Median OS: ", arm_labels[1]),
      paste0("Median OS: ", arm_labels[2]),
      paste0("HR: ", arm_labels[2], " vs ", arm_labels[1])
    ),
    reported = sprintf("%.2f", c(reported_medians, reported_hr)),
    reconstructed = sprintf("%.2f", c(reconstructed_medians, reconstructed_hr)),
    comparison = c(
      sprintf("absolute difference = %.2f months", median_difference),
      sprintf("absolute log ratio = %.3f", log_hr_difference)
    ),
    acceptance_rule = c(
      rep(sprintf("difference <= %.2f months", median_tolerance_months), 2L),
      sprintf("absolute log ratio <= %.2f", log_hr_tolerance)
    ),
    result = ifelse(c(median_pass, hr_pass), "PASS", "FAIL"),
    check.names = FALSE
  )

  display_data <- data.frame(
    check = results$check,
    reported = c(reported_medians, reported_hr),
    reconstructed = c(reconstructed_medians, reconstructed_hr),
    diagnostic = c(median_difference, log_hr_difference),
    acceptance_limit = c(
      median_tolerance_months,
      median_tolerance_months,
      log_hr_tolerance
    ),
    passed = as.integer(c(median_pass, hr_pass)),
    check.names = FALSE
  )

  failed_checks <- results$check[results$result == "FAIL"]
  if (length(failed_checks) == 0L) {
    cat("\nOVERALL RESULT: PASS - all median OS and HR sanity checks passed.\n")
  } else {
    cat("\nOVERALL RESULT: FAIL - the following check(s) did not pass:\n")
    cat(paste0("- ", failed_checks, collapse = "\n"), "\n")
  }

  invisible(
    list(
      results = results,
      display_data = display_data,
      all_pass = length(failed_checks) == 0L
    )
  )
}

display_sanity_check_table <- function(checks, table_title) {
  overall_result <- ifelse(checks$all_pass, "PASS", "FAIL")
  full_title <- paste0(table_title, " - Overall result: ", overall_result)

  table_data <- checks$display_data |>
    dplyr::mutate(
      check = factor(check, levels = check)
    )

  formatted_table <- table_data |>
    gtsummary::tbl_summary(
      by = check,
      include = -check,
      type = dplyr::everything() ~ "continuous",
      statistic = dplyr::everything() ~ "{mean}",
      digits = list(
        reported ~ 2,
        reconstructed ~ 2,
        diagnostic ~ 3,
        acceptance_limit ~ 2,
        passed ~ 0
      ),
      label = list(
        reported ~ "Reported",
        reconstructed ~ "Reconstructed",
        diagnostic ~ "Absolute difference / log ratio",
        acceptance_limit ~ "Acceptance limit (<=)",
        passed ~ "Result"
      ),
      missing = "no"
    ) |>
    gtsummary::modify_header(
      gtsummary::all_stat_cols() ~ "**{level}**"
    ) |>
    gtsummary::modify_footnote_header(
      footnote = NA_character_,
      columns = gtsummary::all_stat_cols()
    ) |>
    gtsummary::modify_footnote_header(
      footnote = paste0(
        "Cross-check and confirm that the original and reconstructed median ",
        "time to survival have been reported in the same unit of months."
      ),
      columns = c(stat_1, stat_2),
      replace = TRUE
    ) |>
    gtsummary::modify_table_body(
      ~ .x |>
        dplyr::mutate(
          dplyr::across(
            dplyr::starts_with("stat_"),
            ~ ifelse(
              variable == "passed",
              ifelse(.x == "1", "PASS", "FAIL"),
              .x
            )
          )
        )
    ) |>
    gtsummary::bold_labels() |>
    gtsummary::modify_caption(paste0("**", full_title, "**"))

  # Quarto renders the gtsummary object as a publication-ready table.
  if (knitr::is_html_output()) {
    return(formatted_table)
  }

  # Printing a gtsummary table interactively opens its HTML output in the
  # RStudio Viewer without calling the native spreadsheet Data Viewer.
  if (interactive()) {
    tryCatch(
      print(formatted_table),
      error = function(error_condition) {
        message("The HTML Viewer could not be opened; using the console table.")
      }
    )
  }

  # Always retain a console fallback and keep the underlying data accessible.
  print(checks$results, row.names = FALSE)
  invisible(checks$results)
}

# ---- scenario-1-sanity-check ----

scenario1_medians <- stats::setNames(
  arm_summary$median_os,
  arm_summary$arm
)

scenario1_sanity_checks <- run_reconstruction_sanity_checks(
  reconstructed_hr = hr,
  reported_hr = reported_hr_arm2_vs_arm1,
  reconstructed_medians = c(
    scenario1_medians[[arm1_label]],
    scenario1_medians[[arm2_label]]
  ),
  reported_medians = c(
    reported_median_os_arm1,
    reported_median_os_arm2
  ),
  arm_labels = c(arm1_label, arm2_label)
)

display_sanity_check_table(
  scenario1_sanity_checks,
  "Scenario 1 reconstruction sanity checks"
)

