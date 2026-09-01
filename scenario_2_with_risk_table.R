# SCENARIO 2: TWO-ARM IPD RECONSTRUCTION WITH A NUMBERS-AT-RISK TABLE
# Generated from the accompanying Quarto tutorial.
# Start in the USER INPUTS section. Change only the clearly marked values.
# Then run this file from top to bottom.

# ---- PACKAGES ----

required_packages <- c(
  "dplyr",
  "purrr",
  "ggplot2",
  "survival",
  "survminer",
  "broom"
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

# ---- RISK TABLE INPUTS ----

# ====================== USER INPUTS: CHANGE THESE ======================

# STEP 1: Keep TRUE while learning with the dummy data.
# Change this to FALSE when using your own digitized curves and risk tables.
risk_use_demo_data <- TRUE

# STEP 2: Enter Arm 1 information. Arm 1 is the reference group.
risk_arm1_file <- "~/Downloads/arm1_digitized_curve.csv" # CHANGE: Arm 1 CSV
risk_arm1_label <- "Control"                             # CHANGE: paper's arm name

# STEP 3: Copy the risk-table time headings and Arm 1 counts from the paper.
# The first time must be 0, and the first count is the baseline number at risk.
risk_arm1_times <- c(0, 6, 12, 18, 24, 30) # CHANGE: risk-table times
risk_arm1_n <- c(62, 52, 39, 28, 20, 14)   # CHANGE: Arm 1 risk counts

# CHANGE only if the paper reports a total number of events for Arm 1.
# Otherwise leave as NA_integer_.
risk_arm1_total_events <- NA_integer_

# STEP 4: Enter Arm 2 information.
risk_arm2_file <- "~/Downloads/arm2_digitized_curve.csv" # CHANGE: Arm 2 CSV
risk_arm2_label <- "Treatment"                           # CHANGE: paper's arm name

# STEP 5: Copy the risk-table time headings and Arm 2 counts from the paper.
# These vectors must have the same length. Usually both arms share risk times,
# but enter them separately in case the paper reports different intervals.
risk_arm2_times <- c(0, 6, 12, 18, 24, 30) # CHANGE: risk-table times
risk_arm2_n <- c(68, 61, 52, 42, 34, 26)   # CHANGE: Arm 2 risk counts
risk_arm2_total_events <- NA_integer_       # CHANGE only if reported

# STEP 6: Review time and plot settings.
risk_arm1_time_max <- NA_real_ # USUALLY LEAVE UNCHANGED
risk_arm2_time_max <- NA_real_ # USUALLY LEAVE UNCHANGED
time_unit <- "Months"
x_axis_max <- 36                 # CHANGE to the largest time to display
x_axis_break_by <- 6             # CHANGE to the spacing used in the paper (spacing between ticks on the x-axis)
arm_colours <- c("#D55E00", "#0072B2") # OPTIONAL

# ==================== END OF VALUES TO CHANGE =========================

# ---- RISK TABLE DEMO CURVES ----

risk_demo_control <- data.frame(
  time = c(0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36),
  survival = c(1.00, 0.94, 0.84, 0.74, 0.64, 0.55, 0.47,
               0.39, 0.32, 0.27, 0.23, 0.19, 0.16)
)

risk_demo_treatment <- data.frame(
  time = c(0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36),
  survival = c(1.00, 0.97, 0.91, 0.85, 0.78, 0.70, 0.63,
               0.56, 0.49, 0.43, 0.38, 0.33, 0.29)
)

# ---- CLEAN KM FUNCTION ----

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

# ---- RECONSTRUCT ARM FUNCTION ----

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

# ---- ARM DIAGNOSTIC FUNCTION ----

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

# ---- RISK TABLE RECONSTRUCTION FUNCTION ----

reconstruct_arm_with_risk_table <- function(
    file = NULL,
    curve_data = NULL,
    risk_times,
    numbers_at_risk,
    arm_label,
    total_events = NA_integer_,
    time_max = NA_real_,
    max_iterations = 1000L) {

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
  if (!is.na(time_max)) raw_time <- raw_time * time_max

  curve <- clean_km(raw_time, raw_survival)
  tS <- curve$time
  S <- curve$survival

  risk_times <- as.numeric(risk_times)
  numbers_at_risk <- as.integer(numbers_at_risk)

  if (length(risk_times) != length(numbers_at_risk)) {
    stop("risk_times and numbers_at_risk must have equal lengths.")
  }
  if (length(risk_times) < 2L) {
    stop("At least two reported risk times are required for this workflow.")
  }
  if (risk_times[1] != 0 || any(diff(risk_times) <= 0)) {
    stop("Risk times must begin at 0 and be strictly increasing.")
  }
  if (any(numbers_at_risk < 0) || any(diff(numbers_at_risk) > 0)) {
    stop("Numbers at risk must be nonnegative and non-increasing.")
  }
  if (max(risk_times) >= max(tS)) {
    stop("The last risk-table time must be earlier than the last digitized time.")
  }

  lower <- vapply(
    risk_times,
    function(x) min(which(tS >= x)),
    integer(1)
  )
  upper <- vapply(
    c(risk_times[-1], Inf),
    function(x) max(which(tS < x)),
    integer(1)
  )

  n_intervals <- length(numbers_at_risk)
  n_t <- upper[n_intervals]
  baseline_n <- numbers_at_risk[1]
  effective_risk <- numbers_at_risk

  n_hat <- rep(baseline_n + 1L, n_t + 1L)
  censored <- integer(n_t)
  events <- integer(n_t)
  km_hat <- rep(1, n_t)
  last_event_index <- rep(1L, n_intervals)
  number_censored <- integer(n_intervals)

  allocate_censoring <- function(start_index, end_time_index, n_censor) {
    allocation <- integer(n_t)
    if (n_censor <= 0L || end_time_index <= start_index) return(allocation)

    censor_times <- tS[start_index] +
      seq_len(n_censor) *
      (tS[end_time_index] - tS[start_index]) / (n_censor + 1L)

    allocation[start_index:(end_time_index - 1L)] <- hist(
      censor_times,
      breaks = tS[start_index:end_time_index],
      plot = FALSE,
      include.lowest = TRUE
    )$counts
    allocation
  }

  run_interval <- function(start_index, end_index, risk_at_start,
                           previous_event_index, censor_allocation,
                           n_hat_in, events_in, km_hat_in) {
    n_local <- n_hat_in
    d_local <- events_in
    km_local <- km_hat_in
    n_local[start_index] <- risk_at_start
    last <- previous_event_index

    for (k in start_index:end_index) {
      if (start_index == 1L && k == 1L) {
        d_local[k] <- 0L
        km_local[k] <- 1
      } else if (n_local[k] <= 0 || km_local[last] <= 0) {
        d_local[k] <- 0L
        km_local[k] <- km_local[last]
      } else {
        proposed <- round(n_local[k] * (1 - S[k] / km_local[last]))
        d_local[k] <- max(0L, min(as.integer(proposed), n_local[k]))
        km_local[k] <- km_local[last] * (1 - d_local[k] / n_local[k])
      }

      n_local[k + 1L] <- n_local[k] - d_local[k] - censor_allocation[k]
      if (n_local[k + 1L] < 0L) {
        censor_allocation[k] <- max(0L, n_local[k] - d_local[k])
        n_local[k + 1L] <- 0L
      }
      if (d_local[k] > 0L) last <- k
    }

    list(
      n_hat = n_local,
      events = d_local,
      km_hat = km_local,
      censoring = censor_allocation,
      last = last
    )
  }

  # Reconcile all intervals that end at another published risk count.
  for (i in seq_len(n_intervals - 1L)) {
    start_index <- lower[i]
    end_index <- upper[i]
    next_start <- lower[i + 1L]

    initial_censor <- round(
      effective_risk[i] * S[next_start] / S[start_index] -
        effective_risk[i + 1L]
    )
    number_censored[i] <- max(0L, as.integer(initial_censor))

    converged <- FALSE
    for (iteration in seq_len(max_iterations)) {
      interval_censoring <- allocate_censoring(
        start_index,
        next_start,
        number_censored[i]
      )
      censored[start_index:end_index] <-
        interval_censoring[start_index:end_index]

      result <- run_interval(
        start_index,
        end_index,
        effective_risk[i],
        last_event_index[i],
        censored,
        n_hat,
        events,
        km_hat
      )

      n_hat <- result$n_hat
      events <- result$events
      km_hat <- result$km_hat
      censored <- result$censoring

      difference <- n_hat[next_start] - effective_risk[i + 1L]
      if (difference == 0L) {
        converged <- TRUE
        break
      }

      revised_censor <- number_censored[i] + difference
      if (revised_censor < 0L && number_censored[i] == 0L) break
      number_censored[i] <- max(0L, revised_censor)
    }

    if (!converged && n_hat[next_start] < effective_risk[i + 1L]) {
      warning(
        "Published risk count was not attainable in interval ", i,
        "; using the largest internally consistent count."
      )
      effective_risk[i + 1L] <- n_hat[next_start]
    }
    last_event_index[i + 1L] <- result$last
  }

  # Estimate censoring after the final reported risk time.
  previous_span <- tS[upper[n_intervals - 1L]] - tS[lower[1L]]
  final_span <- tS[upper[n_intervals]] - tS[lower[n_intervals]]
  initial_final_censor <- if (previous_span > 0) {
    round(sum(number_censored[seq_len(n_intervals - 1L)]) *
            final_span / previous_span)
  } else {
    0L
  }
  initial_final_censor <- min(
    max(0L, as.integer(initial_final_censor)),
    effective_risk[n_intervals]
  )

  base_n_hat <- n_hat
  base_events <- events
  base_km_hat <- km_hat
  base_censored <- censored

  evaluate_final_interval <- function(final_censor_count) {
    final_allocation <- allocate_censoring(
      lower[n_intervals],
      upper[n_intervals],
      final_censor_count
    )
    censor_candidate <- base_censored
    censor_candidate[lower[n_intervals]:upper[n_intervals]] <-
      final_allocation[lower[n_intervals]:upper[n_intervals]]

    run_interval(
      lower[n_intervals],
      upper[n_intervals],
      effective_risk[n_intervals],
      last_event_index[n_intervals],
      censor_candidate,
      base_n_hat,
      base_events,
      base_km_hat
    )
  }

  if (is.na(total_events)) {
    final_result <- evaluate_final_interval(initial_final_censor)
    number_censored[n_intervals] <- initial_final_censor
  } else {
    candidates <- 0:effective_risk[n_intervals]
    candidate_results <- lapply(candidates, evaluate_final_interval)
    candidate_event_totals <- vapply(
      candidate_results,
      function(x) sum(x$events),
      integer(1)
    )
    best <- which.min(abs(candidate_event_totals - total_events))
    final_result <- candidate_results[[best]]
    number_censored[n_intervals] <- candidates[best]
    if (candidate_event_totals[best] != total_events) {
      warning("The reported total event count could only be approximated.")
    }
  }

  n_hat <- final_result$n_hat
  events <- final_result$events
  km_hat <- final_result$km_hat
  censored <- final_result$censoring

  event_times <- rep(tS[seq_len(n_t)], events[seq_len(n_t)])
  censor_times <- unlist(
    lapply(seq_len(n_t - 1L), function(j) {
      rep((tS[j] + tS[j + 1L]) / 2, censored[j])
    }),
    use.names = FALSE
  )

  assigned <- length(event_times) + length(censor_times)
  if (assigned > baseline_n) {
    stop("Events plus censoring exceed the baseline sample size.")
  }

  remaining <- baseline_n - assigned
  ipd <- data.frame(
    time = c(event_times, censor_times, rep(tS[n_t], remaining)),
    status = c(
      rep(1L, length(event_times)),
      rep(0L, length(censor_times) + remaining)
    ),
    arm = arm_label
  )

  reconstructed_fit <- survival::survfit(
    survival::Surv(time, status) ~ 1,
    data = ipd
  )
  reconstructed_risk <- summary(
    reconstructed_fit,
    times = risk_times,
    extend = TRUE
  )$n.risk

  list(
    curve = curve,
    ipd = ipd,
    events_by_point = events[seq_len(n_t)],
    censored_by_point = censored[seq_len(n_t)],
    estimated_risk = data.frame(
      time = risk_times,
      published = numbers_at_risk,
      used_by_algorithm = effective_risk,
      reconstructed = reconstructed_risk
    )
  )
}

# ---- RUN RISK TABLE RECONSTRUCTION ----

risk_arm1_result <- reconstruct_arm_with_risk_table(
  file = if (!risk_use_demo_data) risk_arm1_file else NULL,
  curve_data = if (risk_use_demo_data) risk_demo_control else NULL,
  risk_times = risk_arm1_times,
  numbers_at_risk = risk_arm1_n,
  arm_label = risk_arm1_label,
  total_events = risk_arm1_total_events,
  time_max = risk_arm1_time_max
)

risk_arm2_result <- reconstruct_arm_with_risk_table(
  file = if (!risk_use_demo_data) risk_arm2_file else NULL,
  curve_data = if (risk_use_demo_data) risk_demo_treatment else NULL,
  risk_times = risk_arm2_times,
  numbers_at_risk = risk_arm2_n,
  arm_label = risk_arm2_label,
  total_events = risk_arm2_total_events,
  time_max = risk_arm2_time_max
)

combined_ipd_with_risk <- dplyr::bind_rows(
  risk_arm1_result$ipd,
  risk_arm2_result$ipd
) |>
  dplyr::mutate(
    status = as.integer(status),
    arm = factor(arm, levels = c(risk_arm1_label, risk_arm2_label))
  )

stopifnot(
  nrow(risk_arm1_result$ipd) == risk_arm1_n[1],
  nrow(risk_arm2_result$ipd) == risk_arm2_n[1],
  all(combined_ipd_with_risk$status %in% c(0L, 1L)),
  !anyNA(combined_ipd_with_risk)
)

head(combined_ipd_with_risk)

# ---- RISK COUNT VALIDATION ----

dplyr::bind_rows(
  dplyr::mutate(risk_arm1_result$estimated_risk, arm = risk_arm1_label),
  dplyr::mutate(risk_arm2_result$estimated_risk, arm = risk_arm2_label)
) |>
  dplyr::select(arm, time, published, used_by_algorithm, reconstructed)

# ---- RISK ARM CHECKS ----

plot_arm_check(risk_arm1_result, risk_arm1_label, arm_colours[1])
plot_arm_check(risk_arm2_result, risk_arm2_label, arm_colours[2])

# ---- RISK TABLE ANALYSIS ----

risk_km_fit <- survival::survfit(
  survival::Surv(time, status) ~ arm,
  data = combined_ipd_with_risk
)

risk_logrank <- survival::survdiff(
  survival::Surv(time, status) ~ arm,
  data = combined_ipd_with_risk
)
risk_logrank_p <- stats::pchisq(
  risk_logrank$chisq,
  df = length(risk_logrank$n) - 1L,
  lower.tail = FALSE
)

risk_cox <- survival::coxph(
  survival::Surv(time, status) ~ arm,
  data = combined_ipd_with_risk
)
risk_cox_result <- broom::tidy(
  risk_cox,
  exponentiate = TRUE,
  conf.int = TRUE
)

risk_hr <- risk_cox_result$estimate[1]
risk_hr_low <- risk_cox_result$conf.low[1]
risk_hr_high <- risk_cox_result$conf.high[1]
risk_hr_p <- risk_cox_result$p.value[1]

risk_medians <- survminer::surv_median(risk_km_fit) |>
  dplyr::transmute(
    arm = sub("^arm=", "", strata),
    median_os = median
  )

risk_arm_summary <- combined_ipd_with_risk |>
  dplyr::group_by(arm) |>
  dplyr::summarise(
    participants = dplyr::n(),
    deaths = sum(status),
    event_rate = 100 * mean(status),
    .groups = "drop"
  ) |>
  dplyr::mutate(arm = as.character(arm)) |>
  dplyr::left_join(risk_medians, by = "arm") |>
  dplyr::mutate(
    median_display = ifelse(
      is.na(median_os), "NR", sprintf("%.2f", median_os)
    ),
    event_rate_display = sprintf("%.2f%%", event_rate)
  )

risk_arm_summary |>
  dplyr::select(
    arm,
    participants,
    deaths,
    median_os = median_display,
    mortality_event_rate = event_rate_display
  )

# ---- RISK TABLE FINAL PLOT ----

risk_annotation <- paste(
  sprintf(
    "HR, %s vs %s: %.2f (95%% CI %.2f–%.2f); Cox p%s",
    risk_arm2_label,
    risk_arm1_label,
    risk_hr,
    risk_hr_low,
    risk_hr_high,
    ifelse(
      risk_hr_p < 0.001,
      "<0.001",
      paste0("=", sprintf("%.3f", risk_hr_p))
    )
  ),
  paste0(
    "Log-rank p",
    ifelse(
      risk_logrank_p < 0.001,
      "<0.001",
      paste0("=", sprintf("%.3f", risk_logrank_p))
    )
  ),
  paste0(
    risk_arm_summary$arm,
    ": median OS = ", risk_arm_summary$median_display, " ", time_unit,
    "; event rate = ", risk_arm_summary$event_rate_display,
    collapse = "\n"
  ),
  sep = "\n"
)

risk_final_plot <- survminer::ggsurvplot(
  risk_km_fit,
  data = combined_ipd_with_risk,
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.25,
  conf.int = FALSE,
  palette = arm_colours,
  legend.title = "Arm",
  legend.labs = c(risk_arm1_label, risk_arm2_label),
  xlab = paste0("Time (", time_unit, ")"),
  ylab = "Overall survival probability",
  xlim = c(0, x_axis_max),
  break.time.by = x_axis_break_by,
  surv.scale = "percent",
  ggtheme = ggplot2::theme_classic(base_size = 13),
  tables.theme = survminer::theme_cleantable(base_size = 11)
)

risk_final_plot$plot <- risk_final_plot$plot +
  ggplot2::annotate(
    "label",
    x = 0.98 * x_axis_max,
    y = 0.98,
    label = risk_annotation,
    hjust = 1,
    vjust = 1,
    size = 3.5,
    linewidth = 0.25,
    fill = scales::alpha("white", 0.88)
  ) +
  ggplot2::theme(
    legend.position = "top",
    axis.title = ggplot2::element_text(face = "bold")
  )

print(risk_final_plot)

