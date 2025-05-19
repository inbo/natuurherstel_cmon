sim_linear <- function(
  trend = 0.004, duration = 12, n_sample = 490, initial = 87.8,
  log_sd_location = 0.21, log_sd_error = 0.21
) {
  stopifnot(require("dplyr"), require("INLA"))
  n_cycle <- ceiling(duration / 10)
  if (n_cycle > 1) {
    formula <- measurement ~ year + f(
      location, model = "iid",
      hyper = list(theta = list(prior = "pc.prec", param = c(0.6, 0.01)))
    )
  } else {
    formula <- measurement ~ year
  }
  expand.grid(location = seq_len(n_sample), cycle = seq_len(n_cycle)) |>
    mutate(year = .data$location %% 10 + 10 * (.data$cycle - 1)) |>
    filter(.data$year < duration) -> design
  replicate(100, {
    design |>
      mutate(
        eta = log(initial) + .data$year * log(1 + trend) +
          rnorm(n_sample, mean = 0, sd = log_sd_location)[.data$location] +
          rnorm(n(), mean = 0, sd = log_sd_error),
        measurement = exp(.data$eta)
      ) -> dataset
    m <- inla(
      formula = formula, data = dataset, family = "lognormal",
      control.family = list(
        hyper = list(theta = list(prior = "pc.prec", param = c(0.3, 0.01)))
      )
    )
    p <- try(inla.pmarginal(0, m$marginals.fixed[["year"]]))
    if (inherits(p, "try-error")) {
      browser()
    }
    min(p, 1 - p) * 2
  })
}

display_power <- function(sims, power = 0.9, predicted) {
  stopifnot(
    require("dplyr"), require("ggplot2"), require("purrr"), require("scales")
  )
  p <- ggplot(
    sims,
    aes(x = .data$trend, ymin = .data$lcl, ymax = .data$ucl)
  ) +
    geom_hline(yintercept = power, linetype = 2) +
    geom_errorbar(aes(colour = .data$sims)) +
    geom_point(aes(y = .data$simpower, colour = .data$sims)) +
    scale_x_continuous("trend", limits = c(0, NA), labels = percent) +
    scale_y_continuous("power", limits = c(0, 1), labels = percent) +
    scale_colour_gradient(low = "red", high = "blue", limits = c(0, 1000))
  if (missing(predicted)) {
    print(p)
    return(invisible(NULL))
  }
  pred_high <- min(predicted$trend[power <= predicted$lcl])
  pred_low <- min(predicted$trend[power <= predicted$ucl])
  pred_fit <- min(predicted$trend[power <= predicted$fit])
  p <- p +
    geom_rect(
      xmin = pred_low, xmax = pred_high, ymin = -Inf, ymax = Inf, alpha = 0.05
    ) +
    geom_vline(xintercept = pred_fit, linetype = 3) +
    geom_ribbon(data = predicted, alpha = 0.2, fill = "darkgreen") +
    geom_line(data = predicted, aes(y = .data$fit), colour = "darkgreen") +
    ggtitle(
      sprintf(
        "smallest detectable trend: %.2f%% (%.2f%%; %.2f%%)",
        100 * pred_fit, 100 * pred_low, 100 * pred_high
      )
    )
  print(p)
  return(invisible(NULL))
}

predict_power <- function(sims, step_size = 0.0001) {
  stopifnot(require("dplyr"))
  sims |>
    mutate(
      not = .data$sims - .data$significant
    ) |>
    glm(formula = cbind(significant, not) ~ trend, family = binomial) -> model
  new_data <- data.frame(
    trend = seq(
      step_size, 2 * max(sims$trend[sims$simpower < 1]), by = step_size
    )
  )
  new_data |>
    predict(object = model, se.fit = TRUE) -> preds
  new_data |>
    mutate(
      fit = plogis(preds$fit),
      lcl = qnorm(0.025, preds$fit, preds$se.fit) |>
        plogis(),
      ucl = qnorm(0.975, preds$fit, preds$se.fit) |>
        plogis()
    )
}

estimate_trend <- function(
  duration = 12, n_sample = 490, initial = 87.8, log_sd_location = 0.21,
  log_sd_error = 0.21, alpha = 0.1, power = 0.9, step_size = 0.0001,
  connection = duckdb::dbConnect(
    duckdb::duckdb(), dbdir = "data/power_cmon.duckdb", read_only = FALSE
  )
) {
  stopifnot(
    require("assertthat"), require("DBI"), require("dplyr"), require("duckdb"),
    require("purrr"), require("tidyr")
  )
  assert_that(
    is.count(duration), is.count(n_sample), is.number(initial),
    is.number(log_sd_location), noNA(log_sd_location), log_sd_location > 0,
    is.number(log_sd_error), noNA(log_sd_error), log_sd_error > 0,
    is.number(alpha), noNA(alpha), alpha > 0, alpha < 1,
    is.number(power), noNA(power), power > 0, power < 1
  )
  if (!"trend" %in% dbListTables(conn = connection)) {
    data.frame(
      trend = numeric(0), duration = integer(0), n_sample = integer(0),
      initial = numeric(0), log_sd_location = numeric(0),
      log_sd_error = numeric(0), p = numeric(0)
    ) |>
      dbCreateTable(conn = connection, name = "trend")
  }
  query <- sprintf(
    "SELECT trend, SUM(p < %f) AS significant, COUNT(trend) AS sims
FROM trend
WHERE
  ABS(log_sd_location - %f) < 1e-4 AND ABS(log_sd_error - %f) < 1e-4 AND
    n_sample = %i AND duration = %i AND ABS(initial - %f) < 1e-4
GROUP BY trend
ORDER BY trend",
    alpha, log_sd_location, log_sd_error, n_sample, duration, initial
  )
  dbGetQuery(conn = connection, statement = query) |>
    mutate(simpower = .data$significant / .data$sims) -> sims
  if (nrow(sims) == 0) {
    message("initializing to 0.004")
    data.frame(
      trend = 0.004, duration = duration, n_sample = n_sample,
      initial = initial, log_sd_location = log_sd_location,
      log_sd_error = log_sd_error,
      p = sim_linear(
        trend = 0.004, duration = duration, n_sample = n_sample,
        initial = initial, log_sd_location = log_sd_location,
        log_sd_error = log_sd_error
      )
    ) |>
      dbAppendTable(conn = connection, name = "trend")
    dbGetQuery(conn = connection, statement = query) |>
      mutate(
        simpower = .data$significant / .data$sims,
        lcl = map2(.data$significant, .data$sims, ~ binom.test(.x, .y)) |>
          map(~ .x$conf.int),
        ucl = map_dbl(.data$lcl, ~ .x[2]),
        lcl = map_dbl(.data$lcl, ~ .x[1])
      ) -> sims
    display_power(sims = sims, power = power)
  }
  while (max(sims$simpower) <= power) {
    extra <- 2 * max(sims$trend)
    message("increasing to ", extra)
    data.frame(
      trend = extra, duration = duration, n_sample = n_sample,
      initial = initial, log_sd_location = log_sd_location,
      log_sd_error = log_sd_error,
      p = sim_linear(
        trend = extra, duration = duration, n_sample = n_sample,
        initial = initial, log_sd_location = log_sd_location,
        log_sd_error = log_sd_error
      )
    ) |>
      dbAppendTable(conn = connection, name = "trend")
    dbGetQuery(conn = connection, statement = query) |>
      mutate(
        simpower = .data$significant / .data$sims,
        lcl = map2(.data$significant, .data$sims, ~ binom.test(.x, .y)) |>
          map(~ .x$conf.int),
        ucl = map_dbl(.data$lcl, ~ .x[2]),
        lcl = map_dbl(.data$lcl, ~ .x[1])
      ) -> sims
    display_power(sims = sims, power = power)
  }
  while (
    min(sims$trend) >= 2 * step_size && min(sims$simpower) >= (power / 2)
  ) {
    extra <- round(min(sims$trend) / 2, -floor(log10(step_size)))
    stopifnot(extra > 0)
    message("decreasing to ", extra)
    data.frame(
      trend = extra, duration = duration, n_sample = n_sample,
      initial = initial, log_sd_location = log_sd_location,
      log_sd_error = log_sd_error,
      p = sim_linear(
        trend = extra, duration = duration, n_sample = n_sample,
        initial = initial, log_sd_location = log_sd_location,
        log_sd_error = log_sd_error
      )
    ) |>
      dbAppendTable(conn = connection, name = "trend")
    dbGetQuery(conn = connection, statement = query) |>
      mutate(
        simpower = .data$significant / .data$sims,
        lcl = map2(.data$significant, .data$sims, ~ binom.test(.x, .y)) |>
          map(~ .x$conf.int),
        ucl = map_dbl(.data$lcl, ~ .x[2]),
        lcl = map_dbl(.data$lcl, ~ .x[1])
      ) -> sims
    display_power(sims = sims, power = power)
  }
  dbGetQuery(conn = connection, statement = query) |>
    mutate(
      simpower = .data$significant / .data$sims,
      lcl = map2(.data$significant, .data$sims, ~ binom.test(.x, .y)) |>
        map(~ .x$conf.int),
      ucl = map_dbl(.data$lcl, ~ .x[2]),
      lcl = map_dbl(.data$lcl, ~ .x[1])
    ) -> sims
  stopifnot(nrow(sims) > 1)
  predicted <- predict_power(sims = sims, step_size = step_size)

  predicted |>
    filter(.data$lcl < power, power < .data$ucl) |>
    bind_rows(
      predicted |>
        filter(.data$lcl > power) |>
        slice_min(.data$trend, n = 1),
      sims |>
        filter(.data$lcl < power, power < .data$ucl, .data$sims < 1000) |>
        select("trend")
    ) |>
    left_join(
      sims |>
        select(-"lcl", -"ucl"),
      by = "trend"
    ) |>
    mutate(to_do = 1000 - replace_na(.data$sims, 0)) |>
    filter(.data$to_do > 0) -> candidate
  while (nrow(candidate) > 0) {
    display_power(sims = sims, power = power, predicted = predicted)
    if (nrow(candidate) == 1) {
      extra <- candidate$trend
    } else {
      extra <- sample(candidate$trend, size = 1, prob = candidate$to_do)
    }
    message("interpolate ", extra)
    data.frame(
      trend = extra, duration = duration, n_sample = n_sample,
      initial = initial, log_sd_location = log_sd_location,
      log_sd_error = log_sd_error,
      p = sim_linear(
        trend = extra, duration = duration, n_sample = n_sample,
        initial = initial, log_sd_location = log_sd_location,
        log_sd_error = log_sd_error
      )
    ) |>
      dbAppendTable(conn = connection, name = "trend")
    dbGetQuery(conn = connection, statement = query) |>
      mutate(
        simpower = .data$significant / .data$sims,
        lcl = map2(.data$significant, .data$sims, ~ binom.test(.x, .y)) |>
          map(~ .x$conf.int),
        ucl = map_dbl(.data$lcl, ~ .x[2]),
        lcl = map_dbl(.data$lcl, ~ .x[1])
      ) -> sims
    predicted <- predict_power(sims = sims, step_size = step_size)
    predicted |>
      filter(.data$lcl < power, power < .data$ucl) |>
      bind_rows(
        predicted |>
          filter(.data$lcl > power) |>
          slice_min(.data$trend, n = 1),
        sims |>
          filter(.data$lcl < power, power < .data$ucl, .data$sims < 1000) |>
          select("trend")
      ) |>
      left_join(
        sims |>
          select(-"lcl", -"ucl"),
        by = "trend"
      ) |>
      mutate(to_do = 1000 - replace_na(.data$sims, 0)) |>
      filter(.data$to_do > 0) -> candidate
  }
  dbGetQuery(conn = connection, statement = query) |>
    mutate(
      simpower = .data$significant / .data$sims,
      lcl = map2(.data$significant, .data$sims, ~ binom.test(.x, .y)) |>
        map(~ .x$conf.int),
      ucl = map_dbl(.data$lcl, ~ .x[2]),
      lcl = map_dbl(.data$lcl, ~ .x[1])
    ) -> sims
  predicted <- predict_power(sims = sims, step_size = step_size)
  display_power(sims = sims, power = power, predicted = predicted)
  return(
    c(
      estimate = min(predicted$trend[predicted$fit >= power]),
      lcl = min(predicted$trend[predicted$ucl >= power]),
      ucl = min(predicted$trend[predicted$lcl >= power])
    )
  )
}

# var: akker: 0.085
# n_sample: 794
estimate_trend(
  duration = 6, log_sd_location = 0.21, log_sd_error = 0.21, n_sample = 794
)
estimate_trend(
  duration = 12, log_sd_location = 0.21, log_sd_error = 0.21, n_sample = 794
)
estimate_trend(
  duration = 18, log_sd_location = 0.21, log_sd_error = 0.21, n_sample = 794
)
estimate_trend(
  duration = 20, log_sd_location = 0.21, log_sd_error = 0.21, n_sample = 794
)
estimate_trend(
  duration = 24, log_sd_location = 0.21, log_sd_error = 0.21, n_sample = 794
)
estimate_trend(
  duration = 30, log_sd_location = 0.21, log_sd_error = 0.21, n_sample = 794
)

# var bos: 0.17
# n_sample: 490
estimate_trend(
  duration = 6, log_sd_location = 0.36, log_sd_error = 0.21, n_sample = 490
)
estimate_trend(
  duration = 12, log_sd_location = 0.36, log_sd_error = 0.21, n_sample = 490
)
estimate_trend(
  duration = 18, log_sd_location = 0.36, log_sd_error = 0.21, n_sample = 490
)
estimate_trend(
  duration = 20, log_sd_location = 0.36, log_sd_error = 0.21, n_sample = 490
)
estimate_trend(
  duration = 24, log_sd_location = 0.36, log_sd_error = 0.21, n_sample = 490
)
estimate_trend(
  duration = 30, log_sd_location = 0.36, log_sd_error = 0.21, n_sample = 490
)

# var natuur: 0.38
# n_sample: 446
estimate_trend(
  duration = 6, log_sd_location = 0.58, log_sd_error = 0.21, n_sample = 446
)
estimate_trend(
  duration = 12, log_sd_location = 0.58, log_sd_error = 0.21, n_sample = 446
)
estimate_trend(
  duration = 18, log_sd_location = 0.58, log_sd_error = 0.21, n_sample = 446
)
estimate_trend(
  duration = 20, log_sd_location = 0.58, log_sd_error = 0.21, n_sample = 446
)
estimate_trend(
  duration = 24, log_sd_location = 0.58, log_sd_error = 0.21, n_sample = 446
)
estimate_trend(
  duration = 30, log_sd_location = 0.58, log_sd_error = 0.21, n_sample = 446
)

# var: grasland: 0.085
# n_sample: 406
estimate_trend(
  duration = 6, log_sd_location = 0.21, log_sd_error = 0.21, n_sample = 406
)
estimate_trend(
  duration = 12, log_sd_location = 0.21, log_sd_error = 0.21, n_sample = 406
)
estimate_trend(
  duration = 18, log_sd_location = 0.21, log_sd_error = 0.21, n_sample = 406
)
estimate_trend(
  duration = 20, log_sd_location = 0.21, log_sd_error = 0.21, n_sample = 406
)
estimate_trend(
  duration = 24, log_sd_location = 0.21, log_sd_error = 0.21, n_sample = 406
)
estimate_trend(
  duration = 30, log_sd_location = 0.21, log_sd_error = 0.21, n_sample = 406
)

# var ruimtebeslag: 0.14
# n_sample: 458
estimate_trend(
  duration = 6, log_sd_location = 0.32, log_sd_error = 0.21, n_sample = 458
)
estimate_trend(
  duration = 12, log_sd_location = 0.32, log_sd_error = 0.21, n_sample = 458
)
estimate_trend(
  duration = 18, log_sd_location = 0.32, log_sd_error = 0.21, n_sample = 458
)
estimate_trend(
  duration = 20, log_sd_location = 0.32, log_sd_error = 0.21, n_sample = 458
)
estimate_trend(
  duration = 24, log_sd_location = 0.32, log_sd_error = 0.21, n_sample = 458
)
estimate_trend(
  duration = 30, log_sd_location = 0.32, log_sd_error = 0.21, n_sample = 458
)
