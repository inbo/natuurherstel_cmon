sim_linear <- function(
  trend = 0.004, duration = 12, n_sample = 490, initial = 87.8,
  log_sd_location = 0.21, log_sd_error = 0.21
) {
  stopifnot(require(glmmTMB), require(tidyverse))
  n_cycle <- ceiling(duration / 10)
  if (n_cycle > 1) {
    formula <- measurement ~ year + (1 | location)
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
    m <- try(
      glmmTMB(
        formula = formula, family = lognormal(link = "log"), data = dataset
      )
    )
    if (inherits(m, "try-error") || any(is.na(m$sdr$cov.fixed))) {
      return(numeric(0))
    }
    if (!inherits(m, "glmmTMB")) {
      stop(class(m))
    }
    coef(summary(m))$cond["year", "Pr(>|z|)"]
  }) |>
    unlist() |>
    na.omit()
}

sim_change <- function(
  trend = 0.004, n_period = 2, duration = 6, n_sample = 490, initial = 87.8,
  log_sd_location = 0.21, log_sd_error = 0.21
) {
  stopifnot(require(glmmTMB), require(tidyverse), n_period >= 2)
  n_cycle <- ceiling(duration * n_period / 10)
  if (n_cycle > 1) {
    formula <- measurement ~ period + (1 | location)
  } else {
    formula <- measurement ~ period
  }
  expand.grid(location = seq_len(n_sample), cycle = seq_len(n_cycle)) |>
    mutate(
      year = .data$location %% 10 + 10 * (.data$cycle - 1),
      period = factor(.data$year %/% duration)
    ) |>
    filter(.data$year < duration * n_period) -> design
  replicate(100, {
    design |>
      mutate(
        eta = log(initial) + .data$year * log(1 + trend) +
          rnorm(n_sample, mean = 0, sd = log_sd_location)[.data$location] +
          rnorm(n(), mean = 0, sd = log_sd_error),
        measurement = exp(.data$eta)
      ) -> dataset
    m <- try(
      glmmTMB(
        formula = formula, family = lognormal(link = "log"), data = dataset
      )
    )
    if (inherits(m, "try-error") || any(is.na(m$sdr$cov.fixed))) {
      return(numeric(0))
    }
    if (!inherits(m, "glmmTMB")) {
      stop(class(m))
    }
    tail(coef(summary(m))$cond[, "Pr(>|z|)"], 1)
  }) |>
    unlist() |>
    na.omit()
}

estimate_trend <- function(
  duration = 12, n_sample = 490, initial = 87.8, log_sd_location = 0.21,
  log_sd_error = 0.21,  alpha = 0.1, power = 0.9,
  connection = duckdb::dbConnect(
    duckdb::duckdb(), dbdir = "data/power_cmon.duckdb", read_only = FALSE
  )
) {
  stopifnot(require(assertthat), require(duckdb), require(tidyverse))
  assert_that(
    is.count(duration), is.count(n_sample), is.number(initial),
    is.number(log_sd_location), noNA(log_sd_location), log_sd_location > 0,
    is.number(log_sd_error), noNA(log_sd_error), log_sd_error > 0,
    is.number(alpha), noNA(alpha), alpha > 0, alpha < 1,
    is.number(power), noNA(power), power > 0, power < 1
  )
  if (!"trend" %in% dbListTables(conn = connection)) {
    #   dbRemoveTable(conn = connection, name = "trend")
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
  sims <- dbGetQuery(conn = connection, statement = query)
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
  }
  dbGetQuery(conn = connection, statement = query) |>
    mutate(simpower = .data$significant / .data$sims) -> sims
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
      mutate(simpower = .data$significant / .data$sims) -> sims
  }
  while (min(sims$simpower) >= power) {
    extra <- round(min(sims$trend) / 2, 4)
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
      mutate(simpower = .data$significant / .data$sims) -> sims
  }
  stopifnot(nrow(sims) > 1)
  sims |>
    mutate(
      not = .data$sims - .data$significant
    ) |>
    glm(formula = cbind(significant, not) ~ trend, family = binomial) -> model
  new_data <- data.frame(
    trend = seq(0.0001, 2 * max(sims$trend[sims$simpower < 1]), by = 0.0001)
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
    ) -> candidate
  candidate |>
    filter(.data$lcl < power, power < .data$ucl) |>
    bind_rows(
      candidate |>
        filter(.data$lcl > power) |>
        slice_min(.data$trend, n = 1)
    ) |>
    left_join(sims, by = "trend") |>
    mutate(to_do = 1000 - replace_na(.data$sims, 0)) |>
    filter(.data$to_do > 0) -> candidate
  while (nrow(candidate) > 0) {
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
      mutate(simpower = .data$significant / .data$sims) -> sims
    sims |>
      mutate(
        not = .data$sims - .data$significant
      ) |>
      glm(formula = cbind(significant, not) ~ trend, family = binomial) -> model
    new_data <- data.frame(
      trend = seq(0.0001, 2 * max(sims$trend[sims$simpower < 1]), by = 0.0001)
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
      ) -> candidate
    p <- sims |>
      filter(.data$simpower < 1) |>
      mutate(
        lcl = map2(.data$significant, .data$sims, ~ binom.test(.x, .y)) |>
          map(~ .x$conf.int),
        ucl = map_dbl(.data$lcl, ~ .x[2]),
        lcl = map_dbl(.data$lcl, ~ .x[1])
      ) |>
      ggplot(aes(x = trend, ymin = lcl, ymax = ucl)) +
      geom_hline(yintercept = power, linetype = 2) +
      geom_rect(
        xmin = min(candidate$trend[candidate$ucl >= power]),
        xmax = min(candidate$trend[candidate$lcl >= power]),
        ymin = -Inf, ymax = Inf, alpha = 0.1
      ) +
      geom_vline(
        xintercept = min(candidate$trend[candidate$fit >= power]),
        linetype = 2
      ) +
      geom_errorbar() +
      geom_point(aes(size = sims, y = simpower), alpha = 0.2) +
      geom_ribbon(data = candidate, alpha = 0.1, fill = "blue") +
      geom_line(data = candidate, aes(y = fit), colour = "blue") +
      scale_y_continuous("power", limits = c(0, 1)) +
      scale_size_continuous(limits = c(0, NA))
    print(p)
    candidate |>
      filter(.data$lcl < power, power < .data$ucl) |>
      bind_rows(
        candidate |>
          filter(.data$lcl > power) |>
          slice_min(.data$trend, n = 1)
      ) |>
      left_join(sims, by = "trend") |>
      mutate(to_do = 1000 - replace_na(.data$sims, 0)) |>
      filter(.data$to_do > 0) -> candidate
  }
  sims |>
    mutate(
      not = .data$sims - .data$significant
    ) |>
    glm(formula = cbind(significant, not) ~ trend, family = binomial) -> model
  new_data <- data.frame(
    trend = seq(0.0001, max(sims$trend[sims$simpower < 1]), by = 0.0001)
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
    ) -> new_data
  p <- sims |>
    filter(.data$simpower < 1) |>
    mutate(
      lcl = map2(.data$significant, .data$sims, ~ binom.test(.x, .y)) |>
        map(~ .x$conf.int),
      ucl = map_dbl(.data$lcl, ~ .x[2]),
      lcl = map_dbl(.data$lcl, ~ .x[1])
    ) |>
    ggplot(aes(x = trend, ymin = lcl, ymax = ucl)) +
    geom_hline(yintercept = power, linetype = 2) +
    geom_rect(
      xmin = min(new_data$trend[new_data$ucl >= power]),
      xmax = min(new_data$trend[new_data$lcl >= power]),
      ymin = -Inf, ymax = Inf, alpha = 0.1
    ) +
    geom_vline(
      xintercept = min(new_data$trend[new_data$fit >= power]),
      linetype = 2
    ) +
    geom_errorbar() +
    geom_point(aes(size = sims, y = simpower), alpha = 0.2) +
    geom_ribbon(data = new_data, alpha = 0.1, fill = "blue") +
    geom_line(data = new_data, aes(y = fit), colour = "blue") +
    scale_y_continuous("power", limits = c(0, 1)) +
    scale_size_continuous(limits = c(0, NA))
  print(p)
  return(
    c(
      estimate = min(new_data$trend[new_data$fit >= power]),
      lcl = min(new_data$trend[new_data$ucl >= power]),
      ucl = min(new_data$trend[new_data$lcl >= power])
    )
  )
}

estimate_average <- function(
  n_period = 2, duration = 6, n_sample = 490, initial = 87.8,
  log_sd_location = 0.21, log_sd_error = 0.21, alpha = 0.1, power = 0.9,
  connection = duckdb::dbConnect(
    duckdb::duckdb(), dbdir = "data/power_cmon.duckdb", read_only = FALSE
  )
) {
  stopifnot(require(assertthat), require(duckdb), require(tidyverse))
  assert_that(
    is.count(n_period), n_period >= 2, is.count(duration), is.count(n_sample),
    is.number(initial),
    is.number(log_sd_location), noNA(log_sd_location), log_sd_location > 0,
    is.number(log_sd_error), noNA(log_sd_error), log_sd_error > 0,
    is.number(alpha), noNA(alpha), alpha > 0, alpha < 1,
    is.number(power), noNA(power), power > 0, power < 1
  )
  if (!"average" %in% dbListTables(conn = connection)) {
    #   dbRemoveTable(conn = connection, name = "average")
    data.frame(
      trend = numeric(0), duration = integer(0), n_sample = integer(0),
      initial = numeric(0), log_sd_location = numeric(0), n_period = integer(0),
      log_sd_error = numeric(0), p = numeric(0)
    ) |>
      dbCreateTable(conn = connection, name = "average")
  }
  query <- sprintf(
    "SELECT trend, SUM(p < %f) AS significant, COUNT(trend) AS sims
FROM average
WHERE
  ABS(log_sd_location - %f) < 1e-4 AND ABS(log_sd_error - %f) < 1e-4 AND
    n_sample = %i AND duration = %i AND ABS(initial - %f) < 1e-4 AND
    n_period = %i
GROUP BY trend
ORDER BY trend",
    alpha, log_sd_location, log_sd_error, n_sample, duration, initial, n_period
  )
  sims <- dbGetQuery(conn = connection, statement = query)
  if (nrow(sims) == 0) {
    message("initializing to 0.004")
    data.frame(
      trend = 0.004, duration = duration, n_sample = n_sample,
      initial = initial, log_sd_location = log_sd_location,
      log_sd_error = log_sd_error, n_period = n_period,
      p = sim_change(
        trend = 0.004, duration = duration, n_sample = n_sample,
        initial = initial, log_sd_location = log_sd_location,
        log_sd_error = log_sd_error, n_period = n_period
      )
    ) |>
      dbAppendTable(conn = connection, name = "average")
  }
  dbGetQuery(conn = connection, statement = query) |>
    mutate(simpower = .data$significant / .data$sims) -> sims
  while (max(sims$simpower) < power) {
    extra <- 2 * max(sims$trend)
    message("increasing to ", extra)
    data.frame(
      trend = extra, duration = duration, n_sample = n_sample,
      initial = initial, log_sd_location = log_sd_location,
      log_sd_error = log_sd_error, n_period = n_period,
      p = sim_change(
        trend = extra, duration = duration, n_sample = n_sample,
        initial = initial, log_sd_location = log_sd_location,
        log_sd_error = log_sd_error, n_period = n_period
      )
    ) |>
      dbAppendTable(conn = connection, name = "average")
    dbGetQuery(conn = connection, statement = query) |>
      mutate(simpower = .data$significant / .data$sims) -> sims
  }
  while (min(sims$simpower) > power) {
    extra <- round(min(sims$trend) / 2, 4)
    stopifnot(extra > 0)
    message("decreasing to ", extra)
    data.frame(
      trend = extra, duration = duration, n_sample = n_sample,
      initial = initial, log_sd_location = log_sd_location,
      log_sd_error = log_sd_error, n_period = n_period,
      p = sim_change(
        trend = extra, duration = duration, n_sample = n_sample,
        initial = initial, log_sd_location = log_sd_location,
        log_sd_error = log_sd_error, n_period = n_period
      )
    ) |>
      dbAppendTable(conn = connection, name = "average")
    dbGetQuery(conn = connection, statement = query) |>
      mutate(simpower = .data$significant / .data$sims) -> sims
  }
  stopifnot(nrow(sims) > 1)
  sims |>
    mutate(
      not = .data$sims - .data$significant
    ) |>
    glm(formula = cbind(significant, not) ~ trend, family = binomial) -> model
  new_data <- data.frame(trend = seq(0.0001, 2 * max(sims$trend), by = 0.0001))
  new_data |>
    predict(object = model, se.fit = TRUE) -> preds
  new_data |>
    mutate(
      fit = plogis(preds$fit),
      lcl = qnorm(0.025, preds$fit, preds$se.fit) |>
        plogis(),
      ucl = qnorm(0.975, preds$fit, preds$se.fit) |>
        plogis()
    ) -> candidate
  candidate |>
    filter(.data$lcl < power, power < .data$ucl) |>
    bind_rows(
      candidate |>
        filter(.data$lcl > power) |>
        slice_min(.data$trend, n = 1)
    ) |>
    left_join(sims, by = "trend") |>
    mutate(to_do = 1000 - replace_na(.data$sims, 0)) |>
    filter(.data$to_do > 0) -> candidate
  while (nrow(candidate) > 0) {
    if (nrow(candidate) == 1) {
      extra <- candidate$trend
    } else {
      extra <- sample(candidate$trend, size = 1, prob = candidate$to_do)
    }
    message("interpolate ", extra)
    data.frame(
      trend = extra, duration = duration, n_sample = n_sample,
      initial = initial, log_sd_location = log_sd_location,
      log_sd_error = log_sd_error, n_period = n_period,
      p = sim_change(
        trend = extra, duration = duration, n_sample = n_sample,
        initial = initial, log_sd_location = log_sd_location,
        log_sd_error = log_sd_error, n_period = n_period
      )
    ) |>
      dbAppendTable(conn = connection, name = "average")
    dbGetQuery(conn = connection, statement = query) |>
      mutate(simpower = .data$significant / .data$sims) -> sims
    sims |>
      mutate(
        not = .data$sims - .data$significant
      ) |>
      glm(formula = cbind(significant, not) ~ trend, family = binomial) -> model
    new_data <- data.frame(
      trend = seq(0.0001, 2 * max(sims$trend), by = 0.0001)
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
      ) -> candidate
    p <- sims |>
      mutate(
        lcl = map2(.data$significant, .data$sims, ~ binom.test(.x, .y)) |>
          map(~ .x$conf.int),
        ucl = map_dbl(.data$lcl, ~ .x[2]),
        lcl = map_dbl(.data$lcl, ~ .x[1])
      ) |>
      ggplot(aes(x = trend, ymin = lcl, ymax = ucl)) +
      geom_hline(yintercept = power, linetype = 2) +
      geom_rect(
        xmin = min(candidate$trend[candidate$ucl >= power]),
        xmax = min(candidate$trend[candidate$lcl >= power]),
        ymin = -Inf, ymax = Inf, alpha = 0.1
      ) +
      geom_vline(
        xintercept = min(candidate$trend[candidate$fit >= power]),
        linetype = 2
      ) +
      geom_errorbar() +
      geom_point(aes(size = sims, y = simpower)) +
      geom_ribbon(data = candidate, alpha = 0.1, fill = "blue") +
      geom_line(data = candidate, aes(y = fit), colour = "blue") +
      scale_y_continuous("power", limits = c(0, 1)) +
      scale_size_continuous(limits = c(0, NA))
    print(p)
    candidate |>
      filter(.data$lcl < power, power < .data$ucl) |>
      bind_rows(
        candidate |>
          filter(.data$lcl > power) |>
          slice_min(.data$trend, n = 1)
      ) |>
      left_join(sims, by = "trend") |>
      mutate(to_do = 1000 - replace_na(.data$sims, 0)) |>
      filter(.data$to_do > 0) -> candidate
  }
  sims |>
    mutate(
      not = .data$sims - .data$significant
    ) |>
    glm(formula = cbind(significant, not) ~ trend, family = binomial) -> model
  new_data <- data.frame(trend = seq(0.0001, max(sims$trend), by = 0.0001))
  new_data |>
    predict(object = model, se.fit = TRUE) -> preds
  new_data |>
    mutate(
      fit = plogis(preds$fit),
      lcl = qnorm(0.025, preds$fit, preds$se.fit) |>
        plogis(),
      ucl = qnorm(0.975, preds$fit, preds$se.fit) |>
        plogis()
    ) -> new_data
  p <- sims |>
    mutate(
      lcl = map2(.data$significant, .data$sims, ~ binom.test(.x, .y)) |>
        map(~ .x$conf.int),
      ucl = map_dbl(.data$lcl, ~ .x[2]),
      lcl = map_dbl(.data$lcl, ~ .x[1])
    ) |>
    ggplot(aes(x = trend, ymin = lcl, ymax = ucl)) +
    geom_hline(yintercept = power, linetype = 2) +
    geom_rect(
      xmin = min(new_data$trend[new_data$ucl >= power]),
      xmax = min(new_data$trend[new_data$lcl >= power]),
      ymin = -Inf, ymax = Inf, alpha = 0.1
    ) +
    geom_vline(
      xintercept = min(new_data$trend[new_data$fit >= power]),
      linetype = 2
    ) +
    geom_errorbar() +
    geom_point(aes(size = sims, y = simpower)) +
    geom_ribbon(data = new_data, alpha = 0.1, fill = "blue") +
    geom_line(data = new_data, aes(y = fit), colour = "blue") +
    scale_y_continuous("power", limits = c(0, 1)) +
    scale_size_continuous(limits = c(0, NA))
  print(p)
  return(
    c(
      estimate = min(new_data$trend[new_data$fit >= power]),
      lcl = min(new_data$trend[new_data$ucl >= power]),
      ucl = min(new_data$trend[new_data$lcl >= power])
    )
  )
}

# var bos: 0.17
# log_sd_location: sqrt(0.17 - 0.085 / 2) = 0.36
# log_sd_error = sqrt(0.085 / 2) = 0.21
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
# log_sd_location: sqrt(0.38 - 0.085 / 2) = 0.58
# log_sd_error = sqrt(0.085 / 2) = 0.21
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

# var: akker: 0.085
# log_sd_location = sqrt(0.085 / 2) = 0.21
# log_sd_error = sqrt(0.085 / 2) = 0.21
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

# var: grasland: 0.085
# log_sd_location = sqrt(0.085 / 2) = 0.21
# log_sd_error = sqrt(0.085 / 2) = 0.21
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
# log_sd_location: sqrt(0.14 - 0.085 / 2) = 0.32
# log_sd_error = sqrt(0.085 / 2) = 0.21
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
