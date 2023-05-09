prior_weight <- function(cost_column, potential_column) {
  minmax_scale <- function(column) {
    if (max(column) - min(column) > 0) {
      return(scales::rescale(column, to = c(0, 1)))
    } else {
      return(0)
    }
  }
  cost_norm <- minmax_scale(cost_column)
  potential_norm <- minmax_scale(potential_column)
  return((potential_norm + 1 - cost_norm) / 2)
}

pixel_treatments <- function(points, treatments_raster, treatments_file) {
  raster <- terra::rasterize(points, treatments_raster, background = 0)
  terra::writeRaster(raster, treatments_file, overwrite = TRUE)
}

buffer_treatments <- function(points, treatments_raster, treatments_file, cost_raster) {
  rasterization_factor <- 10
  buffers <- terra::buffer(points, width = points$buffer_size[[1]])
  high_res_raster <- terra::disagg(treatments_raster, fact = rasterization_factor)
  buffers_raster <- terra::rasterize(buffers, high_res_raster)
  count_raster <- terra::aggregate(buffers_raster,
    fact = rasterization_factor,
    fun = function(x) {
      length(x[!is.na(x)])
    }
  )
  treatment_raster <- count_raster / terra::global(count_raster, "max")[[1]]
  treatment_cost_raster <- treatments_raster * cost_raster
  terra::writeRaster(treatments_raster, treatments_file, overwrite = TRUE)
  actual_cost <- terra::global(treatment_cost_raster, "sum", na.rm = T)[[1]]
  return(actual_cost)
}

treatments <- function(points, treatments_raster, treatments_file, cost_raster) {
  if ("buffer_size" %in% colnames(points)) {
    buffer_treatments(
      points,
      treatments_raster,
      treatments_file,
      cost_raster
    )
  } else {
    pixel_treatments(points, treatments_raster, treatments_file)
  }
}

run_pops <- function(parameters, management = TRUE) {
  if (!management) {
    parameters$management <- FALSE
  }
  output <- do.call(PoPS::pops_multirun, parameters)
  # evaluate by area
  area <- output$infected_areas
  area <- area[[1, ncol(area)]]
  return(area)
}

best_guess <- function(points,
                       weight_column,
                       treatments_raster,
                       cost_raster,
                       budget,
                       pops_parameters) {
  sorted <- points[order(points[[weight_column]][[1]], decreasing = TRUE), ]
  candidate <- sorted[cumsum(sorted$cost) <= budget, ]

  treatments(
    candidate,
    treatments_raster,
    pops_parameters$treatments_file,
    cost_raster
  )
  evaluation <- run_pops(pops_parameters)
  return(list(candidate = candidate, evaluation = evaluation))
}

estimate_initial_threshold <- function(points,
                                       weight_column,
                                       treatments_raster,
                                       cost_raster,
                                       budget,
                                       baseline,
                                       pops_parameters) {
  evaluated_list <- c()
  runs <- 10
  for (run in seq(runs)) {
    candidate <- sample_candidate(points, weight_column, budget)
    treatments(
      candidate,
      treatments_raster,
      pops_parameters$treatments_file,
      cost_raster
    )
    evaluated_list[run] <- run_pops(pops_parameters)
  }
  return(quantile(evaluated_list, probs = 0.1) / baseline)
}

sample_candidate <- function(points, weight_column, budget) {
  # to allow sample vector size 1
  sample.vec <- function(x, ...) x[sample(length(x), ...)]
  tmp_budget <- 0
  tries <- 0
  count <- 1
  candidate_cats <- c()
  cats <- points$cat
  weights <- points[[weight_column]][[1]]
  costs <- points$cost
  while (length(cats)) {
    sample_cat <- sample.vec(cats, 1, prob = weights)
    selected_index <- which(cats == sample_cat)
    cost <- costs[[selected_index]]
    if (tmp_budget + cost <= budget) {
      tmp_budget <- tmp_budget + cost
      candidate_cats[count] <- sample_cat
      count <- count + 1
      cats <- cats[-selected_index]
      weights <- weights[-selected_index]
      costs <- costs[-selected_index]
    } else {
      tries <- tries + 1
      if (tries > 50) {
        break
      }
    }
  }
  return(points[points$cat %in% candidate_cats, ])
}


generation <- function(points,
                       weight_column,
                       treatments_raster,
                       cost_raster,
                       budget,
                       min_particles,
                       threshold_percentile,
                       threshold,
                       baseline,
                       pops_parameters) {
  evaluated_list <- c()
  evaluated_index <- 1
  actual_cost <- c()
  tested <- 0
  minim_evaluated <- baseline
  new_weights <- setNames(as.list(rep(0, length(points$cat))), points$cat)
  while (length(evaluated_list) < min_particles) {
    tested <- tested + 1
    candidate <- sample_candidate(points, weight_column, budget)
    treatments(
      candidate,
      treatments_raster,
      pops_parameters$treatments_file,
      cost_raster
    )
    evaluated <- run_pops(pops_parameters)
    if (evaluated / baseline < threshold) {
      for (cat in as.character(candidate$cat[[1]])) {
        new_weights[[cat]] <- new_weights[[cat]] + 1
      }
      evaluated_list[evaluated_index] <- evaluated
      evaluated_index <- evaluated_index + 1
      if (evaluated < minim_evaluated) {
        minim_evaluated <- evaluated
        best_candidate <- candidate
      }
    }
    acceptance_rate <- length(evaluated_list) / tested
    # print intermediate
    if (tested %% 10 == 0) {
      if (length(evaluated_list) != 0) {
        message(
          "Rate: ", acceptance_rate,
          ", particles left: ", min_particles - length(evaluated_list),
          ", min evaluated: ", min(evaluated_list) / baseline,
          ", median evaluated: ", median(evaluated_list) / baseline
        )
      }
    }
  }
  perc <- quantile(evaluated_list, probs = threshold_percentile / 100)
  output <- list(
    weights = new_weights,
    acceptance_rate = acceptance_rate,
    threshold = perc / baseline_evaluated,
    best_candidate = best_candidate,
    best_candidate_evaluation = minim_evaluated,
    actual_cost = actual_cost
  )
  return(output)
}

filter_particles <- function(points, weight_column, iteration, percentile) {
  # filter unsuccessful ones
  filtered <- points[points[[weight_column]] > 0, ]
  weights_percentile <- quantile(filtered[[weight_column]][[1]],
    percentile / 100,
    names = FALSE
  )
  # filter unlikely ones
  points[points[[weight_column]] <= weights_percentile, weight_column] <- 0
  return(points)
}

optimize <- function(infestation_potential_file,
                     cost_file,
                     buffer_size_file = NULL,
                     min_particles,
                     budget,
                     filter_percentile = 15,
                     threshold_percentile = 10,
                     infected_file,
                     host_file,
                     total_populations_file,
                     parameter_means,
                     parameter_cov_matrix,
                     temp = FALSE,
                     temperature_coefficient_file = "",
                     precip = FALSE,
                     precipitation_coefficient_file = "",
                     model_type = "SI",
                     latency_period = 0,
                     time_step = "month",
                     season_month_start = 1,
                     season_month_end = 12,
                     start_date = "2008-01-01",
                     end_date = "2008-12-31",
                     use_survival_rates = FALSE,
                     survival_rate_month = 3,
                     survival_rate_day = 15,
                     survival_rates_file = "",
                     use_lethal_temperature = FALSE,
                     temperature_file = "",
                     lethal_temperature = -12.87,
                     lethal_temperature_month = 1,
                     mortality_on = FALSE,
                     mortality_rate = 0,
                     mortality_time_lag = 0,
                     mortality_frequency = "year",
                     mortality_frequency_n = 1,
                     management = FALSE,
                     treatment_dates = c(""),
                     treatments_file = "",
                     treatment_method = "ratio",
                     natural_kernel_type = "cauchy",
                     anthropogenic_kernel_type = "cauchy",
                     natural_dir = "NONE",
                     anthropogenic_dir = "NONE",
                     number_of_iterations = 100,
                     number_of_cores = NA,
                     pesticide_duration = 0,
                     pesticide_efficacy = 1.0,
                     random_seed = NULL,
                     output_frequency = "year",
                     output_frequency_n = 1,
                     movements_file = "",
                     use_movements = FALSE,
                     start_exposed = FALSE,
                     generate_stochasticity = TRUE,
                     establishment_stochasticity = TRUE,
                     movement_stochasticity = TRUE,
                     dispersal_stochasticity = TRUE,
                     establishment_probability = 0.5,
                     dispersal_percentage = 0.99,
                     quarantine_areas_file = "",
                     use_quarantine = FALSE,
                     use_spreadrates = FALSE,
                     use_overpopulation_movements = FALSE,
                     overpopulation_percentage = 0,
                     leaving_percentage = 0,
                     leaving_scale_coefficient = 1,
                     exposed_file = "",
                     mask = NULL,
                     write_outputs = "None",
                     output_folder_path = "",
                     network_filename = "",
                     network_movement = "walk",
                     use_initial_condition_uncertainty = FALSE,
                     use_host_uncertainty = FALSE) {
  # parameters for pops
  multi_pops_parameters <- as.list(environment())
  multi_pops_parameters$infestation_potential_file <- NULL
  multi_pops_parameters$cost_file <- NULL
  multi_pops_parameters$min_particles <- NULL
  multi_pops_parameters$budget <- NULL
  multi_pops_parameters$filter_percentile <- NULL
  multi_pops_parameters$threshold_percentile <- NULL
  multi_pops_parameters$buffer_size_file <- NULL

  # extract infected points and associated cost and potential
  infected_raster <- terra::rast(infected_file)
  infected_points <- terra::as.points(infected_raster)
  names(infected_points) <- "infected"
  infected_points <- infected_points[infected_points$infected > 0]
  # cat
  infected_points$cat <- 1:nrow(infected_points)
  # cost
  cost_raster <- terra::rast(cost_file)
  names(cost_raster) <- "cost"
  infected_points <- terra::extract(cost_raster, infected_points, bind = TRUE)
  # infestation potential
  potential_raster <- terra::rast(infestation_potential_file)
  names(potential_raster) <- "potential"
  infected_points <- terra::extract(potential_raster, infected_points, bind = TRUE)

  if (!is.null(buffer_size_file)) {
    buffer_size_raster <- terra::rast(buffer_size_file)
    names(buffer_size_raster) <- "buffer_size"
    infected_points <- terra::extract(buffer_size_raster, infected_points, bind = TRUE)
  }
  buffer_size_raster <- NULL

  # prior weights
  iteration <- 1
  weight_column <- paste0("weight_", iteration)
  infected_points[[weight_column]] <- prior_weight(infected_points$cost, infected_points$potential)

  # prepare treatment raster
  temporary_directory <- tempdir()
  treatments_raster <- terra::rast(terra::ext(infected_raster), resolution = terra::res(infected_raster))
  terra::crs(treatments_raster) <- terra::crs(infected_raster)
  multi_pops_parameters$treatments_file <- file.path(temporary_directory, "treatments.tif")
  multi_pops_parameters$management <- TRUE

  # baseline
  baseline <- run_pops(multi_pops_parameters, management = FALSE)
  message("Baseline:", baseline)

  # best guess
  best_guess <- best_guess(
    infected_points,
    weight_column,
    treatments_raster,
    cost_raster,
    budget,
    multi_pops_parameters
  )
  # best_guess$candidate
  message("Best guess:", best_guess$evaluation)

  # initial threshold
  thresholds <- c()
  thresholds[1] <- estimate_initial_threshold(
    infected_points,
    weight_column,
    treatments_raster,
    cost_raster,
    budget,
    baseline,
    multi_pops_parameters
  )

  filtered_points <- infected_points
  tmp_points <- infected_points
  acceptance_rates <- c()
  actual_cost_mean <- c()
  while (TRUE) {
    message("Iteration ", iteration, ", threshold ", thresholds[length(thresholds)])
    results <- generation(
      filtered_points,
      weight_column,
      treatments_raster,
      cost_raster,
      budget,
      min_particles,
      threshold_percentile,
      thresholds[length(thresholds)],
      baseline,
      multi_pops_parameters
    )
    acceptance_rates <- append(acceptance_rates, results$acceptance_rate)
    thresholds <- append(thresholds, results$threshold)
    if (length(results$actual_cost) != 0) {
      actual_cost_mean <- append(
        actual_cost_mean,
        mean(results$actual_cost)
      )
    }
    new_weight_column <- paste0("weight_", iteration + 1)
    new_weights <- results$weights
    filtered_df[[new_weight_column]] <- results$weights[match(
      filtered_points$cat,
      as.integer(names(results$weights))
    )]
    # filtering
    tmp <- filter_particles(
      filtered_points,
      new_weight_column,
      iteration,
      filter_percentile
    )
    before <- filtered_points[filtered_points[[new_weight_column]] > 0, ]
    after <- tmp_df[tmp_df[[new_weight_column]] > 0, ]
    if (sum(after_df$cost) >= budget) {
      message(
        "Filtered ", length(before) - length(after),
        ", remaining: ", length(after)
      )
      filtered_points <- tmp
      weight_column <- new_weight_column
    } else {
      break
    }
    iteration <- iteration + 1
  }

  results <- list(
    best_candidate = results$best_candidate,
    best_evaluation = results$best_candidate_evaluation,
    weights = filtered_points,
    acceptance_rates = acceptance_rates,
    thresholds = thresholds,
    actual_cost_mean = actual_cost_mean
  )
  return(results)
}

library("PoPS")

results <- optimize(
  infestation_potential_file = "/tmp/potential.tif",
  cost_file = "/tmp/walk_cost.tif",
  buffer_size_file = NULL,
  min_particles = 10,
  budget = 20000,
  filter_percentile = 15,
  threshold_percentile = 10,
  infected_file = "/tmp/infected_2019.tif",
  host_file = "/tmp/host.tif",
  total_populations_file = "/tmp/total.tif",
  parameter_means = c(1.6, 25, 0.99, 4000, 1, 1, 100, 200),
  parameter_cov_matrix = matrix(data = 0, nrow = 8, ncol = 8),
  time_step = "week",
  start_date = "2019-01-01",
  end_date = "2019-12-31",
  natural_kernel_type = "exponential",
  anthropogenic_kernel_type = "cauchy",
  treatment_dates = c("2019-12-01"),
  natural_dir = "N",
  anthropogenic_dir = "N",
  random_seed = 1,
  output_frequency = "year",
  output_frequency_n = 1,
  number_of_iterations = 20,
  number_of_cores = 10,
  temp = FALSE,
  temperature_coefficient_file = "/tmp/weather.tif"
)
