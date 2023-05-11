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
  return(sum(points$cost))
}

buffer_treatments <- function(points, treatments_raster, treatments_file, cost_raster) {
  rasterization_factor <- 10
  buffers <- terra::buffer(points, width = points$buffer_size)
  high_res_raster <- terra::disagg(treatments_raster, fact = rasterization_factor)
  buffers_raster <- terra::rasterize(buffers, high_res_raster)
  count_raster <- terra::aggregate(buffers_raster,
    fact = rasterization_factor,
    fun = function(x) {
      length(x[!is.na(x)])
    }
  )
  treatments_raster <- count_raster / terra::global(count_raster, "max")[[1]]
  treatment_cost_raster <- treatments_raster * cost_raster
  terra::writeRaster(treatments_raster, treatments_file, overwrite = TRUE)
  actual_cost <- terra::global(treatment_cost_raster, "sum", na.rm = T)[[1]]
  return(actual_cost)
}

treatments <- function(points, treatments_raster, treatments_file, cost_raster) {
  if ("buffer_size" %in% names(points)) {
    return(buffer_treatments(
      points,
      treatments_raster,
      treatments_file,
      cost_raster
    ))
  }
  return(pixel_treatments(points, treatments_raster, treatments_file))
}

run_rpops <- function(parameters, management = TRUE) {
  if (!management) {
    parameters$management <- FALSE
  }
  output <- do.call(PoPS::pops_multirun, parameters)
  # evaluate by area
  area <- output$infected_areas
  area <- area[[1, ncol(area)]]
  return(area)
}

grass_lazy_init <- function(infected_file, host_file, total_file, weather_file) {
  library("rgrass")
  raster <- terra::rast(infected_file)
  grassbin <- system("grass8 --config path", intern = TRUE)
  # start GRASS GIS from R
  initGRASS(
    gisBase = grassbin,
    home = tempdir(),
    SG = raster,
    gisDbase = tempdir(),
    location = "optimization",
    mapset = "PERMANENT",
    override = TRUE,
    remove_GISRC = TRUE
  )
  host <- "host"
  total <- "total"
  infected <- "infected"
  weather <- "weather"
  average <- "average_R_temp"
  weather_list_file <- file.path(tempdir(), "weather.txt")
  if (!file.exists(weather_list_file)) {
    execGRASS("r.in.gdal",
      input = infected_file,
      output = infected,
      flags = c("quiet", "o")
    )
    execGRASS("r.in.gdal",
      input = host_file,
      output = host,
      flags = c("quiet", "o")
    )
    execGRASS("r.in.gdal",
      input = total_file,
      output = total,
      flags = c("quiet", "o")
    )
    execGRASS("r.in.gdal",
      input = weather_file,
      output = weather,
      flags = c("quiet", "o")
    )
    execGRASS("g.list",
      pattern = "weather*",
      flags = c("quiet", "m"),
      type = "raster",
      output = weather_list_file
    )
  }
  return(list(
    host = host,
    total = total,
    infected = infected,
    average = average,
    weather_file = weather_list_file
  ))
}

run_r.pops.spread <- function(parameters, management = TRUE) {
  grassparams <- grass_lazy_init(
    parameters$infected_file,
    parameters$host_file,
    parameters$total_populations_file,
    parameters$temperature_coefficient_file
  )
  if (management) {
    treatments <- "treatments_R_temp"
    execGRASS("r.in.gdal",
      input = parameters$treatments_file,
      output = treatments,
      flags = c("o", "overwrite", "quiet")
    )
  } else {
    treatments <- ""
  }
  execGRASS(
    "r.pops.spread",
    model_type = parameters$model_type,
    host = grassparams$host,
    total_plants = grassparams$total,
    infected = grassparams$infected,
    start_date = parameters$start_date,
    end_date = parameters$end_date,
    step_unit = parameters$time_step,
    step_num_units = 1,
    treatments = treatments,
    treatment_date = parameters$treatment_dates,
    treatment_length = 0,
    treatment_application = "all_infected_in_cell",
    reproductive_rate = parameters$parameter_means[[1]],
    natural_distance = parameters$parameter_means[[2]],
    anthropogenic_distance = parameters$parameter_means[[4]],
    natural_direction_strength = parameters$parameter_means[[5]],
    anthropogenic_direction_strength = parameters$parameter_means[[6]],
    natural_direction = "N",
    anthropogenic_direction = "N",
    natural_dispersal_kernel = "exponential",
    anthropogenic_dispersal_kernel = "cauchy",
    percent_natural_dispersal = parameters$parameter_means[[3]],
    seasonality = "1,12",
    weather_coefficient_file = grassparams$weather_file,
    average = grassparams$average,
    output_frequency = "every_n_steps",
    output_frequency_n = 400,
    random_seed = parameters$random_seed,
    runs = parameters$number_of_iterations,
    nprocs = parameters$number_of_cores,
    flags = c("overwrite", "quiet")
  )
  infected_area <- execGRASS("r.info", map = grassparams$average, flags = c("g", "e"), intern = TRUE)
  infected_area <- tail(strsplit(infected_area[[26]], " ")[[1]], 1)
  infected_area <- as.numeric(gsub("\"", "", infected_area))

  return(infected_area)
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
    cost <- treatments(
      candidate,
      treatments_raster,
      pops_parameters$treatments_file,
      cost_raster
    )
    evaluated <- run_pops(pops_parameters)
    if (evaluated / baseline < threshold) {
      for (cat in as.character(candidate$cat)) {
        new_weights[[cat]] <- new_weights[[cat]] + 1
      }
      evaluated_list[evaluated_index] <- evaluated
      evaluated_index <- evaluated_index + 1
      actual_cost[evaluated_index] <- cost
      if (evaluated < minim_evaluated) {
        minim_evaluated <- evaluated
        best_candidate <- candidate
      }
    }
    acceptance_rate <- length(evaluated_list) / tested
    particles_left <- min_particles - length(evaluated_list)
    min_evaluated <- minim_evaluated / baseline
    # print intermediate
    if (tested %% 10 == 0) {
      message(
        "Rate: ", acceptance_rate,
        ", particles left: ", particles_left,
        ", min evaluated: ", min_evaluated
      )
    }
  }
  perc <- quantile(evaluated_list, probs = threshold_percentile / 100)
  output <- list(
    weights = new_weights,
    acceptance_rate = acceptance_rate,
    threshold = perc / baseline,
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
                     GRASS = TRUE,
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
  multi_pops_parameters$GRASS <- NULL

  # grass_init()

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

  run_pops <<- run_rpops
  if (GRASS) {
    run_pops <<- run_r.pops.spread
  }

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
    filtered_points[[new_weight_column]] <- results$weights[match(
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
    after <- tmp[tmp[[new_weight_column]] > 0, ]
    if (sum(after$cost) >= budget) {
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

  output <- list(
    best_candidate = results$best_candidate,
    best_evaluation = results$best_candidate_evaluation,
    weights = filtered_points,
    acceptance_rates = acceptance_rates,
    thresholds = thresholds,
    actual_cost_mean = actual_cost_mean
  )
  return(output)
}

library("PoPS")

potential_file <- "/home/akratoc/dev/pops/optimization/potential.tif"
pixel_cost_file <- "/home/akratoc/dev/pops/optimization/walk_cost.tif"
infected_file <- "/home/akratoc/dev/pops/optimization/infected_2019.tif"
host_file <- "/home/akratoc/dev/pops/optimization/host.tif"
total_file <- "/home/akratoc/dev/pops/optimization/total.tif"
weather_file <- "/home/akratoc/dev/pops/optimization/weather.tif"
buffer_size_file <- "/home/akratoc/dev/pops/optimization/buffer_size.tif"
buffer_cost_file <- "/home/akratoc/dev/pops/optimization/buffer_cost.tif"
results <- optimize(
  infestation_potential_file = potential_file,
  cost_file = buffer_cost_file,
  buffer_size_file = buffer_size_file,
  min_particles = 10,
  budget = 80000,
  filter_percentile = 15,
  threshold_percentile = 10,
  infected_file = infected_file,
  host_file = host_file,
  total_populations_file = total_file,
  parameter_means = c(1.6, 25, 0.99, 4000, 1, 1, 100, 200),
  parameter_cov_matrix = matrix(data = 0, nrow = 8, ncol = 8),
  time_step = "week",
  start_date = "2019-01-01",
  end_date = "2023-12-31",
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
  temperature_coefficient_file = weather_file,
  GRASS = TRUE
)
