read_points <- function(name) {
  return(read_VECT(name, layer = "1"))
  # return (as.data.frame(data, geom="XY"))
}

get_norm_column_name <- function(column) {
  return(paste0(column, "_norm"))
}

minmax_scale <- function(df, column) {
  new_column <- get_norm_column_name(column)
  max_val <- max(df[[column]])
  min_val <- min(df[[column]])

  if (max_val - min_val > 0) {
    df[[new_column]] <- scales::rescale(df[[column]][[1]], to = c(0, 1))
  } else {
    df[[new_column]] <- 0
  }
  return(df)
}

prepare_treatments <- function(df, treatment_map, cost_raster) {
  r <- terra::rast(terra::ext(cost_raster), resolution = terra::res(cost_raster))
  terra::crs(r) <- terra::crs(cost_raster)
  write_RAST(terra::rasterize(df, r), treatment_map, flags = "overwrite")
}

prepare_treatments_buffer <- function(df, treatment_map, cost_raster, distance, distance_column) {
  rasterization_factor <- 10
  if (is.null(distance_column)) {
    buffers <- terra::buffer(df, width = distance)
  } else {
    buffers <- terra::buffer(df, width = df[[distance_column]][[1]])
  }
  r <- terra::rast(terra::ext(cost_raster), resolution = terra::res(cost_raster))
  terra::crs(r) <- terra::crs(cost_raster)
  high_res_raster <- terra::disagg(r, fact = rasterization_factor)
  buffers_raster <- terra::rasterize(buffers, high_res_raster)
  count_raster <- terra::aggregate(buffers_raster,
    fact = rasterization_factor,
    fun = function(x) {
      length(x[!is.na(x)])
    }
  )
  treatment_raster <- count_raster / terra::global(count_raster, "max")[[1]]
  treatment_cost_raster <- treatment_raster * cost_raster
  write_RAST(treatment_raster, treatment_map, flags = "overwrite")
  actual_cost <- terra::global(treatment_cost_raster, "sum", na.rm = T)[[1]]
  return(actual_cost)
}

pops <- function(treatments, average, weather_file, nprocs) {
  execGRASS(
    "r.pops.spread",
    model_type = "SI",
    host = "lemma_den100m",
    total_plants = "lemma_max100m",
    infected = "cum_inf_2019",
    start_date = "2019-01-01",
    end_date = "2023-12-31",
    step_unit = "week",
    step_num_units = 1,
    treatments = treatments,
    treatment_date = "2019-03-01",
    treatment_length = 0,
    treatment_application = "all_infected_in_cell",
    reproductive_rate = 1.6,
    natural_distance = 25,
    anthropogenic_distance = 4000,
    natural_direction_strength = 1,
    anthropogenic_direction_strength = 1,
    natural_direction = "N",
    anthropogenic_direction = "N",
    natural_dispersal_kernel = "exponential",
    anthropogenic_dispersal_kernel = "cauchy",
    percent_natural_dispersal = 0.1,
    seasonality = "1,12",
    weather_coefficient_file = weather_file,
    average = average,
    output_frequency = "every_n_steps",
    output_frequency_n = 400,
    random_seed = 1,
    runs = nprocs,
    nprocs = nprocs,
    flags = c("overwrite", "quiet")
  )
  infected_area <- execGRASS("r.info", map = average, flags = c("g", "e"), intern = TRUE)
  infected_area <- tail(strsplit(infected_area[[26]], " ")[[1]], 1)
  infected_area <- as.numeric(gsub("\"", "", infected_area))

  return(infected_area)
}

sample_candidate <- function(df, cat_column, cost_column, weight_column, budget) {
  # to allow sample vector size 1
  sample.vec <- function(x, ...) x[sample(length(x), ...)]

  tmp_budget <- 0
  tries <- 0
  count <- 1
  candidate_cats <- c()
  cats <- df[[cat_column]][[1]]
  weights <- df[[weight_column]][[1]]
  costs <- df[[cost_column]][[1]]
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
  return(df[df[[cat_column]][[1]] %in% candidate_cats, ])
}

generation <- function(df,
                       treatment_map,
                       average_map,
                       baseline_evaluated,
                       threshold,
                       min_particles,
                       cat_column,
                       cost_column,
                       weight_column,
                       budget,
                       percentile,
                       weather_file,
                       buffer_distance,
                       buffer_distance_column,
                       cost_raster,
                       nprocs) {
  evaluated_list <- c()
  actual_cost <- c()
  evaluated_index <- 1
  tested <- 0
  cats <- df[[cat_column]][[1]]
  new_weights <- setNames(as.list(rep(0, length(cats))), cats)
  minim_evaluated <- baseline_evaluated
  best_candidate_df <- NULL
  while (length(evaluated_list) < min_particles) {
    candidate_df <- sample_candidate(df, cat_column, cost_column, weight_column, budget)

    if (!is.null(buffer_distance) | !is.null(buffer_distance_column)) {
      cost <- prepare_treatments_buffer(
        candidate_df,
        treatment_map,
        cost_raster,
        buffer_distance,
        buffer_distance_column
      )
      actual_cost[evaluated_index] <- cost
      evaluated_index <- evaluated_index + 1
    } else {
      prepare_treatments(candidate_df, treatment_map, cost_raster)
    }

    evaluated <- pops(treatment_map, average_map, weather_file, nprocs)
    tested <- tested + 1
    print(evaluated / baseline_evaluated)
    if (evaluated / baseline_evaluated < threshold) {
      for (cat in as.character(candidate_df[[cat_column]][[1]])) {
        new_weights[[cat]] <- new_weights[[cat]] + 1
      }
      evaluated_list[evaluated_index] <- evaluated
      evaluated_index <- evaluated_index + 1
      if (evaluated < minim_evaluated) {
        minim_evaluated <- evaluated
        best_candidate_df <- candidate_df
      }
    }
    acceptance_rate <- length(evaluated_list) / tested
    if (tested %% 10 == 0) {
      minim <- NULL
      median <- NULL
      if (length(evaluated_list) != 0) {
        minim <- min(evaluated_list) / baseline_evaluated
        median <- median(evaluated_list) / baseline_evaluated
        print(paste(
          "Rate:", acceptance_rate, "particles left:", min_particles - length(evaluated_list),
          "min evaluated:", minim, "median evaluated:", median
        ))
      }
    }
  }
  perc <- quantile(evaluated_list, probs = percentile / 100)
  return(list(
    new_weights, acceptance_rate, perc / baseline_evaluated,
    list(minim_evaluated, best_candidate_df), actual_cost
  ))
}


estimate_initial_threshold <- function(df,
                                       treatment_map,
                                       average_map,
                                       baseline_evaluated,
                                       cat_column,
                                       cost_column,
                                       weight_column,
                                       budget,
                                       weather_file,
                                       buffer_distance,
                                       buffer_distance_column,
                                       cost_raster,
                                       runs,
                                       nprocs) {
  evaluated_list <- c()
  for (run in seq(runs)) {
    candidate_df <- sample_candidate(df, cat_column, cost_column, weight_column, budget)

    if (!is.null(buffer_distance) | !is.null(buffer_distance_column)) {
      prepare_treatments_buffer(
        candidate_df,
        treatment_map,
        cost_raster,
        buffer_distance,
        buffer_distance_column
      )
    } else {
      prepare_treatments(candidate_df, treatment_map, cost_raster)
    }
    evaluated_list[run] <- pops(treatment_map, average_map, weather_file, nprocs)
  }
  return(quantile(evaluated_list, probs = 0.1) / baseline_evaluated)
}

best_guess <- function(df,
                       treatment_map,
                       average_map,
                       cat_column,
                       cost_column,
                       weight_column,
                       budget,
                       weather_file,
                       buffer_distance,
                       buffer_distance_column,
                       cost_raster,
                       nprocs) {
  sorted_df <- df[order(df[[weight_column]][[1]], decreasing = TRUE), ]
  candidate_df <- sorted_df[cumsum(sorted_df[[cost_column]]) <= budget, ]

  if (!is.null(buffer_distance) | !is.null(buffer_distance_column)) {
    prepare_treatments_buffer(
      candidate_df,
      treatment_map,
      cost_raster,
      buffer_distance,
      buffer_distance_column
    )
  } else {
    prepare_treatments(candidate_df, treatment_map, cost_raster)
  }
  # TODO print somewhere
  return(pops(treatment_map, average_map, weather_file, nprocs))
}

filter_particles <- function(df, cat_column, weight_column, iteration, percentile) {
  # filter unsuccessful ones
  filtered_df <- df[df[[weight_column]] > 0, ]
  weights_percentile <- quantile(filtered_df[[weight_column]][[1]],
    percentile / 100,
    names = FALSE
  )
  # filter unlikely ones
  df[df[[weight_column]] <= weights_percentile, weight_column] <- 0
  return(df)
}


library("rgrass")
library("terra")

infected <- "infected_2019"
potential_column <- "potential"
cost_column <- "pixel_cost"
cost <- "walk_cost"
cat_column <- "cat"
min_particles <- 10
budget <- 20000
filter_percentile <- 15
threshold_percentile <- 10
buffer_distance <- NULL
buffer_distance_column <- NULL
nprocs <- 10
treatment_map <- "treatment_test"
average_map <- "average_test"


weather_file <- "/tmp/weather.txt"
execGRASS("g.list",
  mapset = "weather",
  flags = c("m", "overwrite"),
  type = "raster",
  output = weather_file
)

df <- read_points(infected)

df <- minmax_scale(df, potential_column)
df <- minmax_scale(df, cost_column)
iteration <- 1
weight_column <- paste0("weight_", iteration)
df[[weight_column]] <- (df[[get_norm_column_name(potential_column)]] + 1 - df[[get_norm_column_name(cost_column)]]) / 2

cost_raster <- read_RAST(cost)

baseline_evaluated <- pops("", average_map, weather_file, nprocs)
print(paste("Baseline:", baseline_evaluated))
best_guess_evaluated <- best_guess(
  df,
  treatment_map,
  average_map,
  cat_column,
  cost_column,
  weight_column,
  budget,
  weather_file,
  buffer_distance,
  buffer_distance_column,
  cost_raster,
  nprocs
)
print(paste("Best guess:", best_guess_evaluated))

thresholds <- c()
thresholds[1] <- estimate_initial_threshold(
  df,
  treatment_map,
  average_map,
  baseline_evaluated,
  cat_column,
  cost_column,
  weight_column,
  budget,
  weather_file,
  buffer_distance,
  buffer_distance_column,
  cost_raster,
  10,
  nprocs
)


acceptance_rates <- c()
actual_cost_mean <- c()
filtered_df <- df
tmp_df <- df
while (TRUE) {
  print(paste0("Iteration ", iteration, ", threshold ", thresholds[length(thresholds)]))
  results <- generation(
    filtered_df,
    treatment_map,
    average_map,
    baseline_evaluated,
    thresholds[length(thresholds)],
    min_particles, cat_column,
    cost_column,
    weight_column,
    budget,
    threshold_percentile,
    weather_file,
    buffer_distance,
    buffer_distance_column,
    cost_raster,
    nprocs
  )
  acceptance_rates <- append(acceptance_rates, results[[2]])
  thresholds <- append(thresholds, results[[3]])
  if (length(results[[5]]) != 0) {
    actual_cost_mean <- append(actual_cost_mean, sum(results[[5]]) / len(results[[5]]))
  }
  new_weight_column <- paste0("weight_", iteration + 1)
  new_weights <- results[[1]]

  filtered_df[[new_weight_column]] <- new_weights[match(
    filtered_df[[cat_column]][[1]],
    as.integer(names(new_weights))
  )]
  tmp_df <- filter_particles(
    filtered_df,
    cat_column,
    new_weight_column,
    iteration,
    filter_percentile
  )
  before_df <- filtered_df[filtered_df[[new_weight_column]] > 0, ]
  after_df <- tmp_df[tmp_df[[new_weight_column]] > 0, ]
  if (sum(after_df[[cost_column]][[1]]) >= budget) {
    print(paste0(
      "Filtered ", length(before_df) - length(after_df),
      ", remaining: ", length(after_df)
    ))
    filtered_df <- tmp_df
    weight_column <- new_weight_column
  } else {
    break
  }
  iteration <- iteration + 1
}
print(paste0("Best evaluated: ", results[[4]][[1]]))
print(paste0("Acceptance rate: ", acceptance_rates))
print(paste0("Thresholds: ", thresholds))
print(paste0("Actual cost: ", actual_cost_mean))
