import os
import sys
import numpy as np
import pandas as pd
import io
import random
import grass.script as gs


def read_points(name):
    data = gs.read_command(
        "v.out.ascii",
        flags="c",
        input=name,
        columns="*",
        format="point",
        separator="comma",
    ).strip()
    return pd.read_csv(io.StringIO(data))


def minmax_scale(df, column):
    if df[column].max() - df[column].min() > 0:
        df[f"{column}_norm"] = (df[column] - df[column].min()) / (
            df[column].max() - df[column].min()
        )
    else:
        df[f"{column}_norm"] = 0
    return f"{column}_norm"


def prepare_treatments(df, treatment_map):
    treatments_string = df.to_csv(None, header=False, index=False)
    gs.write_command(
        "r.in.xyz",
        input="-",
        stdin=treatments_string,
        output=treatment_map,
        method="n",
        separator="comma",
        type="CELL",
        overwrite=True,
        quiet=True,
    )


def prepare_treatments_buffer(
    df, treatment_map, distance, distance_column, cost_raster
):
    treatments_string = df.to_csv(None, header=True, index=False)
    header = treatments_string.splitlines()[0].split(",")

    region = gs.region()
    gs.write_command(
        "v.in.ascii",
        input="-",
        stdin=treatments_string,
        output=f"{treatment_map}_points",
        separator="comma",
        skip=1,
        overwrite=True,
        quiet=True,
    )
    if distance_column:
        table_column_index = header.index(distance_column) + 1
        column = (
            gs.read_command("db.columns", table=f"{treatment_map}_points")
            .strip()
            .splitlines()[table_column_index]
        )
        gs.run_command(
            "v.buffer",
            input=f"{treatment_map}_points",
            output=f"{treatment_map}_buffers",
            column=column,
            layer=1,
            overwrite=True,
            quiet=True,
        )
    else:
        gs.run_command(
            "v.buffer",
            input=f"{treatment_map}_points",
            output=f"{treatment_map}_buffers",
            distance=distance,
            overwrite=True,
            quiet=True,
        )
    gs.use_temp_region()
    gs.run_command("g.region", res=region["nsres"] / 10, flags="a")
    gs.run_command(
        "v.to.rast",
        input=f"{treatment_map}_buffers",
        output=f"{treatment_map}_buffers",
        use="val",
        overwrite=True,
        quiet=True,
    )
    gs.del_temp_region()
    gs.run_command(
        "r.resamp.stats",
        input=f"{treatment_map}_buffers",
        output=f"{treatment_map}_count",
        method="count",
        overwrite=True,
        quiet=True,
    )
    info = gs.raster_info(f"{treatment_map}_count")
    gs.mapcalc(
        f"{treatment_map} = {treatment_map}_count / {info['max']}",
        quiet=True,
        overwrite=True,
    )
    gs.mapcalc(
        f"{treatment_map}_cost = {treatment_map} * {cost_raster}",
        quiet=True,
        overwrite=True,
    )
    actual_cost = float(
        gs.parse_command("r.univar", map=f"{treatment_map}_cost", flags="g")["sum"]
    )
    return actual_cost


def pops(treatments, average, weather_file, nprocs):
    gs.run_command(
        "r.pops.spread",
        host="lemma_den100m",
        total_plants="lemma_max100m",
        infected="cum_inf_2019",
        start_date="2019-01-01",
        end_date="2023-12-31",
        step_unit="week",
        treatments=treatments,
        treatment_date="2019-03-01",
        treatment_length=0,
        treatment_application="all_infected_in_cell",
        reproductive_rate=3.31,
        natural_distance=79.04,
        natural_direction_strength=0.47,
        natural_direction="N",
        natural_dispersal_kernel="exponential",
        seasonality=[1, 12],
        weather_coefficient_file=weather_file,
        average=average,
        output_frequency="every_n_steps",
        output_frequency_n="400",
        random_seed=1,
        runs=nprocs,
        nprocs=nprocs,
        quiet=True,
        overwrite=True,
    )
    infected_area = float(gs.raster_info(average)["description"].strip('"').split()[-1])
    return infected_area


def select_cats_by_cost(weights, costs, budget):
    tmp_budget = 0
    selected = []
    tries = 0
    weight_dict = weights.copy()
    cats = list(weight_dict.keys())
    weights = list(weight_dict.values())
    while True:
        cat = random.choices(cats, weights)[0]
        if tmp_budget + costs[cat] <= budget:
            tmp_budget += costs[cat]
            selected.append(cat)
            del weights[cats.index(cat)]
            cats.remove(cat)
        else:
            tries += 1
            if tries > 50:
                break
    return selected


def generation(
    df,
    treatment_map,
    average_map,
    baseline_evaluated,
    threshold,
    min_particles,
    weights,
    costs,
    budget,
    percentile,
    weather_file,
    buffer_distance,
    buffer_distance_column,
    cost_raster,
    nprocs,
):
    tested = 0
    new_weights = {key: 0 for key in weights}
    evaluated_list = []
    best_candidate_df = None
    minim_evaluated = baseline_evaluated
    actual_cost = []

    while len(evaluated_list) < min_particles:
        candidate_cats = select_cats_by_cost(weights, costs, budget)
        candidate_df = df.loc[df.cat.isin(candidate_cats)]
        if buffer_distance or buffer_distance_column:
            actual_cost.append(
                prepare_treatments_buffer(
                    candidate_df,
                    treatment_map,
                    buffer_distance,
                    buffer_distance_column,
                    cost_raster,
                )
            )
        else:
            prepare_treatments(candidate_df, treatment_map)
        evaluated = pops(treatment_map, average_map, weather_file, nprocs)
        tested += 1
        if evaluated / baseline_evaluated < threshold:
            for cat in candidate_cats:
                new_weights[cat] += 1
            evaluated_list.append(evaluated)
            # save minimum solution
            if evaluated < minim_evaluated:
                minim_evaluated = evaluated
                best_candidate_df = candidate_df

        acceptance_rate = len(evaluated_list) / tested

        if tested % 10 == 0:
            minim = median = None
            if evaluated_list:
                minim = np.min(evaluated_list) / baseline_evaluated
                median = np.median(evaluated_list) / baseline_evaluated
            print(
                f"Rate: {acceptance_rate}, particles left: {min_particles - len(evaluated_list)}, min evaluated: {minim}, median evaluated: {median}"
            )
    perc = np.percentile(evaluated_list, percentile)
    return (
        new_weights,
        acceptance_rate,
        perc / baseline_evaluated,
        (minim_evaluated, best_candidate_df),
        actual_cost,
    )


def estimate_initial_threshold(
    df,
    treatment_map,
    average_map,
    baseline_evaluated,
    runs,
    costs,
    budget,
    weather_file,
    buffer_distance,
    buffer_distance_column,
    cost_raster,
    nprocs,
):
    weights = {key: 1 for key in costs}
    evaluated = []
    for run in range(runs):
        candidate_cats = select_cats_by_cost(weights, costs, budget)
        candidate_df = df.loc[df.cat.isin(candidate_cats)]
        if buffer_distance or buffer_distance_column:
            prepare_treatments_buffer(
                candidate_df,
                treatment_map,
                buffer_distance,
                buffer_distance_column,
                cost_raster,
            )
        else:
            prepare_treatments(candidate_df, treatment_map)
        evaluated.append(pops(treatment_map, average_map, weather_file, nprocs))
    return np.percentile(evaluated, 10) / baseline_evaluated


def filter_particles(df, weights, iteration, percentile):
    new_df = df.copy()
    new_df["weights"] = new_df["cat"].map(weights)
    # filter unsuccessfull ones
    new_df = new_df.loc[new_df.weights > 1]
    weights_percentile = np.percentile(new_df.weights, percentile)
    # filter unlikely ones
    new_df = new_df.loc[new_df.weights >= weights_percentile]
    return new_df


def best_guess(
    df,
    weight_column,
    cost_column,
    budget,
    treatment_map,
    average_map,
    buffer_distance,
    buffer_distance_column,
    cost_raster,
    weather_file,
    nprocs,
):
    sorted_df = df.sort_values(weight_column, ascending=False)
    candidate_df = sorted_df[sorted_df[cost_column].cumsum() <= budget]
    if buffer_distance or buffer_distance_column:
        prepare_treatments_buffer(
            candidate_df,
            treatment_map,
            buffer_distance,
            buffer_distance_column,
            cost_raster,
        )
    else:
        prepare_treatments(candidate_df, treatment_map)
    candidate_df.to_csv(os.path.join(output, "best_guess.csv"))
    return pops(treatment_map, average_map, weather_file, nprocs)


def main(
    infected,
    potential_column,
    cost_column,
    cost_raster,
    buffer_distance,
    buffer_distance_column,
    budget,
    min_particles,
    filter_percentile,
    threshold_percentile,
    output,
    nprocs,
):
    treatment_map = gs.append_node_pid("treatments")
    average_map = gs.append_node_pid("average")
    weather_file = gs.tempfile(create=False)
    gs.run_command(
        "g.list", mapset="weather", flags="m", type="raster", output=weather_file
    )

    df = read_points(infected)
    potential_norm_column = minmax_scale(df, potential_column)
    cost_norm_column = minmax_scale(df, cost_column)
    df["prior_weight"] = (df[potential_norm_column] + 1 - df[cost_norm_column]) / 2

    initial_evaluated = pops(None, average_map, weather_file, nprocs)
    print(f"Initial evaluated: {initial_evaluated}")

    best_guess_evaluated = best_guess(
        df,
        "prior_weight",
        cost_column,
        budget,
        treatment_map,
        average_map,
        buffer_distance,
        buffer_distance_column,
        cost_raster,
        weather_file,
        nprocs,
    )
    print(f"Best guess: {best_guess_evaluated}")
    weights = dict(zip(df.cat, df.prior_weight))
    costs = dict(zip(df.cat, df[cost_column]))

    thresholds = []
    thresholds.append(
        estimate_initial_threshold(
            df,
            treatment_map,
            average_map,
            initial_evaluated,
            10,
            costs,
            budget,
            weather_file,
            buffer_distance,
            buffer_distance_column,
            cost_raster,
            nprocs,
        )
    )
    acceptance_rates = []
    actual_cost_mean = []
    filtered_dfs = []
    tmpdf = df.copy()
    iteration = 0
    while True:
        print(f"Iteration {iteration}: threshold {thresholds[-1]}")
        weights, acceptance_rate, perc_evaluated, best, actual_cost = generation(
            df,
            treatment_map,
            average_map,
            initial_evaluated,
            thresholds[-1],
            min_particles,
            weights,
            costs,
            budget,
            threshold_percentile,
            weather_file,
            buffer_distance,
            buffer_distance_column,
            cost_raster,
            nprocs,
        )
        acceptance_rates.append(acceptance_rate)
        thresholds.append(perc_evaluated)
        actual_cost_mean.append(
            sum(actual_cost) / len(actual_cost) if len(actual_cost) else None
        )
        filtered_df = filter_particles(tmpdf, weights, iteration, filter_percentile)
        if filtered_df[cost_column].sum() >= budget:
            print(f"Filtered {len(tmpdf) - len(filtered_df)} pixels from {len(tmpdf)}")
            tmpdf = filtered_df
            weights = dict(zip(tmpdf.cat, tmpdf.weights))
            filtered_dfs.append(filtered_df)
        else:
            break
        iteration += 1

    print(f"Best evaluated: {best[0]}")
    best[-1].to_csv(os.path.join(output, "best.csv"))
    tmpdf.to_csv(os.path.join(output, "last_filtered.csv"))
    print(f"Acceptance rate: {acceptance_rates}")
    print(f"Thresholds: {thresholds}")
    print(f"Actual cost: {actual_cost_mean}")


if __name__ == "__main__":
    budget = float(sys.argv[1])
    cost = sys.argv[2]
    min_particles = int(sys.argv[3])
    filter_percentile = float(sys.argv[4])
    threshold_percentile = float(sys.argv[5])
    try:
        buffer_distance = float(sys.argv[6])
        buffer_distance = buffer_distance if buffer_distance > 0 else None
        buffer_distance_column = None
    except ValueError:
        buffer_distance = None
        buffer_distance_column = sys.argv[6]
    output = sys.argv[7]
    nprocs = int(sys.argv[8])
    infected = "infected_2019"
    potential = "potential"
    cost_raster = "walk_cost"
    main(
        infected,
        potential,
        cost,
        cost_raster,
        buffer_distance,
        buffer_distance_column,
        budget,
        min_particles,
        filter_percentile,
        threshold_percentile,
        output,
        nprocs,
    )
