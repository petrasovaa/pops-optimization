import os
import sys
import multiprocessing
import grass.script as gs


def run_buffer(params):
    x, y, column, buffer_size, buffer_size_raster, region_coef = params
    cost_name = f"cost_{x}_{y}"
    circle_name = f"circle_{x}_{y}"
    region = gs.region()
    env = os.environ.copy()
    if buffer_size_raster:
        what = gs.raster_what(buffer_size_raster, coord=(x, y))
        buffer = float(what[0][buffer_size_raster]["value"])
    elif buffer_size:
        buffer = buffer_size
    else:
        raise ValueError("No buffer specified")
    env["GRASS_REGION"] = gs.region_env(
        n=y + buffer,
        s=y - buffer,
        e=x + buffer,
        w=x - buffer,
        res=region["nsres"] / region_coef,
        flags="a",
    )
    gs.run_command(
        "r.circle",
        flags="b",
        output=circle_name,
        coordinates=(x, y),
        max=buffer,
        multiplier=region_coef * region_coef,
        env=env,
    )
    gs.mapcalc(f"{cost_name} = {column} / {circle_name}", env=env)
    cost = float(gs.parse_command("r.univar", map=cost_name, flags="g", env=env)["sum"])
    gs.run_command("g.remove", type="raster", name=[cost_name, circle_name], flags="f")
    return (x, y, cost)


if __name__ == "__main__":
    nproc = int(sys.argv[1])
    vector = "infected_2019"
    cost = "walk_cost"
    buffer_size_raster = "infestation_range"
    buffer_size_constant = None
    region_coef = 10

    data = (
        gs.read_command("v.out.ascii", input=vector, separator="comma", quiet=True)
        .strip()
        .splitlines()
    )
    params = []
    for line in data:
        x, y, cat = line.split(",")
        params.append((float(x), float(y), cost, buffer_size_constant, buffer_size_raster, region_coef))
    with multiprocessing.Pool(processes=nproc) as pool:
        results = pool.map_async(run_buffer, params).get()
        # write results as json to file
        with open("buffers.csv", "w") as f:
            for each in results:
                f.write(",".join([str(v) for v in each]))
                f.write("\n")
        gs.run_command(
            "r.in.xyz", input="buffers.csv", output="range_buffer_cost", separator="comma"
        )
