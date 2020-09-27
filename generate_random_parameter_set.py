import argparse
import random as rd

import numpy as np

def run():

    parser = argparse.ArgumentParser()
    parser.add_argument("--save-location", type=str, default="")
    parser.add_argument("--kill-radius", type=int, nargs=2, default=[0,10])
    parser.add_argument("--inhibit-radius", type=int, nargs=2, default=[0,10])
    parser.add_argument("--grow-radius", type=int, nargs=2, default=[0,10])
    parser.add_argument("--mutation", type=float, nargs=2, default=[0.0, 0.1])
    parser.add_argument("--kill-margin", type=float, nargs=2, default=[0.0, 0.5])
    parser.add_argument("--inhib-margin", type=float, nargs=2, default=[0.0, 0.5])
    parser.add_argument("--num-jobs", type=int, default=10)
    parser.add_argument("--seed", type=int, default=12)
    args = parser.parse_args()

    # meta-job.
    file_location=      args.save_location
    kill_radius=        args.kill_radius
    inhib_radius=       args.inhibit_radius
    grow_radius=        args.grow_radius
    mutation=           args.mutation
    kill_margin=        args.kill_margin
    inhib_margin=       args.inhib_margin

    num_jobs=           args.num_jobs

    rd.seed(args.seed)

    params_list = []
    for i in range(num_jobs):
        k = rd.randint(kill_radius[0], kill_radius[1])
        i = rd.randint(inhib_radius[0], inhib_radius[1])
        g = rd.randint(grow_radius[0], grow_radius[1])
        m = rd.uniform(mutation[0], mutation[1])
        km = rd.uniform(kill_margin[0], kill_margin[1])
        im = rd.uniform(inhib_margin[0], inhib_margin[1])
        params_list.append(
                [
                    k, i, g, m, km, im
                    ]
                )

    header = "K-RADIUS, I-RADIUS, G-RADIUS, MUTATE, K-MARGIN, I-MARGIN"
    np.savetxt(file_location, params_list, delimiter=',', header=header)



if __name__ == "__main__":
    run()
