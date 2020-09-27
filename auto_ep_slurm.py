import os
import math
from copy import deepcopy

import argparse

def run():
    class job:
        def __init__(self,
                file_location,
                kill_radius, inhib_radius, grow_radius,
                mutation, kill_margin, inhib_margin,
                versions):
            self.file_location = file_location
            self.params = [
                    kill_radius, inhib_radius, grow_radius,
                    mutation, kill_margin, inhib_margin, 
                    versions #NOTE: if new parameter added, check version_number is always last (because it is the SEED).
                    ]

        def __str__(self):
            return "JOB: {}, {}".format(self.params, self.file_location)

        def get_command(self, executable_name):
            cmd = "./{} {}".format(executable_name, self.file_location)
            for s in self.params:
                for p in s:
                    cmd += " {}".format(p)
            cmd += " {}".format(self.params[-1][0]) # SEED is just the version number.
            return cmd

    parser = argparse.ArgumentParser()
    parser.add_argument("--save-location", type=str, default="")
    parser.add_argument("--kill-radius", type=int, nargs=2, default=[0,10])
    parser.add_argument("--inhibit-radius", type=int, nargs=2, default=[0,10])
    parser.add_argument("--grow-radius", type=int, nargs=2, default=[0,10])
    parser.add_argument("--mutation", type=float, nargs=3, default=[0.0, 0.1, 10])
    parser.add_argument("--kill-margin", type=float, nargs=3, default=[0.0, 0.5, 10])
    parser.add_argument("--inhib-margin", type=float, nargs=3, default=[0.0, 0.5, 10])
    parser.add_argument("--versions", type=int, nargs=2, default=[0,1])
    parser.add_argument("--num-jobs", type=int, default=2)
    args = parser.parse_args()

    # meta-job.
    file_location=      args.save_location
    kill_radius=        args.kill_radius
    inhib_radius=       args.inhibit_radius
    grow_radius=        args.grow_radius
    mutation=           args.mutation
    kill_margin=        args.kill_margin
    inhib_margin=       args.inhib_margin
    versions=           args.versions

    num_jobs=           args.num_jobs

    jobs = [
            job(file_location,
                kill_radius, inhib_radius, grow_radius,
                mutation, kill_margin, inhib_margin,
                versions)
            ]
    max_jobs = 2**len(jobs[0].params)
    if num_jobs > max_jobs:
        raise Exception("Maximum number of jobs is {}".format(max_jobs))

    print("META JOB: {}".format(jobs[0]))

    # creates N 'embarassingly parallel' jobs. Job num must be a power of two.
    if (num_jobs & (num_jobs - 1)) != 0 and num_jobs > 0:
        raise Exception("number of jobs must be a power of 2.");

    def split_int_list(data):
        size = data[1] - data[0]
        if size + 1 <= 3:
            raise Exception("Cannot split when there are only three runs")
        mid = data[0] + int(size / 2)
        return [data[0], mid], [mid + 1, data[1]]

    def split_float_list(data):
        if data[2] <= 3:
            raise Exception("Cannot split when there are only three runs")
        rate = (data[1] - data[0]) / (data[2] - 1)
        c_mid = int(data[2] / 2)
        h1 = data[0] + rate * (c_mid-1)
        return [data[0], h1, c_mid], [h1 + rate, data[1], data[2] - c_mid]

    split_point = 0
    while num_jobs > 1:

        new_jobs = []
        for j in jobs:
            if split_point == len(j.params) - 1: # e.g. we are splitting on the version number.
                js = [[j.params[split_point][0]], [j.params[split_point][1]]]
            elif len(j.params[split_point]) == 2:
                js = split_int_list(j.params[split_point])
            elif len(j.params[split_point]) == 3:
                js = split_float_list(j.params[split_point])
            for s in js:
                new_j = deepcopy(j)
                new_j.params[split_point] = s
                new_jobs.append(new_j)
        jobs = new_jobs

        split_point += 1
        num_jobs /= 2

    print("SUB-JOBS:")
    for j in jobs:
        print(j)
        print(j.get_command('a.out'))

    def run_slurm(exec_command):
        cmd = "" + \
        "#!/bin/sh\n"+ \
        "#SBATCH -n 1\n"+ \
        "#SBATCH --mem-per-cpu=128\n"+ \
        "#SBATCH --share\n\n"+ \
        exec_command + "\n"
        return cmd

    for idx,j in enumerate(jobs):
        print("JOB {}:".format(idx))
        print(run_slurm(j.get_command('a.out')))

    exec_file = 'a.out'

    # create the sbatch bash scripts.
    for i,j in enumerate(jobs):
        cmd=run_slurm(j.get_command(exec_file))
        #os.system("echo {} >> job_{}.sh".format(cmd,i))
        with open("job3_{}.sh".format(i), "w") as f:
            cmd=run_slurm(j.get_command(exec_file))
            f.write(cmd)

    #for j in jobs:
    #    cmd = "sbatch <<{}".format(run_slurm(j.get_command(exec_file)))
    #    os.system(cmd)

if __name__ == "__main__":
    run()
