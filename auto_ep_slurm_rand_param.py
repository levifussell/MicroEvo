import os
import math
from copy import deepcopy

import argparse

def run():
    class job:
        def __init__(self,
                file_location,
                param_file_location,
                row_lower, row_upper,
                seed):
            self.file_location = file_location
            self.param_file_location = param_file_location
            self.row_lower = row_lower
            self.row_upper = row_upper
            self.seed = seed

        def __str__(self):
            return "JOB: {}, {}, {}, {}, {}".format(self.file_location, self.param_file_location, self.row_lower, self.row_upper, self.seed)

        def get_command(self, executable_name):
            cmd = "./{} {} {} {} {} {}".format(executable_name, self.file_location, self.param_file_location, self.row_lower, self.row_upper, self.seed)
            return cmd

    parser = argparse.ArgumentParser()
    parser.add_argument("--save-location", type=str, default="")
    parser.add_argument("--param-location", type=str, default="")
    parser.add_argument("--seed", type=int, default=12)
    parser.add_argument("--num-runs", type=int, default=1000)
    parser.add_argument("--num-jobs", type=int, default=100)
    args = parser.parse_args()

    # meta-job.
    file_location=      args.save_location
    param_location=     args.param_location
    seed=               args.seed
    num_runs=           args.num_runs
    num_jobs=           args.num_jobs

    jobs = []
    runs_per_job = num_runs // num_jobs
    for i in range(num_jobs):
        jobs.append(
                job(file_location, param_location, runs_per_job * i, runs_per_job * (i+1), seed))

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
        with open("job4_{}.sh".format(i), "w") as f:
            cmd=run_slurm(j.get_command(exec_file))
            f.write(cmd)

    #for j in jobs:
    #    cmd = "sbatch <<{}".format(run_slurm(j.get_command(exec_file)))
    #    os.system(cmd)

if __name__ == "__main__":
    run()
