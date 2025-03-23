import argparse
from numpy import arange
import os
import pathlib
import time
import sys
from pathlib import Path
from copy import deepcopy
import numpy as np


def to_list(arg):
    return [float(i) for i in arg.split(",")]


def to_string_list(arg):
    return [str(i).lower() for i in arg.split(",")]


# doc: save run information to run log file, append if necessary
def saveLogs(cmdinput) -> str:
    commit_hash = os.popen("git rev-parse HEAD").read().strip()
    options = f"git branch: {get_active_branch_name()} @{commit_hash}\n"
    options += f"start time: {time.ctime()}\n"
    options += "  + command: python " + f"{' '.join(sys.argv)}\n"
    options += "\n".join("  + {}: {}".format(k, v) for k, v in cmdinput.items())
    print(options)

    savepath = cmdinput["savepath"]
    pathlib.Path(f"{savepath}/logs").mkdir(parents=True, exist_ok=True)
    simInfoFile = f"{savepath}/logs/search_table_information.log"
    with open(simInfoFile, mode="a") as historyFile:
        historyFile.write(options)
        historyFile.write("\n\n\n")
    return simInfoFile


def get_active_branch_name():
    head_dir = Path(".") / ".git" / "HEAD"
    with head_dir.open("r") as f:
        content = f.read().splitlines()
    for line in content:
        if line[0:4] == "ref:":
            return line.partition("refs/heads/")[2]


def create_logtxt(c_opts):
    opts, to_csv = {}, {}

    _options = deepcopy(c_opts)
    args = argparse.Namespace(**_options)

    for arg in _options:
        value = getattr(args, arg)

        if isinstance(value, list) or isinstance(value, np.ndarray):
            opts[arg] = ",".join([str(v) for v in value])

            list_values = getattr(args, arg)
            for i in range(len(list_values)):
                to_csv[str(arg) + f"_{i+1}"] = str(list_values[i])
        else:
            to_csv[str(arg)] = str(getattr(args, arg))
            opts[arg] = str(value)

    opts = " ".join([f"--{k} {v}" for k, v in opts.items() if v != "" and v != " "])
    opts = opts.replace("True", "").replace("False", "")
    print(f"Options:\n\t{opts}\n")

    savePath = c_opts["savepath"]

    # create log folder and write options to file
    logPath = savePath + "/logs/"
    pathlib.Path(logPath).mkdir(parents=True, exist_ok=True)
    with open(logPath + "log.txt", mode="w") as logFile:
        print(opts, file=logFile)

    # create opt file to save input paramters
    simInfoFile = f"{savePath}/logs/opts.csv"
    with open(simInfoFile, mode="w") as optsFile:
        headers, values = zip(*to_csv.items())
        print(",".join(headers), file=optsFile)
        print(",".join(values), file=optsFile)


parser = argparse.ArgumentParser()
parser.add_argument("--environments", type=int, default=1)
parser.add_argument("--numberTrials", type=int, default=1)
parser.add_argument("--numberSamples", type=int, default=50)
parser.add_argument("--dims", required=True, type=to_list)
parser.add_argument("--savepath", required=True)

parser.add_argument("--mutation", required=True, type=float)
parser.add_argument("--selection", required=True, type=float)
parser.add_argument("--compensation", required=False, type=float, default=0.0)
parser.add_argument("--intensity", default=0.0, type=float)
parser.add_argument("--radius", type=int, default=10)
parser.add_argument("--density", type=float, default=0.0)

parser.add_argument("--env_type", type=str, default="uniform")
parser.add_argument("--initial_type", type=str, default="alt")
parser.add_argument("--standing_variation", action="store_true")
parser.add_argument("--detailed_analytics", action="store_true")
parser.add_argument("--overwrite", action="store_true")
parser.add_argument("--heatmap", action="store_true")

parser.add_argument("--parameters", default="selection,mutation", type=to_string_list)
parser.add_argument("--intervals_1", required=True, type=to_list)
parser.add_argument("--intervals_2", required=True, type=to_list)

parser.add_argument("--model", type=str, default="src/base/")
parser.add_argument("--msg", type=str, default="")

parser.add_argument("--model_flags", type=to_string_list, default=[])

if __name__ == "__main__":
    args = parser.parse_args()

    extra_flags = []
    if args.standing_variation:
        extra_flags.append("standing_variation")
    if args.detailed_analytics:
        extra_flags.append("detailed_analytics")
    if args.overwrite:
        extra_flags.append("overwrite")
    if args.heatmap:
        extra_flags.append("heatmap")

    model_flags = ",".join(args.model_flags) if len(args.model_flags) > 0 else ""

    available_parameters = [
        "selection",
        "mutation",
        "intensity",
        "radius",
        "density",
        "compensation",
    ]
    p_1, p_2 = args.parameters
    if p_1 not in available_parameters or p_2 not in available_parameters:
        print(
            "parameters given are not in the set of available options. checking names and spelling."
        )
        exit()

    start, step, stop = args.intervals_1
    mutation = args.mutation
    selection = args.selection
    compensation = args.compensation
    intensity = args.intensity
    radius = args.radius
    density = args.density

    # time will become the unique identifier for the run
    current_datetime = time.strftime("%Y-%m-%d-%H:%M")
    args_dict = vars(args)
    args_dict["savepath"] = os.path.normpath(
        args_dict["savepath"] + "/" + current_datetime
    )
    args_dict["model_flags"] = model_flags

    simFileInfo = saveLogs(args_dict)
    create_logtxt(args_dict)

    for parameterVal in arange(start, stop + step, step=step):
        if p_1 == "mutation":
            mutation = parameterVal
        elif p_1 == "intensity":
            intensity = parameterVal
        elif p_1 == "selection":
            selection = parameterVal
        elif p_1 == "compensation":
            compensation = parameterVal
        elif p_1 == "radius":
            radius = parameterVal
        elif p_1 == "density":
            density = parameterVal
        else:
            print("unable to parse your commands. Examine script/cmd")
            exit()

        cmd = f"""sbatch src/routines/slurm_main.sh \
                --model {args.model} \
                --initial_type {args.initial_type} \
                --savepath {args_dict['savepath']} \
                --environments {args.environments} \
                --numberTrials {args.numberTrials} \
                --numberSamples {args.numberSamples} \
                --dims {",".join(str(x) for x in args.dims)} \
                --selection {selection} \
                --compensation {compensation} \
                --mutation {mutation} \
                --intensity {intensity} \
                --density {density} \
                --radius {radius} \
                --parameter {p_2} \
                --env_type {args.env_type} \
                --intervals {",".join(str(x) for x in args.intervals_2)} \
                --flags {",".join(extra_flags)} \
                --model_flags {model_flags} 
                """
        os.system(cmd)
