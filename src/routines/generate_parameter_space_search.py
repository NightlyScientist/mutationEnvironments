import argparse
from numpy import arange
import os
import pathlib
import time
import sys


def to_list(arg):
    return [float(i) for i in arg.split(",")]


def to_string_list(arg):
    return [str(i).lower() for i in arg.split(",")]


# doc: save run information to run log file, append if necessary
def saveLogs(cmdinput) -> str:
    options = f"start time: {time.ctime()}\n"
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


parser = argparse.ArgumentParser()
parser.add_argument("--environments", type=int, default=1)
parser.add_argument("--numberTrials", type=int, default=1)
parser.add_argument("--numberSamples", type=int, default=50)
parser.add_argument("--dims", required=True, type=to_list)
parser.add_argument("--savepath", required=True)

parser.add_argument("--mutation", required=True, type=float)
parser.add_argument("--selection", required=True, type=float)
parser.add_argument("--intensity", required=True, type=float)
parser.add_argument("--radius", type=int, default=10)
parser.add_argument("--density", type=float, default=0.09)

parser.add_argument("--env_type", type=str, default="uniform")
parser.add_argument("--initial_type", type=str, default="alt")
parser.add_argument("--standing_variation", action="store_true")
parser.add_argument("--detailed_analytics", action="store_true")
parser.add_argument("--overwrite", action="store_true")

parser.add_argument("--parameters", default="selection,mutation", type=to_string_list)
parser.add_argument("--intervals_1", required=True, type=to_list)
parser.add_argument("--intervals_2", required=True, type=to_list)

parser.add_argument("--model", type=str, default="src/base/")


if __name__ == "__main__":
    args = parser.parse_args()

    extra_flags = []
    if args.standing_variation:
        extra_flags.append("standing_variation")
    if args.detailed_analytics:
        extra_flags.append("detailed_analytics")
    if args.overwrite:
        extra_flags.append("overwrite")

    available_parameters = ["selection", "mutation", "intensity", "radius", "density"]
    p_1, p_2 = args.parameters
    if not p_1 in available_parameters or not p_2 in available_parameters:
        print("parameters given are not in the set of available options. checking names and spelling.")
        exit()

    start, step, stop = args.intervals_1
    mutation = args.mutation
    selection = args.selection
    intensity = args.intensity
    radius = args.radius
    density = args.density

    simFileInfo = saveLogs(vars(args))

    for parameterVal in arange(start, stop + step, step=step):
        if p_1 == "mutation":
            mutation = parameterVal
        elif p_1 == "intensity":
            intensity = parameterVal
        elif p_1 == "selection":
            selection = parameterVal
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
                --savepath {args.savepath} \
                --environments {args.environments} \
                --numberTrials {args.numberTrials} \
                --numberSamples {args.numberSamples} \
                --dims {",".join(str(x) for x in args.dims)} \
                --selection {selection} \
                --mutation {mutation} \
                --intensity {intensity} \
                --density {density} \
                --radius {radius} \
                --parameter {p_2} \
                --env_type {args.env_type} \
                --intervals {",".join(str(x) for x in args.intervals_2)} \
                --flags {",".join(extra_flags)}
                """
        os.system(cmd)
