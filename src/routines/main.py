import argparse
from numpy import arange
import pathlib
import os
import time
import multiprocessing as mp
import sys
import copy

class KeyboardInterruptError(Exception):
    pass


class Parameters:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def collectParameters(cmdinput) -> "list[Parameters]":
    parameterList = []
    start, step, stop = cmdinput["intervals"]

    seeds = list(range(1, cmdinput["environments"] + 1))

    # .step through range of environments, paramaterized by seed value
    for seed in seeds:
        # .step through range of parameter values
        for parameterVal in arange(start, stop + step, step=step):
            cdict = copy.deepcopy(cmdinput)

            # .modify scanned parameter and savepath
            cdict[cmdinput["parameter"]] = parameterVal
            parameters = Parameters(**cdict)
            parameters.savepath = f'{cdict["savepath"]}'
            parameters.rngSeed = seed

            parameterList.append((parameters))
    return parameterList


def to_list(arg):
    return [float(i) for i in arg.split(",")]


def runJulia(prms: Parameters):
    cmd_overwrite = "--rewrite" if prms.overwrite else ""
    cmd_standingVar = "--standing_variation" if prms.standing_variation else ""
    cmd_detailed_anal = "--detailed_analytics" if prms.detailed_analytics else ""

    command = f"""julia --project=@. -O3 {prms.model}/main.jl \
            --numberTrials {prms.numberTrials} \
            --numberSamples {int(prms.numberSamples)} \
            --radius {int(prms.radius)} \
            --density {float(prms.density)} \
            --width {int(prms.dims[0])} \
            --height {int(prms.dims[1])} \
            --rngSeed {int(prms.rngSeed)} \
            --outputPath {prms.savepath} \
            --selection {round(prms.selection, 4)} \
            --intensity {round(prms.intensity, 4)} \
            --mutation {round(prms.mutation, 4)} \
            --env_type {prms.env_type} \
            --separation {int(prms.separation)} \
            --initial_type {prms.initial_type} \
            {cmd_overwrite} \
            {cmd_standingVar} \
            {cmd_detailed_anal} \
            1> /dev/null
        """
    os.system(command)


def main(nworkers, params, simInfoFile):
    pool = mp.Pool(processes=nworkers)
    t1 = time.time()
    try:
        print("starting the pool map")
        pool.map(runJulia, params)
        pool.close()
        print("pool map complete")
    except KeyboardInterrupt:
        print("got ^C while pool mapping, terminating the pool")
        pool.terminate()
        print("pool is terminated")
    except Exception as e:
        print("got exception: %r, terminating the pool" % (e,))
        pool.terminate()
        print("pool is terminated")
    finally:
        print("joining pool processes")
        pool.join()
        print("join complete")

    msg = f"\nfinished. Elapsed wall time: {round((time.time()-t1)/3600, 3)} hour(s), at {time.ctime()}"
    print(msg)

    with open(simInfoFile, mode="a") as historyFile:
        historyFile.write(msg)


# doc: save run information to run log file, append if necessary
def saveLogs(cmdinput) -> str:
    options = f"start time: {time.ctime()}\n"
    options += "  + command: python " + f"{' '.join(sys.argv)}\n"
    options += "\n".join("  + {}: {}".format(k, v) for k, v in cmdinput.items())
    options += "\n\n"
    print(options)

    savepath = cmdinput["savepath"]
    pathlib.Path(f"{savepath}/logs").mkdir(parents=True, exist_ok=True)
    simInfoFile = f"{savepath}/logs/summary_simulation_info.log"
    with open(simInfoFile, mode="a") as historyFile:
        historyFile.write(options)
    return simInfoFile


# doc: parse cmd input generate
def parseCmds(args):
    cmdinput = vars(args)

    cmdinput[args.parameter] = "xxx"
    pc = Parameters(**cmdinput)

    names = [
        f"XY:{int(pc.dims[0])},{int(pc.dims[1])}",
        f"nENV:{int(pc.environments)}",
        f"s:{round(pc.selection, 3) if type(pc.selection) != str else pc.selection}",
        f"m:{round(pc.mutation, 3) if type(pc.mutation) != str else pc.mutation}",
        f"v:{round(pc.intensity, 2) if type(pc.intensity) != str else pc.intensity}",
        f"d:{round(pc.density, 3) if type(pc.density) != str else pc.density}",
        f"r:{int(pc.radius) if type(pc.radius) != str else pc.radius}",
        f'intervals:{",".join([str(round(s, 3)) for s in pc.intervals])}',
        f"trials:{int(pc.numberTrials)}",
    ]

    if args.env_type != "uniform":
        names.append(f"EV:{args.env_type}")
    if args.standing_variation:
        names.append(args.initial_type)
        names.append("sv")
    if args.detailed_analytics:
        names.append("da")

    # .update parameter dictionary with new savepath
    cmdinput["savepath"] = f"{args.savepath}/" + "_".join(names)
    return cmdinput


parser = argparse.ArgumentParser()
parser.add_argument("--environments", type=int, default=1)
parser.add_argument("--numberTrials", type=int, default=1)
parser.add_argument("--numberSamples", type=int, default=50)

parser.add_argument("--dims", required=True, type=to_list)

parser.add_argument("--rngSeed", type=int, default=1)
parser.add_argument("--overwrite", action="store_true")
parser.add_argument("--savepath", required=True)

parser.add_argument("--separation", type=int, default=100)
parser.add_argument("--env_type", type=str, default="uniform")
parser.add_argument("--initial_type", type=str, default="alt")
parser.add_argument("--standing_variation", action="store_true")
parser.add_argument("--detailed_analytics", action="store_true")

parser.add_argument("--mutation", required=True, type=float)
parser.add_argument("--selection", required=True, type=float)
parser.add_argument("--intensity", required=True, type=float)
parser.add_argument("--radius", type=int, default=10)
parser.add_argument("--density", type=float, default=0.09)

parser.add_argument("--parameter", default="mutation")
parser.add_argument("--intervals", required=True, type=to_list)

parser.add_argument("--num_threads", type=int, default=10)
parser.add_argument("--background", action="store_true")

parser.add_argument("--model", type=str, default="src/base/")

if __name__ == "__main__":
    args = parser.parse_args()
    cmdinput = parseCmds(args)
    simInfoFile = saveLogs(cmdinput)

    if not args.background:
        print(f"do you wish to proceed?")
        proceed = input("[y/N]: ")

        if proceed.lower() == "n":
            print("you decided to not proceed. exiting.")
            exit()

    print("submitting your jobs.")

    mp.freeze_support()
    main(args.num_threads, collectParameters(cmdinput), simInfoFile)
