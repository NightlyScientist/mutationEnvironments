import os
import pandas as pd
from dataclasses import dataclass


@dataclass
class Paths:
    base: str
    img: str


@dataclass
class EnsemblePaths:
    all_paths: list[Paths]

    def __getitem__(self, index):
        return self.all_paths[index]

    def __iter__(self):
        return iter(self.all_paths)

def fetchSingleRunPath(input_path: str) -> EnsemblePaths:
    directories = []
    for sub_path in os.listdir(input_path):
        abs_path = os.path.join(input_path, sub_path)

        # check if summary file exist to verify it's a valid ensemble rn
        info_file_path = os.path.join(abs_path, "logs", "search_table_information.log")

        if os.path.isdir(abs_path) and os.path.isfile(info_file_path):
            directories.append(abs_path)

    allPaths = [setOutputPaths(x) for x in directories]
    return EnsemblePaths(allPaths)


def fetchEnsemblePaths(input_path: str) -> EnsemblePaths:
    directories = []
    for sub_path in os.listdir(input_path):
        abs_path = os.path.join(input_path, sub_path)

        # check if summary file exist to verify it's a valid ensemble rn
        info_file_path = os.path.join(abs_path, "logs", "search_table_information.log")

        if os.path.isdir(abs_path) and os.path.isfile(info_file_path):
            directories.append(abs_path)

    allPaths = [setOutputPaths(x) for x in directories]
    return EnsemblePaths(allPaths)


def setOutputPaths(input_path=None) -> Paths:
    basePath = input("base path") if input_path is None else input_path

    # verify path exists and is none-empty
    if not os.path.exists(basePath):
        raise FileNotFoundError(f"Path {basePath} does not exist")

    if not os.listdir(basePath):
        raise FileNotFoundError(f"Path {basePath} is empty")

    # create images directory if it does not exist
    imgPath = os.path.join(basePath, "images")
    if not os.path.exists(imgPath):
        os.makedirs(imgPath)
    return Paths(basePath, imgPath)


def with_iter(iterable):
    with iterable as iter:
        yield from iter


def fetchWorkspaceEnv(rel_path: str) -> dict[str, str]:
    """easily set .paths file to point to relavent path, otherwise fallback to cli input, and lastly resort to input()"""
    options = dict()

    file_name = os.path.realpath(os.path.join(rel_path, ".workspace_env"))
    if not os.path.isfile(file_name):
        print(".workspace_env does not exist or not found.")
        options["top_level_path"] = input("top_level_path: ")
        options["example_ensemble_path"] = input("example_ensemble_path: ")
        return options

    for line in with_iter(open(file_name, "r")):
        if line.startswith("#") or line.strip() == "":
            continue
        key_value = [x.strip().replace('"', "") for x in line.split("=")]
        assert len(key_value) == 2, f"parsing error in .workspace_env: {key_value}"
        options[key_value[0]] = key_value[1]
        print(f"key {key_value[0]}: value {key_value[1]}")

    assert "top_level_path" in options, "missing top_level in .workspace_env"
    return options


def ensembleTableInfo(input_path: str):
    """from an ensemble folder, specificed by input_path, read all options and generate a table summarizing the entire parameter space"""

    dictList: list[dict[str, str]] = []
    for root, dirs, _ in os.walk(input_path):
        for directory in dirs:
            sub_path = os.path.join(root, directory)
            # if "env" not in sub_path and "ENV" not in sub_path:
            #     continue
            options_file = os.path.join(sub_path, "options.csv")
            if not os.path.isfile(options_file):
                continue

            opts = pd.read_csv(options_file, sep="\t", header=0).to_dict(
                orient="index"
            )[0]

            opts["path"] = sub_path

            dictList.append(opts)
    return pd.DataFrame.from_dict(dictList)
    # return pd.DataFrame.from_dict(dictList)
