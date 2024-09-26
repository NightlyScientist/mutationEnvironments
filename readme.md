
# Code Execution

A search through two parameters can be generated using the following command

> python src/routines/generate_parameter_space_search.py --numberTrials 100 --dims 1000,1000 --savepath workspace/experiments/phase_diagram_no_mutations_phi_0.25 --mutation 0.000 --selection 0.0 --intensity 10 --radius 10 --density 0.25 --overwrite --parameters selection,intensity --intervals_1 0,0.01,0.1 --intervals_2 0,0.5,6 --model src/base/ --environments 20 --standing_variation
To execute a simulation through a single parameter search:

A single simulation across one variable can be perfomed using 

> sbatch src/scripts/slurm_main.sh --model src/base/ --initial_type alt --numberTrials 500 --dims 500,500 --mutation 0.01 --intensity 4 --density 0  --selection 0.1 --radius 5 --savepath workspace/simulation_data/stable_density --parameter selection --intervals 0,0.025,0.5 --env_type uniform --flags standing_variation

A list of simulation options and their descriptions can be found using

> python src/processing/generate_parameter_space_search.py --help

Source code used to generate article figures are contained in the [figures/](./figures/) directory. Analysis is performed by providing a list of directories (input variable) containing the simluation data from [generate_parameter_space_search.py](./src/routines/generate_parameter_space_search.py). If all data are saved in, for example, workarea/experiments/, the figure scripts will automatically scan and collect all simulation data and proceed with the analysis. Commands used to generate the datasets are located in the `datasets.md` file. These commmands will source the parameter space search routine to submit many jobs to the slurm queue. The output location is given by `--savepath` option, and are named based on which figure they're used to generate. 


Additionally, individual simulations can be executed using main.jl; a list of command-line options can be displayed using

> julia src/base/main.jl --help

This code requires both Julia (v1.10.4) and python (v3.11.5) to be installed. Our Julia environment is contained in the [Manifest.toml](Manifest.toml) and [Project.toml](Project.toml) files. Our Python environment is provided in [requirements.txt](requirements.txt) and in [environment.yml](environment.yml).

Julia package dependencies can be downloaded and installed for the current project using the [Project.toml](Project.toml) file through the following command:

> julia --project=. -e 'using Pkg; Pkg.instantiate()'

If using the Conda package manager for python, the python dependencies can be installed using 

> conda create --name <env_name> --file requirements.txt

While the script [generate_parameter_space_search.py](./src/routines/generate_parameter_space_search.py) to generate simulations works best with the [slurm workload manager](https://slurm.schedmd.com/overview.html) installed, the script will check for an existing slurm installion and will fallback to executing sequentially via bash if no slurm installation is found.
