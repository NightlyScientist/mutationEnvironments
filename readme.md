Public Repository for the associated publication, “Geometry-Induced Competitive Release in a Meta-Population Model of Range Expansions in Disordered Environments", Jimmy Gonzalez Nuñez and Daniel A. Beller, Journal of the Royal Society Interface (2025). 

# Code Execution

A search through two parameters can be generated using the following command

>  python src/routines/explore_parameter_space.py --environments 5 --numberTrials 500 --numberSamples 50 --dims 1000.0,1000.0 --mutation 0.0 --selection 0.1 --compensation 0.0 --intensity 19.0 --radius 10 --density 0.25 --env_type uniform --initial_type uniform --standing_variation --overwrite --heatmap --parameters intensity,selection --intervals_1 0,1,7 --intervals_2 0.0,0.02,0.1 --model src/base/ --savepath workspace/revised_simulations/

A single simulation across one variable can be perfomed using 

>  python src/routines/main.py --env_type circle --separation 100 --model src/base/ --initial_type alt --environments 1 --numberTrials 500 --numberSamples 50 --dims 1000.0,1000.0 --selection 0.04 --compensation 0.0 --intensity 6.0 --mutation 0.0 --savepath workspace/revised_individual_sims/ --parameter separation --intervals=10,20,200 --num_threads 4 --density 0.25 --radius 10 --background --overwrite --heatmap

A list of simulation options and their descriptions can be found using

> python src/routines/explore_parameter_space.py --help

Source code used to generate article figures are contained in the [figures/](./figures/) directory. Analysis is performed by providing a list of directories (input variable) containing the simluation data from [explore_parameter_space.py](./src/routines/explore_parameter_space.py). If all data are saved in, for example, workarea/experiments/, the figure scripts will automatically scan and collect all simulation data and proceed with the analysis. Commands used to generate the datasets are located in the `datasets.md` file. These commmands will source the parameter space search routine to submit many jobs to the slurm queue. The output location is given by `--savepath` option, and are named based on which figure they're used to generate. 


Additionally, individual simulations can be executed using main.jl; a list of command-line options can be displayed using

> julia src/base/main.jl --help

This code requires both Julia (v1.10.4) and python (v3.11.5) to be installed. Our Julia environment is contained in the [Manifest.toml](Manifest.toml) and [Project.toml](Project.toml) files. Our Python environment is provided in [requirements.txt](requirements.txt) and in [environment.yml](environment.yml).

Julia package dependencies can be downloaded and installed for the current project using the [Project.toml](Project.toml) file through the following command:

> julia --project=. -e 'using Pkg; Pkg.instantiate()'

If using the Conda package manager for python, the python dependencies can be installed using 

> conda create --name <env_name> --file requirements.txt

While the script [explore_parameter_space.py](./src/routines/explore_parameter_space.py) to generate simulations works best with the [slurm workload manager](https://slurm.schedmd.com/overview.html) installed, the script will check for an existing slurm installion and will fallback to executing sequentially via bash if no slurm installation is found.
