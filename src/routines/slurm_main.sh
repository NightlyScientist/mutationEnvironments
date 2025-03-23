#!/bin/bash
#SBATCH --job-name=mut
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --output=workspace/logs/slurm/submission_%j.txt

pwd; hostname; date

pararser() {
    # default values
    density=${density:-"0.07"}
    radius=${radius:-"10"}
    intensity=${intensity:-"0.0"}
    
    compensation=${compensation:-"0.0"}
    selection=${selection:-"0.0"}
    mutation=${mutation:-"0.0"}
    
    parameter=${parameter:-"mutation"}
    intervals=${intervals:-"0,0.1,1"}
    
    dims=${dims:-"500,1000"}
    numberTrials=${numberTrials:-"50"}
    numberSamples=${numberSamples:-"50"}
    environments=${environments:-"1"}
    
    initial_type=${initial_type:-"alt"}
    env_type=${env_type:-"uniform"}
    separation=${separation:-"100"}

    savepath=${savepath:-"workspace/simulations"}
    model=${model:-"src/base/"}
    flags=${flags:-""}
    model_flags=${model_flags:-" "}

    # Assign the values given by the user
    while [ $# -gt 0 ]; do
        if [[ $1 == *"--"* ]]; then
            param="${1/--/}"
            declare -g $param="$2"
        fi
        shift
    done
}

# get cli options
pararser $@

echo "Running your script now"
echo

# parse extra flags that evaluate to 'store_true'
IFS=',' read -ra ADDR <<< "$flags"
for flg in "${ADDR[@]}"; do
  extra_flags+=" --$flg" 
done

echo "Extra flags: $extra_flags"
echo "model_flags: $model_flags"

if [ "$model_flags" == "" ]; then
  model_flags_cmd=""
else
  model_flags_cmd="--model_flags $model_flags"
fi

python src/routines/main.py --env_type $env_type --separation $separation --model $model --initial_type $initial_type --environments $environments --numberTrials $numberTrials --numberSamples $numberSamples --dims $dims --selection $selection --compensation $compensation --intensity $intensity --mutation $mutation --savepath $savepath --parameter $parameter --intervals=$intervals --num_threads 4 --density $density --radius $radius --background $extra_flags $model_flags_cmd

date
