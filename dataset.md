python src/routines/generate_parameter_space_search.py --numberTrials 100 --dims 1000,1000 --savepath workspace/experiments/phase_diagram_no_mutations_phi_0.25 --mutation 0.000 --selection 0.0 --intensity 10 --radius 10 --density 0.25 --overwrite --parameters selection,intensity --intervals_1 0,0.01,0.1 --intervals_2 0,0.5,6 --model src/base/ --environments 20 --standing_variation

python src/routines/generate_parameter_space_search.py --numberTrials 100 --dims 2000,2000 --savepath workspace/experiments/phase_diagram_no_mutations_phi_0.25_doublesize --mutation 0.000 --selection 0.0 --intensity 10 --radius 10 --density 0.25 --overwrite --parameters selection,intensity --intervals_1 0,0.01,0.1 --intervals_2 0,0.5,6 --model src/base/ --environments 20 --standing_variation

python src/routines/generate_parameter_space_search.py --numberTrials 100 --dims 2000,2000 --savepath workspace/experiments/phase_diagram_no_mutations_phi_0.25_doublesize_3R --mutation 0.000 --selection 0.0 --intensity 10 --radius 30 --density 0.25 --overwrite --parameters selection,intensity --intervals_1 0,0.01,0.1 --intervals_2 0,0.5,6 --model src/base/ --environments 20 --standing_variation

python src/routines/generate_parameter_space_search.py --numberTrials 100 --dims 1000,2000 --savepath workspace/experiments/phase_diagram_no_mutations_phi_0.25_twice_height --mutation 0.000 --selection 0.0 --intensity 10 --radius 10 --density 0.25 --overwrite --parameters selection,intensity --intervals_1 0,0.01,0.1 --intervals_2 0,0.5,6 --model src/base/ --environments 20 --standing_variation

python src/routines/generate_parameter_space_search.py --numberTrials 100 --dims 1000,2000 --savepath workspace/experiments/phase_diagram_no_mutations_phi_0.25_twice_scales --mutation 0.000 --selection 0.0 --intensity 10 --radius 20 --density 0.25 --overwrite --parameters selection,intensity --intervals_1 0,0.01,0.1 --intervals_2 0,0.5,6 --model src/base/ --environments 20 --standing_variation

python src/routines/generate_parameter_space_search.py --numberTrials 200 --dims 1000,1000 --savepath workspace/experiments/phase_diagram_no_mutations_phi_0p1 --mutation 0.000 --selection 0.0 --intensity 10 --radius 10 --density 0.1 --overwrite --parameters selection,intensity --intervals_1 0,0.01,0.1 --intervals_2 0,0.5,6 --model src/base/ --environments 20 --standing_variation

python src/routines/generate_parameter_space_search.py --numberTrials 200 --dims 2000,2000 --savepath workspace/experiments/phase_diagram_no_mutations_phi_0p1_doubled --mutation 0.000 --selection 0.0 --intensity 10 --radius 10 --density 0.1 --overwrite --parameters selection,intensity --intervals_1 0,0.01,0.1 --intervals_2 0,0.5,6 --model src/base/ --environments 20 --standing_variation

python src/routines/generate_parameter_space_search.py --numberTrials 200 --dims 1000,1000 --savepath workspace/experiments/phase_diagram_no_mutations_phi_0p2 --mutation 0.000 --selection 0.0 --intensity 10 --radius 10 --density 0.2 --overwrite --parameters selection,intensity --intervals_1 0,0.01,0.1 --intervals_2 0,0.5,6 --model src/base/ --environments 20 --standing_variation

python src/routines/generate_parameter_space_search.py --numberTrials 200 --dims 1000,1000 --savepath workspace/experiments/phase_diagram_no_mutations_phi_0p15 --mutation 0.000 --selection 0.0 --intensity 10 --radius 10 --density 0.15 --overwrite --parameters selection,intensity --intervals_1 0,0.01,0.1 --intervals_2 0,0.5,6 --model src/base/ --environments 20 --standing_variation

python src/routines/generate_parameter_space_search.py --numberTrials 200 --dims 1000,2000 --savepath workspace/experiments/phase_diagram_no_mutations_twice_height --mutation 0.000 --selection 0.0 --intensity 10 --radius 10 --density 0.1 --overwrite --parameters selection,intensity --intervals_1 0,0.01,0.1 --intervals_2 0,0.5,6 --model src/base/ --environments 20 --standing_variation

python src/routines/generate_parameter_space_search.py --numberTrials 200 --dims 1000,2000 --savepath workspace/experiments/phase_diagram_no_mutations_twice_height --mutation 0.000 --selection 0.0 --intensity 10 --radius 10 --density 0.1 --overwrite --parameters selection,intensity --intervals_1 0,0.01,0.1 --intervals_2 0,0.5,6 --model src/base/ --environments 20 --standing_variation

python src/routines/generate_parameter_space_search.py --numberTrials 100 --dims 1000,1000 --savepath workspace/experiments/size_distributions --mutation 0.001 --selection 0.0 --intensity 10 --radius 10 --density 0.1 --overwrite --parameters selection,intensity --intervals_1 0,0.01,0.1 --intervals_2 0,0.5,6 --model src/base/ --environments 20

python src/routines/generate_parameter_space_search.py --numberTrials 100 --dims 2000,2000 --savepath workspace/experiments/size_distributions_2000x2000 --mutation 0.001 --selection 0.0 --intensity 10 --radius 10 --density 0.1 --overwrite --parameters selection,intensity --intervals_1 0,0.01,0.1 --intervals_2 0,0.5,6 --model src/base/ --environments 20

python src/routines/generate_parameter_space_search.py --numberTrials 100 --dims 1000,1000 --savepath workspace/experiments/phase_diagram_no_mutations --mutation 0.000 --selection 0.0 --intensity 10 --radius 10 --density 0.08 --overwrite --parameters selection,intensity --intervals_1 0,0.01,0.1 --intervals_2 0,0.5,6 --model src/base/ --environments 20 --standing_variation

python src/routines/generate_parameter_space_search.py --numberTrials 100 --dims 1000,1000 --savepath workspace/experiments/percolation/percolation_limits --mutation 0.000 --selection 0.0 --intensity 10 --radius 10 --density 0.08 --overwrite --parameters selection,density --intervals_1 0,0.01,0.1 --intervals_2 0.005,0.005,0.04 --model src/base/ --environments 25 --standing_variation

python src/routines/generate_parameter_space_search.py --numberTrials 100 --dims 1000,1000 --savepath workspace/experiments/percolation/percolation_limits --mutation 0.000 --selection 0.0 --intensity 10 --radius 10 --density 0.08 --overwrite --parameters selection,density --intervals_1 0,0.01,0.1 --intervals_2 0,0.05,0.9 --model src/base/ --environments 25 --standing_variation

python src/routines/generate_parameter_space_search.py --numberTrials 200 --dims 1000,1000 --savepath workspace/experiments/phase_diaggram_mutations_phi_0.25 --mutation 0.0005 --selection 0.0 --intensity 10 --radius 10 --density 0.25 --overwrite --parameters selection,intensity --intervals_1 0,0.01,0.1 --intervals_2 0,0.5,6 --model src/base/ --environments 20
