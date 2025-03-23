"""This script unlike the 'rough front' all the simulation is run here. This simulation uses parallel update rules.
    Here the number of genetic crossings i.e. how many sites can compete for a given site can be entered directly into cline.
    For even # of 'gc' it runs using hexagonal lattice, for even #'gc' it uses square lattice to avoid bias. 
"""
#Please be sure to install following packages before running! 
from os import sep
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import random
from matplotlib.colors import ListedColormap
import tools
import argparse
import math
import statistics
from numpy.core.arrayprint import printoptions
from numpy.core.fromnumeric import mean
import pandas as pd
from pandas.core.indexes.base import Index
 

# paramters
parser = argparse.ArgumentParser()
parser.add_argument("--gc", type=int, default=2)
parser.add_argument("--b", type=float, default=0.14)
parser.add_argument("--s", type=float, default=0.3)
parser.add_argument("--mu", type=float, default=0.03)
parser.add_argument("--phase", action= "store_true") #if called: creates '.data' which can be plotted in analysis.py 
parser.add_argument("--split", action="store_true") #populates first half of colony with bystander population and second half with fast invader 
parser.add_argument("--survival", action="store_true") #if called: created '.data' file of survival probability of a single fast invader surviving at frontier 
parser.add_argument("--mix", action="store_true") #populates initial colony with an even mix of bystander and fast invader 
parser.add_argument("--plot", action="store_true") #creates a image of simulation at the end 
parser.add_argument("--boundary", action="store_true") #Traces out boundary and stores boundary values into a csv of x position and y time
parser.add_argument("--phi", action="store_true") 
parser.add_argument("--rho", action="store_true")
parser.add_argument("--lx", type= int, default=200)
parser.add_argument("--ly", type= int, default=200)
args = parser.parse_args()
###################################
#Sets directory path in which to save files, WILL HAVE TO CHANGE ACCORDING TO USER'S PATH. Simply change next line and adjust path. 
programDir = tools.SetDir(f"/Users/ernie/Documents/Github/Summer 2021")
outFiles = tools.SetDir(f"{programDir}/Data")
rates = [ 1 - args.s + args.b, 1.0 , 1.0 - args.s ] #0: bystander, 1:fast-invader, 2: slow-invader
gc = args.gc
nsteps = args.ly
nsites = args.lx


# create one-D array
colony = np.zeros(nsites, dtype=int)


# simple keeping track of sub-populations
populations = defaultdict(list)


# Creates lattice
final_lattice = []

####################################################################
"""re-initialize colony with any of the following inital populations"""
if args.mix:
    for x in range(nsites):
        colony[x] = random.choice([0,1])

if args.split:
    for x in range(nsites):
        colony[x] = 1 

    for i in range(nsites//2):
        colony[i] = 0 

if args.survival: 
    for x in range(nsites): 
        if x == (nsites//2): 
            colony[x] = 1 
        else: 
            colony[x] = 0 
####################################################################

"""Populations stores the site's id at the front, used to calculate rho_m, phi_f, ..."""
for idx in range(len(colony)):
    populations[colony[idx]].append(idx)


####################################################################
"""This is the heart of the simulation, it goes through each site for each row, picks which site competes for a spot in 
    'new' which is an empty array 1 above 'colony' array. It populates 'new' array piking a 'winner' from already populated 
    'colony' using a simplified G.A which picks based on the rates of the populations. 
"""
to_graph = np.full((nsteps+1, nsites), -1, dtype=int)

for time in range(1, nsteps):
    new = np.empty(colony.size, dtype=int)
    final_lattice.append(new)
    new_populations = defaultdict(list) 
    if (time%2 == 0):
        shifted_line = True
    else: 
        shifted_line = False
    
    for i in range(len(new)):
        competitors = []
        competitors_index = [] 
        if gc % 2 != 0:
            first_competitor_index = i - math.floor(gc / 2)
        elif gc % 2 == 0 and shifted_line == False:
            first_competitor_index =   i - math.floor((gc/2) - 1)
        else: 
            if gc % 2 == 0 and shifted_line == True: 
                first_competitor_index = i - math.floor((gc/2))
            else:
                raise Exception("First competitor index cannot be calcualted")
        

        for j in range(gc):
            current_competitor = first_competitor_index + j
            if current_competitor < 0 or current_competitor > len(new) - 1:
                competitors.append(0)
                competitors_index.append(current_competitor)
            else:
                competitors.append(rates[colony[current_competitor]])
                competitors_index.append(current_competitor)

        R = sum(competitors)
        eta = random.random() * R
        
        
        current_sum = 0
        picked = 0
        for f in range(len(competitors)): 
            current_sum += competitors[f]
            if eta <= current_sum:
                picked = competitors_index[f]
                break
            
   
        pickedID = colony[picked]
        if pickedID == 1 and random.random() < args.mu:
            pickedID = 2
        
        new[i] = pickedID
        new_populations[pickedID].append(i)
        
    populations = new_populations       
    colony = new
    to_graph[time, :] = colony



"""Creates Phase Diagram"""
if args.phase: 
    path = tools.SetDir(f"{outFiles}/Phase Diagram")
    name = f"phase_gc {args.gc}"
    rho_m = 1 - (len(new_populations[0]) ) / (sum(len(new_populations[i]) for i in range(3) )  )
    with open(f"{path}/{name}.data" , "a") as file: 
        file.write(f"{args.b} {args.mu} {rho_m}\n")


"""Calculates the fraction of fast-invader population"""
if args.phi:
    path = tools.SetDir(f"{outFiles}/Phi")
    name = f"phi_gc {args.gc}"
    phi_f = (len(new_populations[1]) ) / (sum(len(new_populations[i]) for i in range(1,3) ) )
    with open (f"{path}/{name}.data", "a") as file: 
        file.write(f"{args.b} {args.mu} {phi_f}\n")
    print("phi of f:", phi_f)

"""Calculates survival probability of a single fast-growing cell at the front"""
if args.survival: 
    path = tools.SetDir(f"{outFiles}/Survival")
    name = f"p_surv_gc {args.gc}"
    p_surv = (len(new_populations[1]) ) / (sum(len(new_populations[i]) for i in range(0,3) ) )
    print(p_surv)
    # with open (f"{path}/{name}.data", "a") as file: 
    #     file.write(f"{args.b}, {args.mu} {p_surv}\n")


"""Traces out boundary and stores it as a csv file"""
if args.boundary:
    import csv 
    for y in range(nsteps):
        for x in range(nsites): 
            x = x + (x % 2) * 0.5

    boundary = []
    for i in range(nsteps): 
        for j in range(nsites): 
            if to_graph[i,j] != 0:
                if j == 0 or j == nsites:
                    break
                else: 
                    boundary.append((i, j)) 
                    break
        
    output_array = np.array(boundary)
    np.savetxt("test.csv", output_array, delimiter=",")
    with open("test.csv", "a") as f:
        writer = csv.writer(f, delimiter = ",", lineterminator = "\n")
        writer.writerows(boundary) 


if args.plot:
    fig, ax = plt.subplots()
    colors = {
        -1: "white",
        0: "yellow",
        1: "red",
        2: "black"
    }

    colors = ListedColormap([ "white", "yellow", "red", "black"])
    ax.imshow(to_graph, origin = "lower", cmap = colors)
    plt.show()


