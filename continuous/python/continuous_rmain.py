# This is the realtime version of the script that runs the simulation of 
# the model for the game using parameters specified by the user. The results 
# of simulation are not saved in this case, but are rendered as plots on the 
# screen as the simulation is taking place.
# 
# This script can also be used to replay an already-run simulation. In this 
# case, the script loads the results.pkl file from the simulation and 
# renders the plots on the screen.

import math, os, pickle, pylab, random
import matplotlib.cm as cm
from continuous import *

def run(params, payoff):
    """
    Run simulation of the model for the game using the parameters specified 
    in params map, and the given payoff function. If a results.pkl file is 
    present in the current folder, then the script replays the simulation. 
    Otherwise, the simulation is run anew. 
    """

    # Canvas for drawing.
    pylab.figure(1, figsize = (12, 6))
    pylab.ion()
    pylab.draw()

    trait_matrix = []

    # Replay of an experiment.
    if os.path.isfile("results.pkl"):
        fh = open("results.pkl", "r")
        pickle.load(fh) 
        population = pickle.load(fh)
        fh.close()
        for time_step in range(0, params["generations"] / \
                                   params["report_freq"] + 1):
            print("Generation %d of %d" %(time_step * params["report_freq"], 
                                          params["generations"]))
            traits = []
            for p in population:
                traits.append(p.get_trait_list()[time_step])
            plot_trait_distribution(params, trait_matrix, traits)
            plot_trait_histogram(params, traits)
    
    # New experiment.
    else:
        start = 0
        
        # Create a population.
        population = [Player(i) for i in range(0, params["population"])]

        # Create a network.
        net = network.build_network(population, params["network_topology"], 
                                    params["network_params"])
    
        # Seed the random number generator for the experiment.
        random.seed(params["seed"])
    
        # Assign a trait for each individual in the population. 
        for p in population:
            p.inherit_trait(params["init_trait"])
            p.commit_inheritance()

        # Create a dynamics module based on network (complete or other) 
        # type and update rule selected. 
        dynamics = dynamics_module(params["network_topology"], 
                                   params["update_rule"])(net, params, payoff)

        # The dynamics.
        for time_step in range(0, params["generations"]):

            # Pre interaction.
            dynamics.pre_interaction()

            # Plot results at report_freq.
            if time_step % params["report_freq"] == 0:
                print("Generation %d of %d" %(time_step, params["generations"]))
                traits = []
                for p in population:
                    traits.append(p.get_trait())
                plot_trait_distribution(params, trait_matrix, traits)
                plot_trait_histogram(params, traits)

            # Interact.
            for count in range(0, params["population"]):
                dynamics.interact()

            # Post interaction.
            dynamics.post_interaction()
            
            # Update.
            for count in range(0, params["population"]):
                dynamics.update()

            # Post update; commit trait inheritance.
            dynamics.post_update()

    # Keep the final plot window around until the user shuts it.
    pylab.show()
    pylab.close(1)
    
def plot_trait_distribution(params, trait_matrix, traits):
    """
    Plot the time evolution of trait distribution given the params map 
    containing the experimental parameters, the trait matrix containing 
    history of traits, and current list of traits.
    """

    l = bins(traits, 0, params["max_trait"], 100)
    r = []
    for bin in l:
        trait = bin / (1.0 * params["population"])
        r.append(1.0 - trait)
    trait_matrix.append(r)
    pylab.subplot(121).clear()
    pylab.xlabel(r"$x$")
    pylab.ylabel(r"$t$")
    pylab.imshow(trait_matrix, interpolation = "bilinear", 
                 origin = "l", cmap = cm.gray, 
                 extent = [0, params["max_trait"], 1, \
                               len(trait_matrix) * params["report_freq"]])
    pylab.axis("tight")
    ax = pylab.gca()
    ax.yaxis.major.formatter.set_powerlimits((0,0))
    pylab.draw()

def plot_trait_histogram(params, traits):
    """
    Plot a histogram of number of individuals versus trait values given 
    the params map containg the experimental parameters, and the list of traits 
    for a particular time step.
    """

    pylab.subplot(122).clear()
    pylab.xlabel(r"$x$")
    pylab.ylabel(r"$I$")
    pylab.hist(traits, 100, range = (0, params["max_trait"]), normed = True, 
               facecolor = "black")
    pylab.xlim(0, params["max_trait"])
    pylab.ylim(0, params["population"])
    ax = pylab.gca()
    ax.yaxis.major.formatter.set_powerlimits((0,0))
    pylab.draw()
