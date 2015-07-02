# This is the offline version of the script that runs the simulation of 
# the model for the game using parameters specified by the user. The results 
# of simulation are saved in results.pkl file, the plots produced from them 
# are saved in plot*.pdf files, and the mean trait value over the last 10 
# percent of the generations is written to STDOUT. The script is offline in 
# the sense that the plots are generated after the simulation is over.
# 
# This script can also be used to extend an existing experiment to include 
# more generations. That can be done by updating the "generations" parameter 
# in the experiment.py file for the experiment to the new value and running 
# the modified file. If a results.pkl file is present in the folder, the 
# script continues the experiment from where left off. Note that when 
# extending an experiment, it is very important that any graph_seed specified 
# in params map be the value from the original experiment, so that the network 
# structure (if any) in the extended experiment is the same as that in the 
# original experiment.

import math, numpy, os, pickle, pylab, random
import matplotlib.cm as cm
from matplotlib.colors import normalize as Norm
from continuous import *

def run(params, payoff):
    """
    Run simulation of the model for the game using the parameters specified 
    in params map, and the given payoff function. Once the simulation is 
    over call gen_plots() to generate plots from the results of simulation.  
    """

    # Extension of an experiment.
    if os.path.isfile("results.pkl"):
        fh = open("results.pkl", "r")
        old_params = pickle.load(fh) 
        population = pickle.load(fh)
        fh.close()
        start = old_params["generations"] + 1
        net = network.build_network(population, params["network_topology"], 
                                    params["network_params"]) 

        # Seed the random number generator for the experiment.
        random.seed(params["seed"])
    
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
    for time_step in range(start, params["generations"] + 1):

        # Pre interaction.
        dynamics.pre_interaction()

        # Save state at report_freq.
        if time_step % params["report_freq"] == 0:
            print("Generation %d of %d" %(time_step, params["generations"]))
            for p in population:
                p.save()

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

    # Pickle (save) the results of the experiment.
    output = open("results.pkl", "w")
    pickle.dump(params, output)
    pickle.dump(population, output)
    output.close()

    # Generate plots.
    gen_plots()

def gen_plots():
    """
    Create plots and save them in plot*.pdf files, and write the mean trait 
    value over the last 10 percent of the generations to STDOUT.
    """

    # Open and load the stuff from the results file.
    fh = open("results.pkl", "r")
    params = pickle.load(fh)
    population = pickle.load(fh)
    fh.close()

    # generations axis ticks.
    y = [i * params["report_freq"] \
             for i in range(0, params["generations"] / \
                                params["report_freq"] + 1)]
         
    # Plot 1 (mean trait vs generations)
    mean = []
    for i in range(0, len(y)):
        l = []
        for p in population:
            trait = p.get_trait_list()[i]
            payoff = p.get_payoff_list()[i]
            l.append(trait)
        mean.append(numpy.average(l))
    last = int(params["generations"] * 0.1)
    i = (params["generations"] / params["report_freq"]) - \
        (last / params["report_freq"])
    final_trait_mean = numpy.average(mean[i:])
    pylab.figure(1, figsize = (7, 4.5), dpi = 500)
    pylab.xlabel(r"$t$")
    pylab.ylabel(r"$\bar{x}$")
    pylab.plot(y, mean, "#000000", alpha = 0.6, linewidth = 2.0)
    pylab.figtext(0.82, 0.85, r"$\bar{x}_\infty = %4.3f$" %(final_trait_mean), 
                  ha = 'center', va = 'center', 
                  bbox = dict(facecolor = 'white', edgecolor = 'black'))
    pylab.xlim(0, params["generations"])
    pylab.ylim(0, params["max_trait"])
    ax = pylab.gca()
    ax.xaxis.major.formatter.set_powerlimits((0,0))
    pylab.savefig("plot1.pdf", format = "pdf")
    pylab.close(1)

    # Plot 2 (trait distribution)
    m = []
    vmin = 1.0
    vmax = 0.0
    for i in range(0, len(y)):
        l = []
        for p in population:
            l.append(p.get_trait_list()[i])
        l = bins(l, 0, params["max_trait"], 100)
        r = []
        for bin in l:
            trait = bin / (1.0 * params["population"])
            c = 1.0 - trait
            if c < vmin and c > 0:
                vmin = c
            if c > vmax and c < 1:
                vmax = c
            r.append(c)
        m.append(r)
    pylab.figure(1, figsize = (7, 4.5), dpi = 500)
    pylab.xlabel(r"$x$")
    pylab.ylabel(r"$t$")
    pylab.imshow(m, interpolation = "bilinear", origin = "l", 
                 cmap = cm.gray, norm = Norm(vmin, vmax), 
                 extent = [0, params["max_trait"], 1, params["generations"]])
    pylab.axis("tight")
    ax = pylab.gca()
    ax.yaxis.major.formatter.set_powerlimits((0,0))
    pylab.savefig("plot2.pdf", format = "pdf")
    pylab.close(1)

