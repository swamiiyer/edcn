# This is the offline version of the script that runs the simulation of 
# the model for the game using parameters specified by the user. The results 
# of simulation are saved in results.pkl file, and the plots produced 
# from them are saved in plot.pdf file. The script is offline in the sense 
# that the plots are generated after the simulation is over.
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

import math, os, pickle, pylab, random
from discrete import *

def run(params):
    """
    Run simulation of the model for the game using the parameters specified 
    in params map. Once the simulation is over call gen_plots() to generate 
    plots from the results of simulation.
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
    
        # Start out with a population with the specified initial fraction of 
        # the population playing trait C and the remaining fraction 
        # playing trait D.
        for p in population:
            if random.random() < params["init_cooperators"]:
                p.inherit_trait(Player.C)
            else:
                p.inherit_trait(Player.D)
            p.commit_inheritance()

    # Create a dynamics module based on network (complete or other) 
    # type and the type of update rule selected. 
    dynamics = dynamics_module(params["network_topology"], 
                               params["update_rule"])(net, params)

    # The dynamics.
    for time_step in range(start, params["generations"] + 1):

        # Pre interaction.
        dynamics.pre_interaction()

        # Save state at report_freq and stop simulation if population 
        # has fixated to a trait.
        if time_step % params["report_freq"] == 0:
            print("Generation %d of %d" %(time_step, params["generations"]))
            c_count = 0
            d_count = 0
            for p in population:
                p.save()
                if p.get_trait() == Player.C:
                    c_count += 1
                else:
                    d_count += 1
            if c_count == 0 or d_count == 0:
                break

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
    Create plots and save them in plot.pdf file, and print the mean fraction 
    of cooperators over the last 10 percent of the generations to STDOUT.
    """

    # Open and load the stuff from the results file.
    fh = open("results.pkl", "r")
    params = pickle.load(fh)
    population = pickle.load(fh)
    fh.close()

    # generations axis ticks.
    x = [i * params["report_freq"] \
             for i in range(0, params["generations"] / \
                                params["report_freq"] + 1)]
    
    y1 = []
    fixation_time_in_ticks = len(population[0].get_trait_list()) - 1
    fixation_time = fixation_time_in_ticks * params["report_freq"]
    fixated_to_c = False
    for t in range(0, params["generations"] / params["report_freq"] + 1):
        freq_c = 0
        if t <= fixation_time_in_ticks:
            for p in population:
                if p.get_trait_list()[t] == Player.C:
                    freq_c += 1
            if freq_c == params["population"]:
                fixated_to_c = True
        else:
            if fixated_to_c:
                freq_c = params["population"]
        freq_c = 1.0 * freq_c / params["population"]
        y1.append(freq_c)
    last = int(params["generations"] * 0.1)
    i = (params["generations"] / params["report_freq"]) - \
        (last / params["report_freq"])
    l = y1[i:]
    avg_c = sum(l) / len(l) 

    # Plot 1 (generations versus frequency of C players).
    pylab.figure(1, figsize = (7, 4.5), dpi = 500)
    pylab.xlabel(r"$t$")
    pylab.ylabel(r"$x$")
    pylab.plot(x, y1, "#000000", alpha = 0.6, linewidth = 2.0)
    pylab.figtext(0.82, 0.85, r"$x_\infty = %4.3f$" %(avg_c), 
                  ha = 'center', va = 'center', 
                  bbox = dict(facecolor = 'white', edgecolor = 'black'))
    pylab.ylim(0.0, 1.0)
    ax = pylab.gca()
    ax.xaxis.major.formatter.set_powerlimits((0,0))
    pylab.savefig("plot.pdf", format = "pdf")
    pylab.close(1)
