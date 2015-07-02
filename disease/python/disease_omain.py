# This is the offline version of the script that runs the simulation of 
# the model for the game using parameters specified by the user. The results 
# of simulation are saved in results.pkl file, the plots produced from them 
# are saved in plot*.pdf files, and the mean virulence value over the last 10 
# percent generations is printed to STDOUT. The script is offline in the 
# sense that the plots are generated after the simulation is over.
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
import matplotlib.cm as cm
from matplotlib.colors import normalize as Norm
from disease import *

def run(params):
    """
    Run simulation of the model for the game using the parameters specified 
    in params map. Once the simulation is over call gen_plots() to
    generate plots from the results of simulation.
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

        # Pick a person to infect. This person is a random neighbor of a 
        # person picked at random.
        r = net.get_random_vertex()
        p = net.get_random_neighbor(r)
        while p == None:
            r = net.get_random_vertex()
            p = net.get_random_neighbor(r)
        p.set_state(Player.INFECTED)
        p.set_infection(params["infection"])
        p.set_c(fitness(params["infection"], params["b"], params["a"]))

    # Disease dynamics.
    for time_step in range(start, params["generations"] + 1):

        # Save state at report_freq.
        if time_step % params["report_freq"] == 0:
            print("Generation %d of %d" %(time_step, params["generations"]))
            for p in population:
                p.save()

        for count in range(0, params["population"]):
            # Pick an individual c from the population at random.
            c = net.get_random_vertex()

            # c is a susceptible individual.
            if c.get_state() == Player.SUSCEPTIBLE:

                # c can die of natural causes. If c dies, then c is reborn 
                # as a susceptible individual; nothing changes here, so we 
                # simply continue with the next individual.
                if random.random() < params["death"]:
                    continue

                # If c does not die, then c can become infected by any 
                # infected neighbor c may have. 
                else:
                    prob_not = 1.0 # probability that c is not infected
                    total_infection = 0.0
                    neighbors = net.get_neighbors(c)
                    if neighbors == None:
                        continue
                    for vertex in neighbors:
                        if vertex.get_state() == Player.INFECTED:
                            prob_not = prob_not  * \
                                (1.0 - vertex.get_infection())
                            total_infection += vertex.get_infection()
                    prob = 1 - prob_not # probability that c is infected

                    # If c does get infected then c takes the transmission 
                    # rate of one of the infected neighbors picked in 
                    # proportion to that individual's transmission rate.
                    new_infection = 0.0
                    if random.random() < prob:
                        r = random.random() * total_infection
                        total = 0.0
                        for vertex in neighbors:
                            if vertex.get_state() == Player.INFECTED:
                                new_infection = vertex.get_infection()
                                total += new_infection
                                if r < total:
                                    break

                        # The pathogen causing the disease could mutate.
                        if random.random() < params["mutation"]:
                            new_infection = max(0.0, 
                                           min(random.gauss(new_infection, 
                                                            params["stddev"]), 
                                               1.0))

                        c.set_state(Player.INFECTED)
                        c.set_infection(new_infection)
                        c.set_c(fitness(new_infection, params["b"], 
                                        params["a"]))

            # c is an infected individual.
            elif c.get_state() == Player.INFECTED:

                # c can die of natural causes or from the disease. If c  
                # dies, then c is reborn as a susceptible individual.
                death_rate = 1 - (1 - params["death"]) * (1 - c.get_c())
                if random.random() < death_rate:
                    c.set_state(Player.SUSCEPTIBLE)
                    c.set_infection(0.0)
                    c.set_c(0.0)

                # If c does not die, then c can recover with the recovery 
                # rate.
                else:
                    if random.random() < params["recovery"]:
                        c.set_state(Player.RECOVERED)
                        c.set_infection(0.0)
                        c.set_c(0.0)

            # c is a recovered individual.
            elif c.get_state() == Player.RECOVERED:

                # c can die of natural causes. If c does die, then c is reborn 
                # as a susceptible individual.
                if random.random() < params["death"]:
                    c.set_state(Player.SUSCEPTIBLE)
                    c.set_infection(0.0)
                    c.set_c(0.0)

            # c is vaccinated, in which case nothing happens, so we continue 
            # with the next individual.
            else:
                continue
        
    # Pickle (save) the results of the experiment.
    output = open("results.pkl", "w")
    pickle.dump(params, output)
    pickle.dump(population, output)
    output.close()

    # Generate plots.
    gen_plots()

def gen_plots():
    """
    Create plots and save them in plots*.pdf file, and print the mean 
    virulence value over the last 10 percent generations to STDOUT.
    """

    # Open and load the stuff from the results file.
    fh = open("results.pkl", "r")
    params = pickle.load(fh)
    population = pickle.load(fh)
    fh.close()

    # x (generations) axis ticks.
    x = [i * params["report_freq"] \
             for i in range(0, params["generations"] / \
                                params["report_freq"] + 1)]

    # Calculate number of susceptible, infected, recovered invididuals, 
    # and mean disease-induced death rate per generation.
    y1 = [] # number of susceptible individuals
    y2 = [] # number of infected individuals
    y3 = [] # number of recovered individuals
    y4 = [] # mean disease-induced death rate
    for i in range(0, len(x)): 
        freq_s = 0
        freq_i = 0
        freq_r = 0
        freq_v = 0
        c_sum = 0.0
        c_mean = 0.0
        for p in population:
            if p.get_state_list()[i] == Player.INFECTED:
                freq_i += 1
                c_sum += p.get_c_list()[i]
            elif p.get_state_list()[i] == Player.RECOVERED:
                freq_r += 1
            elif p.get_state_list()[i] == Player.VACCINATED:
                freq_s += 1
            else:
                freq_v += 1
        if freq_i != 0:
            c_mean = c_sum / freq_i
        y1.append(freq_s / (1.0 * params["population"]))
        y2.append(freq_i / (1.0 * params["population"]))
        y3.append(freq_r / (1.0 * params["population"]))
        y4.append(c_mean)

    # Calculate mean virulence over the last 10% of the generations.
    last = int(params["generations"] * 0.1)
    i = (params["generations"] / params["report_freq"]) - \
        (last / params["report_freq"])
    l = y4[i:]
    final_c_mean = sum(l) / len(l) 

    # Plot 1 (number of susceptible, infected, and recovered individuals 
    # against generation).
    pylab.figure(1, figsize = (7, 4.5), dpi = 500)
    pylab.xlabel(r"$t$")
    pylab.ylabel(r"$s, x, r$")
    pylab.plot(x, y1, "#000000", alpha = 0.6, linewidth = 2.0)
    pylab.plot(x, y2, "#FF0000", alpha = 0.6, linewidth = 2.0)
    pylab.plot(x, y3, "#0000FF", alpha = 0.6, linewidth = 2.0)
    pylab.legend((r"$s$", r"$x$", r"$r$"), 'best', shadow = False)
    pylab.ylim(0, 1.0)
    ax = pylab.gca()
    ax.xaxis.major.formatter.set_powerlimits((0,0))
    pylab.savefig("plot1.pdf", format = "pdf")
    pylab.close(1)

    # Plot 2 (Mean disease-induced death rate versus generation).
    pylab.figure(2, figsize = (7, 4.5), dpi = 500)
    pylab.xlabel(r"$t$")
    pylab.ylabel(r"$\bar{c}$")
    pylab.plot(x, y4, "#000000", alpha = 0.6, linewidth = 2.0)
    pylab.figtext(0.82, 0.85, r"$\bar{c}_\infty = %4.3f$" %(final_c_mean), 
                  ha = 'center', va = 'center', bbox = dict(facecolor = 'white', edgecolor = 'black'))
    pylab.ylim(0, 1.0)
    ax = pylab.gca()
    ax.xaxis.major.formatter.set_powerlimits((0,0))
    pylab.savefig("plot2.pdf", format = "pdf")
    pylab.close(2)

