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
from disease import *

def run(params):
    """
    Run simulation of the model for the game using the parameters specified 
    in params map. If a results.pkl file is present in the current folder, 
    then the script replays the simulation. Otherwise, the simulation is 
    run anew. 
    """
    
    # Canvas for drawing.
    pylab.figure(1, figsize = (12, 6))
    pylab.ion()
    pylab.draw()

    s_list = []
    i_list = []
    r_list = []
    mean_c_list = []
    time_steps = []
    
    # Replay of an experiment.
    if os.path.isfile("results.pkl"):
        fh = open("results.pkl", "r")
        old_params = pickle.load(fh) 
        population = pickle.load(fh)
        fh.close()
        for time_step in range(0, params["generations"] / \
                                   params["report_freq"] + 1):
            print("Generation %d of %d" %(time_step * params["report_freq"], 
                                          params["generations"]))
            time_steps.append(time_step * params["report_freq"])
            c_list = []
            s_count = 0
            i_count = 0
            r_count = 0
            v_count = 0
            for p in population:
                state = p.get_state_list()[time_step]
                if state == Player.SUSCEPTIBLE:
                    s_count += 1
                elif state == Player.INFECTED:
                    i_count += 1
                    c_list.append(p.get_c_list()[time_step])
                elif state == Player.RECOVERED:
                    r_count += 1
                else: # vaccinated
                    v_count += 1
            s_list.append(s_count / (1.0 * params["population"]))
            i_list.append(i_count / (1.0 * params["population"]))
            r_list.append(r_count / (1.0 * params["population"]))
            if (i_count == 0):
                mean_c_list.append(0.0)
            else:
                mean_c_list.append(sum(c_list) / i_count)
            plot_sir_curves(time_steps, s_list, i_list, r_list)
            plot_mean_virulence(time_steps, mean_c_list)
    
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
        # Pick a person to infect. This person is a random neighbor of a 
        # person picked at random.
        r = net.get_random_vertex()
        p = net.get_random_neighbor(r)
        while p == None:
            r = net.get_random_vertex()
            p = net.get_random_neighbor(r)
        p.set_state(Player.INFECTED)
        p.set_infection(params["infection"])
        p.set_c(fitness(params["infection"], params["b"],  params["a"]))

        # Disease dynamics.
        for time_step in range(start, params["generations"] + 1):

            if time_step % params["report_freq"] == 0:
                print("Generation %d of %d" %(time_step, params["generations"]))
                c_list = []
                s_count = 0
                i_count = 0
                r_count = 0
                v_count = 0
                time_steps.append(time_step)
                for p in population:
                    state = p.get_state()
                    if state == Player.SUSCEPTIBLE:
                        s_count += 1
                    elif state == Player.INFECTED:
                        i_count += 1
                        c_list.append(p.get_c())
                    elif state == Player.RECOVERED:
                        r_count += 1
                    else: # vaccinated
                        v_count += 1
                s_list.append(s_count / (1.0 * params["population"]))
                i_list.append(i_count / (1.0 * params["population"]))
                r_list.append(r_count / (1.0 * params["population"]))
                if (i_count == 0):
                    mean_c_list.append(0.0)
                else:
                    mean_c_list.append(sum(c_list) / i_count)
                plot_sir_curves(time_steps, s_list, i_list, r_list)
                plot_mean_virulence(time_steps, mean_c_list)

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
                        if random.random() < prob:
                            r = random.random() * total_infection
                            total = 0.0
                            for vertex in neighbors:
                                new_infection = vertex.get_infection()
                                total += new_infection
                                if r < total:
                                    break

                            # The pathogen causing the disease could mutate.
                            if random.random() < params["mutation"]:
                                new_infection = \
                                    max(0.0, min(random.gauss(new_infection, 
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

                    # c can die of natural causes. If c does die, then c is 
                    # reborn as a susceptible individual.
                    if random.random() < params["death"]:
                        c.set_state(Player.SUSCEPTIBLE)
                        c.set_infection(0.0)
                        c.set_c(0.0)

                # c is vaccinated, in which case nothing happens, so we 
                # continue with the next individual.
                else:
                    continue

    # Keep the final plot window around until the user shuts it.
    pylab.show()
    pylab.close(1)
        
def plot_sir_curves(time_steps, s_list, i_list, r_list):
    """
    Plot the SIR curves given the list of x values 
    specifying the generations, and the lists of y values specifying the 
    number of susceptible, infected and recovered individuals.
    """

    pylab.subplot(121).clear()
    pylab.xlabel(r"$t$")
    pylab.ylabel(r"$s, x, r$")
    pylab.plot(time_steps, s_list, "#000000", alpha = 0.6, linewidth = 2.0)
    pylab.plot(time_steps, i_list, "#FF0000", alpha = 0.6, linewidth = 2.0)
    pylab.plot(time_steps, r_list, "#0000FF", alpha = 0.6, linewidth = 2.0)
    pylab.legend((r"$s$", r"$x$", r"$r$"), 'best', shadow = False)
    pylab.ylim(0, 1)
    ax = pylab.gca()
    ax.xaxis.major.formatter.set_powerlimits((0,0))
    pylab.draw()

def plot_mean_virulence(time_steps, mean_c_list):
    """
    Plot the time evolution of mean virulence given given the list 
    of x values specifying the generations, and the lists of y values 
    specifying the mean virulence levels.
    """

    pylab.subplot(122).clear()
    pylab.xlabel(r"$t$")
    pylab.ylabel(r"$\bar{c}$")
    pylab.plot(time_steps, mean_c_list, "#000000", alpha = 0.6, linewidth = 2.0)
    pylab.ylim(0, 1.0)
    ax = pylab.gca()
    ax.xaxis.major.formatter.set_powerlimits((0,0))
    pylab.draw()

