# This is the realtime version of the script that runs the simulation of 
# the model for the game using parameters specified by the user. The results 
# of simulation are not saved in this case, but are rendered as plots on the 
# screen as the simulation is taking place.
# 
# This script can also be used to replay an already-run simulation. In this 
# case, the script loads the results.pkl file from the simulation and 
# renders the plots on the screen.

import math, os, pickle, pylab, random
from discrete import *

def run(params):
    """
    Run simulation of the model for the game using the parameters specified 
    in params map. If a results.pkl file is present in the current folder, 
    then the script replays the simulation. Otherwise, the simulation 
    is run anew. 
    """

    # Canvas for drawing.
    pylab.figure(1, figsize = (7, 4.5))
    pylab.ion()
    pylab.draw()

    c_list = []
    time_steps = []
    
    # Replay of an experiment.
    if os.path.isfile("results.pkl"):
        fh = open("results.pkl", "r")
        pickle.load(fh) 
        population = pickle.load(fh)
        fh.close()
        done = False
        for time_step in range(0, params["generations"] / \
                                   params["report_freq"] + 1):
            print("Generation %d of %d" %(time_step * params["report_freq"], 
                                          params["generations"]))
            time_steps.append(time_step * params["report_freq"])
            c_count = 0
            d_count = 0
            for p in population:
                if time_step >= len(p.get_trait_list()):
                    done = True
                else: 
                    if p.get_trait_list()[time_step] == Player.C:
                        c_count += 1
                    else:
                        d_count += 1
            if done == True:
                break
            c_list.append(1.0 * c_count / params["population"])
            plot(time_steps, c_list)

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
        for time_step in range(0, params["generations"]):

            # Pre interaction.
            dynamics.pre_interaction()

            # Plot results at report_freq and stop simulation if 
            # population has fixated to a trait.
            if time_step % params["report_freq"] == 0:
                print "Generation %d of %d" %(time_step, params["generations"])
                time_steps.append(time_step)
                c_count = 0
                d_count = 0
                for p in population:
                    if p.get_trait() == Player.C:
                        c_count += 1
                    else:
                        d_count += 1
                c_list.append(1.0 * c_count / params["population"])
                plot_c_d_curves(time_steps, c_list)
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

    # Keep the final plot window around until the user shuts it.
    pylab.show()
    pylab.close(1)
    
def plot(time_steps, c_list):
    """
    Plot the frequency of individuals playing the C strategy.
    """

    pylab.subplot(111).clear()
    pylab.xlabel(r"$t$")
    pylab.ylabel(r"$x$")
    pylab.plot(time_steps, c_list, "#000000", alpha = 0.6, linewidth = 2.0)
    pylab.ylim(0, 1)
    ax = pylab.gca()
    ax.xaxis.major.formatter.set_powerlimits((0,0))
    pylab.draw()
