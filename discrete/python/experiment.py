# This script can be used to simulate discrete 2x2 games by setting the model 
# parameters to the desired values. 
#
# Usage: python experiment.py [--realtime]
# 
# where the optional --realtime argument runs the simulation in real-time mode 
# instead of the default off-line mode.
# 
# In the off-line mode, the results of simulation are saved in results.pkl 
# file, and a plot produced from them is saved in plot.pdf file. 
# In the real-time mode, results of the simulation are plotted in real time, 
# and are not saved to a file. 
#
# Note that the off-line version can also be used to extend a previously 
# run experiment to include more generations. Also, the real-time version 
# can be used to visualize the results saved from running the off-line version.
# Check discrete_omain.py and discrete_rmain.py for documentation of these 
# features.

import discrete, discrete_omain, discrete_rmain, sys, time

def main(args):
    """
    Entry point.
    """

    params = {
        # A 2x2 payoff matrix.
        "payoff" : None, 

        # Selection strength.
        "selection_strength" : None, 

        # Total number of players.
        "population" : None, 

        # Initial fraction of population playing the cooperative trait C.
        "init_cooperators" : None, 
        
        # Number of generations.
        "generations" : None,

        # Frequency (in generations) at which results will be saved.
        "report_freq" : None, 

        # The probability with which a player interacts with another player 
        # having the same trait as itself. This setting is relevant only for 
        # complete networks.
        "assortativity" : None, 

        # Update rule (RE1, RE2, FE1, FE2, IM, BD, DB).
        "update_rule" : None, 

        # Random number seed; default value is milliseconds since midnight 
        # Jan 1, 1970.
        "seed" : int(time.time()), 

        # Network topology and parameters. See build_network() in 
        # edcn/util/network.py.
        "network_topology" : None, 
        "network_params" : {}
        }

    if len(args) == 1 and args[0] == "--realtime":
        discrete_rmain.run(params)
    else:
        discrete_omain.run(params)

if __name__ == "__main__":
    main(sys.argv[1:])
