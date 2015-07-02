# This script can be used to simulate continuous pairwise games by setting 
# the model parameters to the desired values. 
#
# Usage: python experiment.py [--realtime]
# 
# where the optional --realtime argument runs the simulation in real-time mode 
# instead of the default off-line mode.
# 
# In the off-line mode, the results of simulation are saved in results.pkl 
# file, the plots produced from them are saved in plot*.pdf files. 
# In the real-time mode, results of the simulation are plotted in real time, 
# and are not saved to a file. 
#
# Note that the off-line version can also be used to extend a previously 
# run experiment to include more generations. Also, the real-time version 
# can be used to visualize the results saved from running the off-line version.
# Check continuous_omain.py and continuous_rmain.py for documentation of 
# these features.

import continuous, continuous_omain, continuous_rmain, sys, time

def main(args):
    """
    Entry point. 
    """

    # The payoff function.
    def B(x):
        pass
    def C(x):
        pass
    # payoff = lambda x, y: B(y) - C(x) # prisoner's dilemma
    # payoff = lambda x, y: B(x + y) - C(x) # snowdrift game
    # payoff = lambda x, y: B(x) - C(x + y) # tragedy of the commons

    params = {
        # Selection strength.
        "selection_strength" : None, 

        # Total number of players.
        "population" : None, 
        
        # Number of generations.
        "generations" : None,

        # Frequency (in generations) at which results will be saved.
        "report_freq" : None, 

        # Initial trait for each player.
        "init_trait" : None, 

        # Maximum value that the trait of any player can attain.
        "max_trait" : None, 

        # The probability with which a player interacts with another player 
        # having the same trait as itself. This setting is relevant only for 
        # complete networks.
        "assortativity" : None, 

        # group specifies the number of neighbors each player 
        # interacts with during each round of interactions. If None, each 
        # player will interact with its entire neighborhood. The average of 
        # the neighbors' traits is used in calculating the payoff. Note 
        # that this setting is not applicable in "FE2" and "RE2" dynamics,  
        # and must be set to 1 if "assortativity" is positive.
        "group" : None, 

        # Update rule (RE1, RE2, FE1, FE2, IM, BD, DB).
        "update_rule" : None, 

        # Rate at which Gaussian mutations of traits are carried out.
        "mutation" : None, 

        # Standard deviation for the Gaussian mutations.
        "stddev" : None,  

        # Random number seed; default value is milliseconds since midnight 
        # Jan 1, 1970.
        "seed" : int(time.time()), 

        # Network topology and parameters. See build_network() in 
        # edcn/util/network.py.
        "network_topology" : None, 
        "network_params" : {}
        }
    
    if len(args) == 1 and args[0] == "--realtime":
        continuous_rmain.run(params, payoff)
    else:
        continuous_omain.run(params, payoff)

if __name__ == "__main__":
    main(sys.argv[1:])
