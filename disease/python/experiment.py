# This script can be used to simulate virus dynamics on networks by setting 
# the model parameters to the desired values. 
#
# Usage: python experiment.py [--realtime]
# 
# where the optional --realtime argument runs the simulation in real-time mode 
# instead of the default off-line mode.
# 
# In the off-line mode, the results of simulation are saved in results.pkl 
# file, and the plots produced from them are saved in plot*.pdf file.
# In the real-time mode, results of the simulation are plotted in real-time, 
# and are not saved to a file. 
#
# Note that the off-line version can also be used to extend a previously 
# run experiment to include more generations. Also, the real-time version 
# can be used to visualize the results saved from running the off-line version.
# Check disease_omain.py and disease_rmain.py for documentation of these 
# features.

import disease, disease_omain, disease_rmain, sys, time

def main(args):
    """
    Entry point. 
    """

    params = {
        # Total number of players.
        "population" : None, 
        
        # Number of generations.
        "generations" : None,

        # Frequency (in generations) at which results will be saved.
        "report_freq" : None, 

        # Infection probability for the starting infected individual. 
        # This may vary for different individuals in the infected state.
        "infection" : None,

        # Virulence level c of a pathogen is calculated from its infection 
        # probability as c = pow(infection / b, 1.0 / a).
        "b" : 0.0, 
        "a" : 0.0, 

        # Recovery probability.
        "recovery" : None,        

        # Natural death probability.
        "death" : None,

        # Rate at which Gaussian mutations of infection probability is 
        # carried out.
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
        disease_rmain.run(params)
    else:
        disease_omain.run(params)

if __name__ == "__main__":
    main(sys.argv[1:])
