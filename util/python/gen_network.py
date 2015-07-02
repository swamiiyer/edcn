#!/usr/bin/python

# Generates a GraphML file representing the network with the specified 
# topology and properties.
#
# Usage: python gen_network.py <topology> <size> <params> <outfile>
#
# The allowed topologies are:
#
#      Complete, Erdos_Renyi, Barabasi_Albert, Powerlaw_Homophilic, 
#      Powerlaw_Clustered, Random_Regular, Grid_2D, KNN, 
#      Watts_Strogatz
#
# The params argument is a map specifying the parameters for the 
# requested network topology. Allowed parameters are:
#
#      Complete: {}
#      Erdos_Renyi: {"p": ..., "seed": ...}
#      Barabasi_Albert: {"m": ..., "seed": ...}
#      Powerlaw_Homophilic: {"m": ..., "h": ..., "maxiter": ..., 
#                             "seed": ...}
#      Powerlaw_Clustered: {"m": ..., "p": ..., "seed": ...}
#      Random_Regular: {"degree": ..., "seed": ...}
#      Grid_2D: {"m": ..., "n": ..., "periodic": ...}
#      KNN: {"n": ..., "k": ..., "periodic": ...}
#      Watts_Strogatz: {"k": ..., "p": ..., "seed": ...}

import networkx, sys
from network import *

def main(args):
    if len(args) != 4:
        print \
            "Usage: python gen_network.py <topology> <size> <params> <outfile>"
        sys.exit()
    topology = args[0]
    size = int(args[1])
    params = eval(args[2])
    outfile = args[3]
    vertices = [Vertex(i) for i in range(0, size)]
    net = build_network(vertices, topology, params)
    g = net.get_structure()
    networkx.write_graphml(g, outfile)

if __name__ == "__main__":
    main(sys.argv[1:])
