#!/usr/bin/python
#
# This script prints as output the values of some of the basic properties of 
# the input network in GraphML format.
#
# Usage: python network_stats.py <infile>

import networkx, numpy, sys

def main(args):
    if len(args) == 0:
        print "Usage: python network_stats.py <infile>"
        sys.exit()
    infile = args[0]
    g = networkx.read_graphml(infile)
    n = len(g.nodes())
    m = len(g.edges())
    print "number of vertices = %d" %(n)
    print "number of edges = %d" %(m)
    print "connectance = %f" %(1.0 * m / (n * (n - 1) / 2))
    print "average degree = %f" %(numpy.average(g.degree().values()))
    components = networkx.connected_components(g)
    print "number of components = %d" %(len(components))
    print "fractional size of largest component = %f" \
        %(max([len(i) * 1. for i in components]) / \
              len(g.nodes()))
    if len(components) == 1:
        print "diameter = %f" %(networkx.average_shortest_path_length(g)) 
    else:
        print "diameter = infinity"
    print "average clustering coefficient = %f" \
        %(networkx.average_clustering(g))
    print "homophily = %f" %(networkx.degree_assortativity_coefficient(g))

if __name__ == "__main__":
    main(sys.argv[1:])
