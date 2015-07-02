# This script defines an abstraction for a player in the game, and several 
# helper functions.

import math, network, networkx, operator, os, pickle, random

class Player(network.Vertex):
    """
    Abstraction for an individual in the network. Each individual has an id, 
    a infection (infection probability), c (disease-specific death 
    probability), and state (susceptible, infected, recovered, vaccinated) 
    value, and also maintains a history of these values.
    """ 

    # Any player is in one of these states.
    SUSCEPTIBLE = 0
    INFECTED    = 1
    RECOVERED   = 2
    VACCINATED  = 3

    def __init__(self, id):
        """
        Construct a Player given the id. 
        """

        self.id = id
        self.infection = []
        self.c = []
        self.state = []
        self.set_infection(0.0)
        self.set_c(0.0)
        self.set_state(Player.SUSCEPTIBLE)

    def set_infection(self, v):
        """
        Set the transmission rate for this individual. 
        """

        self.infection_current = v
        
    def get_infection(self):
        """
        Return the transmission rate for this individual.
        """

        return self.infection_current

    def get_infection_list(self):
        """
        Return the transmission rate history of this individual.
        """

        return self.infection

    def set_c(self, v):
        """
        Set the disease-induced death rate to the specified value.
        """
        
        self.c_current = v

    def get_c(self):
        """
        Return the disease-induced death rate for this individual.
        """

        return self.c_current
    
    def get_c_list(self):
        """
        Return the disease-induced death rate history of this individual.
        """

        return self.c

    def set_state(self, v):
        """
        Set the state of this individual.
        """

        self.state_current = v

    def get_state(self):
        """
        Return the state of this individual.
        """

        return self.state_current

    def get_state_list(self):
        """
        Return the state history of this individual.
        """

        return self.state

    def save(self):
        """
        Record the transmission rate, disease-induced death rate, and the state 
        of this individual in the respective histories.
        """

        self.infection.append(self.infection_current)
        self.c.append(self.c_current)
        self.state.append(self.state_current)

    def __str__(self):
        """
        Return a string representation for this individual.
        """

        return "Player (id = %d, infection = %f, c = %f, state = %d)" \
            %(self.id, self.infection_current, self.c_current, 
              self.state_current)

def state_stats(infile):
    """
    Lists the number of individuals in different states in the results file 
    specified by infile.
    """

    fh = open(infile, "r")
    params = pickle.load(fh)
    population = pickle.load(fh)
    fh.close()
    s = 0
    i = 0
    r = 0
    v = 0
    for p in population:
        if p.get_state() == Player.SUSCEPTIBLE:
            s += 1
        elif p.get_state() == Player.INFECTED:
            i += 1
        elif p.get_state() == Player.RECOVERED:
            r += 1
        elif p.get_state() == Player.VACCINATED:
            v += 1
    print "N = %d, S = %d, I = %d, R = %d, V = %d" %(len(population), 
                                                     s, i, r, v)

def degree_vaccination(infile, outfile, k):
    """
    Change the state (to VACCINATED) of the top k (%) highly-connected 
    individuals in population contained in infile, and save the altered 
    population into outfile. 
    """

    fh = open(infile, "r")
    params = pickle.load(fh)
    population = pickle.load(fh)
    fh.close()
    net = network.build_network(population, params["network_topology"], 
                                params["network_params"])
    L = networkx.degree_centrality(net.get_structure()).items()
    L = sorted(L, key = operator.itemgetter(1), reverse = True)
    n = int(k / 100.0 * len(population))
    for i in range(0, n):
        person = population[L[i][0]]
        person.set_state(Player.VACCINATED)
        sl = person.get_state_list()
        sl[len(sl) - 1] = Player.VACCINATED
        person.set_infection(0.0)
        il = person.get_infection_list()
        il[len(il) - 1] = 0.0
        person.set_c(0.0)
        cl = person.get_c_list()
        cl[len(cl) - 1] = 0.0
    fh = open(outfile, "w")
    pickle.dump(params, fh)
    pickle.dump(population, fh)
    fh.close()

def referral_vaccination(infile, outfile, k):
    """
    Change the state (to VACCINATED) of k (%) individuals in population 
    contained in infile, and save the altered population into outfile. Each of 
    these individuals is a randomly picked neighbor of a randomly picked 
    person. 
    """

    fh = open(infile, "r")
    params = pickle.load(fh)
    population = pickle.load(fh)
    fh.close() 
    net = network.build_network(population, params["network_topology"], 
                                params["network_params"])
    m = int(k / 100.0 * len(population))
    i = 0
    while i < m:
        person = net.get_random_vertex()
        person = net.get_random_neighbor(person)
        if person == None or person.get_state() == Player.VACCINATED: 
            continue
        person.set_state(Player.VACCINATED)
        sl = person.get_state_list()
        sl[len(sl) - 1] = Player.VACCINATED
        person.set_infection(0.0)
        il = person.get_infection_list()
        il[len(il) - 1] = 0.0
        person.set_c(0.0)
        cl = person.get_c_list()
        cl[len(cl) - 1] = 0.0
        i += 1
    fh = open(outfile, "w")
    pickle.dump(params, fh)
    pickle.dump(population, fh)
    fh.close()

def random_vaccination(infile, outfile, k):
    """
    Change the state (to VACCINATED) of k (%) individuals in population 
    contained in infile, and save the altered population into outfile. Each of 
    these individuals is randomly picked. 
    """

    fh = open(infile, "r")
    params = pickle.load(fh)
    population = pickle.load(fh)
    fh.close()
    net = network.build_network(population, params["network_topology"], 
                                params["network_params"])
    m = int(k / 100.0 * len(population))
    i = 0
    while i < m:
        person = net.get_random_vertex()
        if person == None or person.get_state() == Player.VACCINATED:
            continue
        person.set_state(Player.VACCINATED)
        sl = person.get_state_list()
        sl[len(sl) - 1] = Player.VACCINATED
        person.set_infection(0.0)
        il = person.get_infection_list()
        il[len(il) - 1] = 0.0
        person.set_c(0.0)
        cl = person.get_c_list()
        cl[len(cl) - 1] = 0.0
        i += 1
    fh = open(outfile, "w")
    pickle.dump(params, fh)
    pickle.dump(population, fh)
    fh.close()

def fitness(infection, b, a):
    """
    Return the virulence given the infection probability, and constants 
    b and a.
    """
    
    return math.pow(infection / b, 1.0 / a)

def bins(items, low, high, n_bins):
    """
    Bin the elements in the items list into n_bins and return the list of 
    bins. The low and high values specifiy the minimum and maximum values for 
    the elements.
    """

    l = []
    delta = (high - low) / n_bins
    a = low
    b = a + delta
    for i in range(0, n_bins):
        n = len(filter(lambda x: x >= a and x < b, items))
        l.append(n)
        a = b
        b += delta
    return l
