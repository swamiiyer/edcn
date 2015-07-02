# This script defines an abstraction for a player in the game, modules 
# for carrying out the dynamics of the game, and several helper functions.

import math, network, random

class Player(network.Vertex):
    """
    Abstraction for an individual in the network. Each individual has an id, 
    a trait, a payoff, and maintains a  history of the trait and payoff values.
    In addition, each individual also has a fitness value that is used instead 
    of payoff in case of certain kinds (imitation and birth-death) of dynamics, 
    and an inherited trait -- a trait that is inherited from another 
    individual during the update round of a game -- which becomes the current 
    trait once the update round for a generation is complete.
    """

    def __init__(self, id):
        """
        Construct a Player object given the id. 
        """

        self.id = id
        self.trait_current = 0.0
        self.payoff_current = 0.0
        self.trait = []
        self.payoff = []
        self.fitness_current = 0.0
        self.trait_inherited = 0.0

    def set_trait(self, v):
        """
        Set trait of this invdividual to the specified value.
        """

        self.trait_current = v

    def get_trait(self):
        """
        Return trait of this individual.
        """

        return self.trait_current

    def get_trait_list(self):
        """
        Return trait history of this individual.
        """

        return self.trait

    def inherit_trait(self, v):
        """
        Set v as the trait to be inherited by this individual.
        """
        
        self.trait_inherited = v

    def commit_inheritance(self):
        """
        Set the trait of this individual to what was inherited from another 
        individual.
        """

        self.trait_current = self.trait_inherited

    def set_payoff(self, v):
        """
        Set payoff of this individual to the specified value.
        """

        self.payoff_current = v

    def get_payoff(self):
        """
        Return payoff of this individual.
        """

        return self.payoff_current

    def get_payoff_list(self):
        """
        Return payoff history of this individual.
        """

        return self.payoff

    def get_fitness(self):
        """
        Return fitness of this individual.
        """
        
        return self.fitness_current

    def set_fitness(self, v):
        """
        Set fitness of this individual to the specified value.
        """
        
        self.fitness_current = v
    
    def save(self):
        """
        Record the current trait and payoff values in the respective histories.
        """

        self.trait.append(self.trait_current)
        self.payoff.append(self.payoff_current)

    def __str__(self):
        """
        Return a string representation for this individual.
        """

        return "Player (id = %d, trait = %f, payoff = %f)" \
            %(self.id, self.trait_current, self.payoff_current)

class Dynamics:
    """
    Abstraction for a module for carrying out the game dynamics. 
    """

    def __init__(self, net, params, payoff):
        """
        Construct a game dynamics module given the network structure which 
        is an object of type netpy.Network, the game parameters, and the 
        payoff function.
        """

        self.net = net
        self.params = params
        self.payoff = payoff
        self.total_fitness = 0.0
        self.fitness_cdf = []
        self.min_payoff = float("inf")
        self.max_payoff = float("-inf")

    def pre_interaction(self):
        """
        Invoked immediately before the round of interactions. Reinitialize 
        the min and max payoff of the population.
        """

        self.min_payoff = float("inf")
        self.max_payoff = float("-inf")

    def interact(self):
        """
        Carry out an interaction.
        """

        pass

    def post_interaction(self):
        """
        Computations that need to be performed immediately after the round of 
        interactions must be carried out here.
        """

        pass

    def update(self):
        """
        Carry out an update.
        """

        pass

    def post_update(self):
        """
        Called after the round of updates. Commit inheritance of trait for 
        each individual in the population.
        """

        for vertex in self.net.get_vertices():
            vertex.commit_inheritance()

    def fitness(self, p):
        """
        Return the fitness given the payoff p.
        """
        
        # f = 1.0 - self.params["selection_strength"] + \
        #     self.params["selection_strength"] * p
        # assert f >= 0.0
        f = math.exp(self.params["selection_strength"] * p)
        return f

    def replicate(self, pi, pj):
        """ 
        Returns the probability that i will inherit j's trait given their 
        payoffs pi and pj.
        """
        
        if pi >= pj:
            return 0.0
        return self.params["selection_strength"] * \
            (pj - pi) / (self.max_payoff - self.min_payoff)

    def fermi(self, pi, pj):
        """ 
        Returns the probability that i will inherit j's trait. pi and 
        pj are the payoffs of i and j.
        """

        return 1.0 / (1.0 + math.e ** (-self.params["selection_strength"] * \
                                            (pj - pi))) 

    def roulette(self, l, total_fitness):
        """
        Given the cumulative distribution l of fitness values and the total 
        fitness of the population, using binary search pick and return an 
        individual from l in proportion to fitness.
        """

        r = random.random() * total_fitness
        lo, hi, mid = 0, len(l) - 1, -1
        while lo <= hi:
            mid = lo + (hi - lo) / 2
            if r > l[mid][0]:
                lo = mid + 1
            elif r < l[mid][0]:
                hi = mid - 1
        return l[mid][1]

class Unstructured_Replication1(Dynamics):
    """
    Replication dynamics on an unstructured network.
    """

    def interact(self):
        """
        Carry out an interaction.
        """

        i = self.net.get_random_vertex()
        if random.random() < self.params["assortativity"]:
            payoff = self.payoff(i.get_trait(), i.get_trait())
        else:
            group_trait = 0.0
            if self.params["group"] == 0:
                for j in self.net.get_neighbors(i):
                    group_trait += j.get_trait()
                group_trait -= i.get_trait()
                # Comment the following line if the sum of the traits of 
                # neighbors is desired instead of average.
                group_trait /= (self.params["population"] - 1)
            else:
                for count in range(0, self.params["group"]):
                    j = self.net.get_random_neighbor(i)
                    group_trait += j.get_trait()
                # Comment the following line if the sum of the traits of 
                # neighbors is desired instead of average.
                group_trait /= self.params["group"]
            payoff = self.payoff(i.get_trait(), group_trait)
        i.set_payoff(payoff)

    def post_interaction(self):
        """
        Calculate the min and max payoff of the population.
        """
        
        for vertex in self.net.get_vertices():
            if vertex.get_payoff() < self.min_payoff:
                self.min_payoff = vertex.get_payoff()
            if vertex.get_payoff() > self.max_payoff:
                self.max_payoff = vertex.get_payoff()

    def update(self):
        """
        Carry out an update.
        """

        i = self.net.get_random_vertex()
        j = self.net.get_random_neighbor(i)
        p = self.replicate(i.get_payoff(), j.get_payoff()) 
        trait = i.get_trait()
        if random.random() < p: 
            trait = j.get_trait()
        if random.random() < self.params["mutation"]: 
            trait = max(0.0, min(random.gauss(trait, 
                                              self.params["stddev"]), 
                                 self.params["max_trait"]))
        i.inherit_trait(trait)

class Unstructured_Replication2(Dynamics):
    """
    Replication dynamics on an unstructured network.
    """

    def interact(self):
        """
        Carry out an interaction followed by an update.
        """

        # Interact.
        i = self.net.get_random_vertex()
        j = self.net.get_random_neighbor(i)
        if random.random() < self.params["assortativity"]:
            pi = self.payoff(i.get_trait(), i.get_trait())
        else:
            k = self.net.get_random_neighbor(i)
            if k == None:
                return
            pi = self.payoff(i.get_trait(), k.get_trait())
        if random.random() < self.params["assortativity"]:
            pj = self.payoff(j.get_trait(), j.get_trait())
        else:
            l = self.net.get_random_neighbor(j)
            if l == None:
                return
            pj = self.payoff(j.get_trait(), l.get_trait())
        i.set_payoff(pi)
        j.set_payoff(pj)
        self.min_payoff = min(self.min_payoff, pi, pj)
        self.max_payoff = max(self.max_payoff, pi, pj)

        # Update.
        p = self.replicate(pi, pj)
        trait = i.get_trait()
        if random.random() < p: 
            trait = j.get_trait()
        if random.random() < self.params["mutation"]: 
            trait = max(0.0, min(random.gauss(trait, 
                                              self.params["stddev"]), 
                                 self.params["max_trait"]))
        i.inherit_trait(trait)

    def post_interaction(self):
        """
        Nothing to do here.
        """

        pass

    def update(self):
        """
        Nothing to do here.
        """

        pass

class Unstructured_Fermi1(Dynamics):
    """
    Fermi dynamics on an unstructured network.
    """

    def interact(self):
        """
        Carry out an interaction.
        """

        i = self.net.get_random_vertex()
        if random.random() < self.params["assortativity"]:
            payoff = self.payoff(i.get_trait(), i.get_trait())
        else:
            group_trait = 0.0
            if self.params["group"] == 0:
                for j in self.net.get_neighbors(i):
                    group_trait += j.get_trait()
                group_trait -= i.get_trait()
                # Comment the following line if the sum of the traits of 
                # neighbors is desired instead of average.
                group_trait /= (self.params["population"] - 1)
            else:
                for count in range(0, self.params["group"]):
                    j = self.net.get_random_neighbor(i)
                    group_trait += j.get_trait()
                # Comment the following line if the sum of the traits of 
                # neighbors is desired instead of average.
                group_trait /= self.params["group"]
            payoff = self.payoff(i.get_trait(), group_trait)
        i.set_payoff(payoff)

    def post_interaction(self):
        """
        Nothing to do here.
        """

        pass

    def update(self):
        """
        Carry out an update.
        """

        i = self.net.get_random_vertex()
        j = self.net.get_random_neighbor(i)
        p = self.fermi(i.get_payoff(), j.get_payoff())
        if random.random() < p: 
            trait = j.get_trait()
            if random.random() < self.params["mutation"]: 
                trait = max(0.0, min(random.gauss(trait, 
                                                  self.params["stddev"]), 
                                     self.params["max_trait"]))
            i.inherit_trait(trait)

class Unstructured_Fermi2(Dynamics):
    """
    Fermi dynamics on an unstructured network.
    """

    def interact(self):
        """
        Carry out an interaction followed by an update.
        """

        # Interact.
        i = self.net.get_random_vertex()
        j = self.net.get_random_neighbor(i)
        if random.random() < self.params["assortativity"]:
            pi = self.payoff(i.get_trait(), i.get_trait())
        else:
            k = self.net.get_random_neighbor(i)
            pi = self.payoff(i.get_trait(), k.get_trait())
        if random.random() < self.params["assortativity"]:
            pj = self.payoff(j.get_trait(), j.get_trait())
        else:
            l = self.net.get_random_neighbor(j)
            pj = self.payoff(j.get_trait(), l.get_trait())
        i.set_payoff(pi)
        j.set_payoff(pj)

        # Update.
        p = self.fermi(pi, pj)
        if random.random() < p: 
            trait = j.get_trait()
            if random.random() < self.params["mutation"]: 
                trait = max(0.0, min(random.gauss(trait, 
                                                  self.params["stddev"]), 
                                     self.params["max_trait"]))
            i.inherit_trait(trait)

    def post_interaction(self):
        """
        Nothing to do here.
        """

        pass

    def update(self):
        """
        Nothing to do here.
        """

        pass

class Unstructured_Imitation(Dynamics):
    """
    Imitation dynamics on an unstructured network.
    """

    def interact(self):
        """
        Carry out an interaction.
        """

        i = self.net.get_random_vertex()
        if random.random() < self.params["assortativity"]:
            payoff = self.payoff(i.get_trait(), i.get_trait())
        else:
            group_trait = 0.0
            if self.params["group"] == 0:
                for j in self.net.get_neighbors(i):
                    group_trait += j.get_trait()
                group_trait -= i.get_trait()
                # Comment the following line if the sum of the traits of 
                # neighbors is desired instead of average.
                group_trait /= (self.params["population"] - 1)
            else:
                for count in range(0, self.params["group"]):
                    j = self.net.get_random_neighbor(i)
                    group_trait += j.get_trait()
                # Comment the following line if the sum of the traits of 
                # neighbors is desired instead of average.
                group_trait /= self.params["group"]
            payoff = self.payoff(i.get_trait(), group_trait)
        i.set_payoff(payoff)
        i.set_fitness(self.fitness(payoff))

    def post_interaction(self):
        """
        Calculate the fitness CDF and the total fitness of the population.
        """

        self.total_fitness = 0.0
        self.fitness_cdf = []
        for vertex in self.net.get_vertices():
            self.total_fitness += vertex.get_fitness()
            self.fitness_cdf.append([self.total_fitness, vertex])

    def update(self):
        """
        Carry out an update.
        """

        i = self.net.get_random_vertex()
        j = self.roulette(self.fitness_cdf, self.total_fitness)
        trait = j.get_trait()
        if random.random() < self.params["mutation"]: 
            trait = max(0.0, min(random.gauss(trait, self.params["stddev"]), 
                                 self.params["max_trait"]))
            i.inherit_trait(trait)

class Unstructured_Birth_Death(Dynamics):
    """
    Birth-death dynamics on an unstructured network.
    """

    def interact(self):
        """
        Carry out an interaction.
        """

        i = self.net.get_random_vertex()
        if random.random() < self.params["assortativity"]:
            payoff = self.payoff(i.get_trait(), i.get_trait())
        else:
            group_trait = 0.0
            if self.params["group"] == 0:
                for j in self.net.get_neighbors(i):
                    group_trait += j.get_trait()
                group_trait -= i.get_trait()
                # Comment the following line if the sum of the traits of 
                # neighbors is desired instead of average.
                group_trait /= (self.params["population"] - 1)
            else:
                for count in range(0, self.params["group"]):
                    j = self.net.get_random_neighbor(i)
                    group_trait += j.get_trait()
                # Comment the following line if the sum of the traits of 
                # neighbors is desired instead of average.
                group_trait /= self.params["group"]
            payoff = self.payoff(i.get_trait(), group_trait)
        i.set_payoff(payoff)
        i.set_fitness(self.fitness(payoff))

    def post_interaction(self):
        """
        Calculate the fitness CDF and total fitness of the population.
        """

        self.total_fitness = 0.0
        self.fitness_cdf = []
        for vertex in self.net.get_vertices():
            self.total_fitness += vertex.get_fitness()
            self.fitness_cdf.append([self.total_fitness, vertex])

    def update(self):
        """
        Carry out an update.
        """

        j = self.roulette(self.fitness_cdf, self.total_fitness)
        i = self.net.get_random_neighbor(j)
        trait = j.get_trait()
        if random.random() < self.params["mutation"]: 
            trait = max(0.0, min(random.gauss(trait, self.params["stddev"]), 
                                 self.params["max_trait"]))
            i.inherit_trait(trait)

class Unstructured_Death_Birth(Dynamics):
    """
    Death-birth dynamics on an unstructured network.
    """

    def interact(self):
        """
        Carry out an interaction.
        """

        i = self.net.get_random_vertex()
        if random.random() < self.params["assortativity"]:
            payoff = self.payoff(i.get_trait(), i.get_trait())
        else:
            group_trait = 0.0
            if self.params["group"] == 0:
                for j in self.net.get_neighbors(i):
                    group_trait += j.get_trait()
                group_trait -= i.get_trait()
                # Comment the following line if the sum of the traits of 
                # neighbors is desired instead of average.
                group_trait /= (self.params["population"] - 1)
            else:
                for count in range(0, self.params["group"]):
                    j = self.net.get_random_neighbor(i)
                    group_trait += j.get_trait()
                # Comment the following line if the sum of the traits of 
                # neighbors is desired instead of average.
                group_trait /= self.params["group"]
            payoff = self.payoff(i.get_trait(), group_trait)
        i.set_payoff(payoff)
        i.set_fitness(self.fitness(payoff))

    def post_interaction(self):
        """
        Calculate the fitness CDF and total fitness of the population.
        """

        self.total_fitness = 0.0
        self.fitness_cdf = []
        for vertex in self.net.get_vertices():
            self.total_fitness += vertex.get_fitness()
            self.fitness_cdf.append([self.total_fitness, vertex])

    def update(self):
        """
        Carry out an update.
        """

        i = self.net.get_random_vertex()
        j = self.roulette(fitness_cdf, total_fitness)
        trait = j.get_trait()
        if random.random() < self.params["mutation"]: 
            trait = max(0.0, min(random.gauss(trait, self.params["stddev"]), 
                                 self.params["max_trait"]))
            i.inherit_trait(trait)

class Structured_Replication1(Dynamics):
    """
    Replication dynamics on structured network.
    """

    def interact(self):
        """
        Carry out an interaction.
        """

        i = self.net.get_random_vertex()
        group_trait = 0.0
        if self.params["group"] == 0:
            for j in self.net.get_neighbors(i):
                group_trait += j.get_trait()
            # Comment the following line if the sum of the traits of 
            # neighbors is desired instead of average.
            group_trait /= len(self.net.get_neighbors(i))
        else:
            for count in range(0, self.params["group"]):
                j = self.net.get_random_neighbor(i)
                if j == None:
                    return
                group_trait += j.get_trait()
            # Comment the following line if the sum of the traits of 
            # neighbors is desired instead of average.
            group_trait /= self.params["group"]
        payoff = self.payoff(i.get_trait(), group_trait)
        i.set_payoff(payoff)

    def post_interaction(self):
        """
        Calculate the min and max payoff of the population.
        """
        
        for vertex in self.net.get_vertices():
            if vertex.get_payoff() < self.min_payoff:
                self.min_payoff = vertex.get_payoff()
            if vertex.get_payoff() > self.max_payoff:
                self.max_payoff = vertex.get_payoff()

    def update(self):
        """
        Carry out an update.
        """

        i = self.net.get_random_vertex()
        j = self.net.get_random_neighbor(i)
        if j == None:
            return
        p = self.replicate(i.get_payoff(), j.get_payoff())
        trait = i.get_trait()
        if random.random() < p: 
            trait = j.get_trait()
        if random.random() < self.params["mutation"]: 
            trait = max(0.0, min(random.gauss(trait, 
                                              self.params["stddev"]), 
                                 self.params["max_trait"]))
        i.inherit_trait(trait)

class Structured_Replication2(Dynamics):
    """
    Replication dynamics on a structured network.
    """

    def interact(self):
        """
        Carry out an interaction followed by an update.
        """

        # Interact.
        i = self.net.get_random_vertex()
        j = self.net.get_random_neighbor(i)
        if j == None:
            return
        k = self.net.get_random_neighbor(i)
        l = self.net.get_random_neighbor(j)
        if k == None or l == None or i == l or j == k:
            return
        pi = self.payoff(i.get_trait(), k.get_trait())
        pj = self.payoff(j.get_trait(), l.get_trait())
        i.set_payoff(pi)
        j.set_payoff(pj)
        self.min_payoff = min(self.min_payoff, pi, pj)
        self.max_payoff = max(self.max_payoff, pi, pj)

        # Update.
        p = self.replicate(pi, pj)
        trait = i.get_trait()
        if random.random() < p: 
            trait = j.get_trait()
        if random.random() < self.params["mutation"]: 
            trait = max(0.0, min(random.gauss(trait, 
                                              self.params["stddev"]), 
                                 self.params["max_trait"]))
        i.inherit_trait(trait)

    def post_interaction(self):
        """
        Nothing to do here.
        """

        pass

    def update(self):
        """
        Nothing to do here.
        """

        pass
        
class Structured_Fermi1(Dynamics):
    """
    Fermi dynamics on a structured network.
    """

    def interact(self):
        """
        Carry out an interaction.
        """

        i = self.net.get_random_vertex()
        group_trait = 0.0
        if self.params["group"] == 0:
            for j in self.net.get_neighbors(i):
                group_trait += j.get_trait()
            # Comment the following line if the sum of the traits of 
            # neighbors is desired instead of average.
            group_trait /= len(self.net.get_neighbors(i))
        else:
            for count in range(0, self.params["group"]):
                j = self.net.get_random_neighbor(i)
                if j == None:
                    return
                group_trait += j.get_trait()
            # Comment the following line if the sum of the traits of 
            # neighbors is desired instead of average.
            group_trait /= self.params["group"]
        payoff = self.payoff(i.get_trait(), group_trait)
        i.set_payoff(payoff)

    def post_interaction(self):
        """
        Nothing to do here.
        """

        pass

    def update(self):
        """
        Carry out an update.
        """

        i = self.net.get_random_vertex()
        j = self.net.get_random_neighbor(i)
        if j == None:
            return
        p = self.fermi(i.get_payoff(), j.get_payoff()) 
        if random.random() < p: 
            trait = j.get_trait()
            if random.random() < self.params["mutation"]: 
                trait = max(0.0, min(random.gauss(trait, 
                                                  self.params["stddev"]), 
                                     self.params["max_trait"]))
            i.inherit_trait(trait)

class Structured_Fermi2(Dynamics):
    """
    Fermi dynamics on a structured network.
    """

    def interact(self):
        """
        Carry out an interaction followed by an update.
        """

        # Interact.
        i = self.net.get_random_vertex()
        j = self.net.get_random_neighbor(i)
        if j == None:
            return
        k = self.net.get_random_neighbor(i)
        l = self.net.get_random_neighbor(j)
        if k == None or l == None or i == l or j == k:
            return
        pi = self.payoff(i.get_trait(), k.get_trait())
        pj = self.payoff(j.get_trait(), l.get_trait())
        i.set_payoff(pi)
        j.set_payoff(pj)

        # Update.
        p = self.fermi(pi, pj)
        if random.random() < p: 
            trait = j.get_trait()
            if random.random() < self.params["mutation"]: 
                trait = max(0.0, min(random.gauss(trait, 
                                                  self.params["stddev"]), 
                                     self.params["max_trait"]))
            i.inherit_trait(trait)

    def post_interaction(self):
        """
        Nothing to do here.
        """

        pass

    def update(self):
        """
        Nothing to do here.
        """

        pass

class Structured_Imitation(Dynamics):
    """
    Imitation dynamics on a structured network.
    """

    def interact(self):
        """
        Carry out an interaction.
        """

        i = self.net.get_random_vertex()
        group_trait = 0.0
        if self.params["group"] == 0:
            for j in self.net.get_neighbors(i):
                group_trait += j.get_trait()
            # Comment the following line if the sum of the traits of 
            # neighbors is desired instead of average.
            group_trait /= len(self.net.get_neighbors(i))
        else:
            for count in range(0, self.params["group"]):
                j = self.net.get_random_neighbor(i)
                if j == None:
                    return
                group_trait += j.get_trait()
            # Comment the following line if the sum of the traits of 
            # neighbors is desired instead of average.
            group_trait /= self.params["group"]
        payoff = self.payoff(i.get_trait(), group_trait)
        i.set_payoff(payoff)
        i.set_fitness(self.fitness(payoff))

    def post_interaction(self):
        """
        Nothing to do here.
        """

        pass

    def update(self):
        """
        Carry out an update.
        """

        i = self.net.get_random_vertex()
        total_fitness = 0.0
        fitness_cdf = []
        for vertex in self.net.get_neighbors(i):
            total_fitness += vertex.get_fitness()
            fitness_cdf.append([self.total_fitness, vertex])
        total_fitness += i.get_fitness()
        fitness_cdf.append([total_fitness, i])
        j = self.roulette(fitness_cdf, total_fitness)
        trait = j.get_trait()
        if random.random() < self.params["mutation"]: 
            trait = max(0.0, min(random.gauss(trait, self.params["stddev"]), 
                                 self.params["max_trait"]))
        i.inherit_trait(trait)

class Structured_Birth_Death(Dynamics):
    """
    Birth-death dynamics on a structured network.
    """

    def interact(self):
        """
        Carry out an interaction.
        """

        i = self.net.get_random_vertex()
        group_trait = 0.0
        if self.params["group"] == 0:
            for j in self.net.get_neighbors(i):
                group_trait += j.get_trait()
            # Comment the following line if the sum of the traits of 
            # neighbors is desired instead of average.
            group_trait /= len(self.net.get_neighbors(i))
        else:
            for count in range(0, self.params["group"]):
                j = self.net.get_random_neighbor(i)
                if j == None:
                    return
                group_trait += j.get_trait()
            # Comment the following line if the sum of the traits of 
            # neighbors is desired instead of average.
            group_trait /= self.params["group"]
        payoff = self.payoff(i.get_trait(), group_trait)
        i.set_payoff(payoff)
        i.set_fitness(self.fitness(payoff))

    def post_interaction(self):
        """
        Calculate the fitness CDF and the total fitness of the population.
        """
        
        self.total_fitness = 0.0
        self.fitness_cdf = []
        for vertex in self.net.get_vertices():
            self.total_fitness += vertex.get_fitness()
            self.fitness_cdf.append([self.total_fitness, vertex])

    def update(self):
        """
        Carry out an update.
        """

        j = self.roulette(self.fitness_cdf, self.total_fitness)
        i = self.net.get_random_neighbor(j)
        if i == None:
            return
        trait = j.get_trait()
        if random.random() < self.params["mutation"]: 
            trait = max(0.0, min(random.gauss(trait, self.params["stddev"]), 
                                 self.params["max_trait"]))
        i.inherit_trait(trait)

class Structured_Death_Birth(Dynamics):
    """
    Death-birth dynamics on a structured network.
    """

    def interact(self):
        """
        Carry out an interaction.
        """

        i = self.net.get_random_vertex()
        group_trait = 0.0
        if self.params["group"] == 0:
            for j in self.net.get_neighbors(i):
                group_trait += j.get_trait()
            # Comment the following line if the sum of the traits of 
            # neighbors is desired instead of average.
            group_trait /= len(self.net.get_neighbors(i))
        else:
            for count in range(0, self.params["group"]):
                j = self.net.get_random_neighbor(i)
                if j == None:
                    return
                group_trait += j.get_trait()
            # Comment the following line if the sum of the traits of 
            # neighbors is desired instead of average.
            group_trait /= self.params["group"]
        payoff = self.payoff(i.get_trait(), group_trait)
        i.set_payoff(payoff)
        i.set_fitness(self.fitness(payoff))

    def post_interaction(self):
        """
        Nothing to do here.
        """

        pass

    def update(self):
        """
        Carry out an update.
        """

        i = self.net.get_random_vertex()
        total_fitness = 0.0
        fitness_cdf = []
        for vertex in self.net.get_neighbors(i):
            total_fitness += vertex.get_fitness()
            fitness_cdf.append([total_fitness, vertex])
        j = self.roulette(fitness_cdf, total_fitness)
        trait = j.get_trait()
        if random.random() < self.params["mutation"]: 
            trait = max(0.0, min(random.gauss(trait, self.params["stddev"]), 
                                 self.params["max_trait"]))
        i.inherit_trait(trait)

def dynamics_module(topology, update_rule):
    """
    Return the appropriate game dynamics module given the network 
    topology and the update rule.
    """

    dm = None
    if topology == "Complete":
        if update_rule == "RE1":
            dm = lambda a, b, c: Unstructured_Replication1(a, b, c)
        elif update_rule == "RE2":
            dm = lambda a, b, c: Unstructured_Replication2(a, b, c)
        elif update_rule == "FE1":
            dm = lambda a, b, c: Unstructured_Fermi1(a, b, c)
        elif update_rule == "FE2":
            dm = lambda a, b, c: Unstructured_Fermi2(a, b, c)
        elif update_rule == "IM":
            dm = lambda a, b, c: Unstructured_Imitation(a, b, c)
        elif update_rule == "BD":
            dm = lambda a, b, c: Unstructured_Birth_Death(a, b, c)
        elif update_rule == "DB":
            dm = lambda a, b, c: Unstructured_Death_Birth(a, b, c)
    else:
        if update_rule == "RE1":
            dm = lambda a, b, c: Structured_Replication1(a, b, c)
        elif update_rule == "RE2":
            dm = lambda a, b, c: Structured_Replication2(a, b, c)
        elif update_rule == "FE1":
            dm = lambda a, b, c: Structured_Fermi1(a, b, c)
        elif update_rule == "FE2":
            dm = lambda a, b, c: Structured_Fermi2(a, b, c)
        elif update_rule == "IM":
            dm = lambda a, b, c: Structured_Imitation(a, b, c)
        elif update_rule == "BD":
            dm = lambda a, b, c: Structured_Birth_Death(a, b, c)
        elif update_rule == "DB":
            dm = lambda a, b, c: Structured_Death_Birth(a, b, c)
    return dm

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
