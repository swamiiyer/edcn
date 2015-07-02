# This script defines a higher-level abstraction for vertices in a network and 
# different kinds networks.

import networkx, operator, random, time

class Vertex(object):
    """
    Abstraction of a vertex in a network. A richer representation for a vertex 
    capturing project-specific details can be obtained by inheriting 
    from this class.
    """

    def __init__(self, id):
        """
        Construct a Vertex given an id for it. The id must be an integer 
        in the range [0 - n] where n is the number of vertices in the network.
        """

        self.id = id

    def get_id(self):
        """
        Return the id of this vertex.
        """

        return self.id

    def __str__(self):
        """
        Return a string representation for this vertex.
        """
        return "Vertex (id = %d)" %(self.id)

class Network(object):
    """
    Abstraction for a network. Concrete networks (Barabasi Albert network 
    for example) inherit from this class. 
    """

    def get_vertices(self):
        """
        Return the list of vertices in the network.
        """

        return self.vertices

    def get_structure(self):
        """
        Return networkx graph object representing this network's structure.
        """

        return self.g
    
    def get_neighbors(self, vertex):
        """
        Return a list of neighbors of the specified vertex, or None.
        """

        l = None
        if self.__adj_list_map__.has_key(vertex.id):
            l = self.__adj_list_map__[vertex.id]
        else:
            l = [self.vertices[i] for i in self.g.neighbors(vertex.id)]
            self.__adj_list_map__[vertex.id] = l
        return l

    def get_random_neighbor(self, vertex):
        """
        Return a random neighbor of the specified vertex in the network, or 
        None.
        """

        if len(self.get_neighbors(vertex)) != 0:
            return random.choice(self.get_neighbors(vertex))
        return None

    def get_random_vertex(self):
        """
        Return a random vertex in the network.
        """

        return random.choice(self.vertices)

class Complete_Network(Network):
    """
    Representation of a Complete network, where every vertex is connected to 
    every other vertex.
    """
    
    def __init__(self, vertices):
        """
        Construct a Complete network, where every vertex is connected to 
        every other vertex. The vertices argument contains objects of 
        type Vertex representing the vertices of the network. 
        Note: Since we know that every vertex is connected to every other 
        vertex, for efficiency, we don't build a NetworkX graph.
        """

        self.vertices = vertices
        self.g = None
        self.__adj_list_map__ = {}

    def get_neighbors(self, vertex):
        """
        Return a list of neighbors of the specified vertex. Since every vertex 
        is connected to every other vertex, we return the whole population. 
        Note that we are assuming that a vertex is also connected to itself.
        """

        return self.vertices

    def __str__(self):
        """
        Return a string representation for this network.
        """

        return "Complete Network (n = %d)" %(len(self.vertices))

class Erdos_Renyi_Network(Network):
    """
    Representation of an Erdos Renyi network.
    """
    
    def __init__(self, vertices, p, seed = int(time.time())):
        """
        Construct an Erdos Renyi network. The vertices argument contains 
        objects of type Vertex representing the vertices of the network. The 
        p argument specifies that two randomly picked vertices have an edge 
        between them. The seed argument is the seed for the random number 
        generator used to build the network.
        """

        self.vertices = vertices
        self.p = p
        self.seed = seed
        self.g = networkx.erdos_renyi_graph(len(vertices), p, seed = seed)
        self.__adj_list_map__ = {}

    def __str__(self):
        """
        Return a string representation for this network.
        """
        
        return "Erdos Renyi Network (n = %d, p = %f, seed = %d)" \
            %(len(self.vertices), self.p, self.seed)

class Barabasi_Albert_Network(Network):
    """
    Representation of a Barabasi Albert network.
    """
    
    def __init__(self, vertices, m, seed = int(time.time())):
        """
        Construct a Barabasi Albert network. The vertices argument contains 
        objects of type Vertex representing the vertices of the network. The 
        m argument specifies the number of attachments made from the 
        new vertex to the existing vertices during the network construction 
        process. The seed argument is the seed for the random number 
        generator used to build the network.
        """

        self.vertices = vertices
        self.m = m
        self.seed = seed
        self.g = networkx.barabasi_albert_graph(len(vertices), m, seed = seed)
        self.__adj_list_map__ = {}

    def __str__(self):
        """
        Return a string representation for this network.
        """

        return "Barabasi Albert Network (n = %d, m = %d, seed = %d)" \
            %(len(self.vertices), self.m, self.seed)

class Powerlaw_Homophilic_Network(Network):
    """
    Representation of a Powerlaw Homophilic network.
    """
    
    def __init__(self, vertices, m, r, maxiter, seed = int(time.time())):
        """
        Construct an Powerlaw Homophilic network. The vertices argument 
        contains objects of type Vertex representing the vertices of the 
        network. The m argument specifies the number of attachments made from 
        the new vertex to the existing vertices during the network 
        construction process. The r value specifies the desired value 
        for homophily. The maxiter value specifies the maximum 
        number of iterations (rewiring attempts) made until the desired 
        homophily is reached. The seed argument is the seed 
        for the random number generator used to build the network. 

        The algorithm used to build the homophilic network is discussed 
        in the following article:

        Construction and properties of homophilic random networks by 
        R. Xulvi-Brunet, I.M. Sokolov. [arXiv:cond-mat/0405095v1]

        The algorithm terminates if the desired homophily is 
        attained or if the maximum number of iterations (maxiter) is reached, 
        which ever comes first. Raises an exception if maxiter is reached 
        before achieving the desired level of homophily. The desired 
        homophily is attained to an accuracy of three decimal places.
        """

        self.vertices = vertices
        self.m = m
        self.r = r
        self.maxiter = maxiter
        self.seed = seed
        self.g = networkx.barabasi_albert_graph(len(vertices), m, seed = seed)
        self.__adj_list_map__ = {}

        # Modify self.g to make it homophilic.
        rc = 0.0
        success = False
        for t in range(0, self.maxiter):
            e1, e2 = random.choice(self.g.edges()), \
                random.choice(self.g.edges())
            n1, n2, n3, n4 = e1[0], e1[1], e2[0], e2[1]
            if e1 == e2:
                continue
            l = [(n1, self.g.degree(n1)), (n2, self.g.degree(n2)), 
                 (n3, self.g.degree(n3)), (n4, self.g.degree(n4))]
            l = sorted(l, key = operator.itemgetter(1))
            n1, n2, n3, n4 = l[0][0], l[1][0], l[2][0], l[3][0]
            if self.r >= 0.0 and not self.g.has_edge(n1, n2) and \
                    not self.g.has_edge(n3, n4):
                self.g.remove_edges_from([(e1[0], e1[1]), (e2[0], e2[1])])
                self.g.add_edges_from([(n1, n2), (n3, n4)])
                rc = networkx.degree_assortativity_coefficient(self.g)
            elif self.r < 0.0 and not self.g.has_edge(n1, n4) and \
                    not self.g.has_edge(n2, n3):
                self.g.remove_edges_from([(e1[0], e1[1]), (e2[0], e2[1])])
                self.g.add_edges_from([(n1, n4), (n2, n3)])
                rc = networkx.degree_assortativity_coefficient(self.g)
            #print t, rc # for tracing purposes
            if (self.r > 0 and rc > self.r) or (self.r < 0 and rc < self.r):
                success = True
                break
        if not success:
            raise Exception("maxiter reached before achieving desired homophily")

    def __str__(self):
        """
        Return a string representation for this network.
        """

        return ("Powerlaw Homophilic Network (n = %d, m = %d, p = %f, " + \
            "r = %f, maxiter = %d, seed = %d)") \
            %(len(self.vertices), self.m, self.p, self.r, self.maxiter, self.seed)

class Powerlaw_Clustered_Network(Network):
    """
    Representation of a Powerlaw Clustered network.
    """
    
    def __init__(self, vertices, m, p, seed = int(time.time())):
        """
        Construct a Powerlaw Clustered network. The vertices argument contains 
        objects of type Vertex representing the vertices of the network. The 
        m argument specifies the number of attachments made from the 
        new vertex to the existing vertices during the network construction 
        process. The p argument specifies the probability of adding a triangle 
        after adding a random edge. The seed argument is the seed for 
        the random number generator used to build the network.
        """

        self.vertices = vertices
        self.m = m
        self.p = p
        self.seed = seed
        self.g = networkx.powerlaw_cluster_graph(len(vertices), m, p, seed = seed)
        self.__adj_list_map__ = {}

    def __str__(self):
        """
        Return a string representation for this network.
        """

        return \
            "Powerlaw Clustered Network (n = %d, m = %d, p = %f, seed = %d)" \
            %(len(self.vertices), self.m, self.p, self.seed)

class Random_Regular_Network(Network):
    """
    Representation of a random regular network.
    """

    def __init__(self, vertices, degree, seed = int(time.time())):
        """
        Construct a Random Regular Network. The vertices argument contains 
        objects of type Vertex representing the vertices of the network. The 
        degree argument specifies the degree of each vertex. The seed argument 
        is the seed for the random number generator used to build the network. 
        Note that an random regular network has no self loops or 
        parallel edges.
        """

        self.vertices = vertices
        self.degree = degree
        self.seed = seed
        self.g = networkx.random_regular_graph(degree, len(vertices), seed = seed)
        self.__adj_list_map__ = {}

    def __str__(self):
        """
        Return a string representation for this network.
        """

        return "Random Regular Network (n = %d, degree = %d)" \
            %(len(self.vertices), self.degree)

class Grid_2D_Network(Network):
    """
    Representation of a 2D grid network with 4-cell Neumann neighborhood.
    """

    def __init__(self, vertices, m, n, periodic = False):
        """
        Construct a 2D Grid Network. The vertices argument contains 
        objects of type Vertex representing the vertices of the network. The m 
        and n arguments indicate the number of rows and columns in the grid. 
        The periodic argument determines if the boundary vertices are 
        connected via periodic boundary conditions. 
        """

        self.vertices = vertices
        self.m = m
        self.n = n
        self.periodic = periodic
        g = networkx.grid_2d_graph(m, n, periodic)
        self.g = networkx.convert_node_labels_to_integers(g)
        self.__adj_list_map__ = {}

    def __str__(self):
        """
        Return a string representation for this network.
        """

        return "Grid 2D Network (n = %d, m = %d, n = %d, periodic = %s)" \
            %(len(self.vertices), self.m, self.n, self.periodic)

class KNN_Network(Network):
    """
    Representation of a KNN (K Nearest Neighbors) network.
    """

    def __init__(self, vertices, n, k, periodic = False):
        """
        Construct a KNN (K Nearest Neighbors) Network. The vertices argument 
        contains objects of type Vertex representing the vertices of the 
        network. The n argument specifies size (n^2) of the network. The k 
        argument specifies the neighborhood size of a vertex; for example 8 
        implies 8 a 3x3 neighborhood minus the vertex itself, 24 implies a 5x5 
        neighborhood minus the vertex itself, and so on. The periodic 
        argument determines if the boundary vertices are connected via 
        periodic boundary conditions. 
        """

        self.vertices = vertices
        self.n = n
        self.k = k
        self.periodic = periodic
        self.__adj_list_map__ = {}

        # Add n^2 vertices to an empty graph.
        self.g = networkx.Graph()
        for i in range(0, n ** 2):
            self.g.add_node(i)
        
        # Connect up each vertex to it nearest neighbors up to depth k.
        l = (int((k + 1) ** 0.5) - 1) / 2
        if self.periodic == True: 
            for r in range(0, self.n):
                for c in range(0, self.n):
                    from_id = r * self.n + c
                    for i in range(r - l, r + l + 1):
                        for j in range(c - l, c + l + 1):
                            to_id = i % self.n * self.n  + j % self.n
                            if from_id != to_id:
                                self.g.add_edge(from_id, to_id)
        else:
            for r in range(0, self.n):
                for c in range(0, self.n):
                    from_id = r * self.n + c
                    for i in range(r - l, r + l + 1):
                        if i < 0 or i >= self.n:
                            continue
                        for j in range(c - l, c + l + 1):
                            if j < 0 or j >= self.n:
                                continue
                            to_id = i * self.n  + j
                            if from_id != to_id:
                                self.g.add_edge(from_id, to_id)

    def __str__(self):
        """
        Return a string representation for this network.
        """

        return "KNN Network (n = %d, n = %d, k = %d, periodic = %s)" \
            %(len(self.vertices), self.n, self.k, self.periodic)

class Watts_Strogatz_Network(Network):
    """
    Representation of a Watts Strogatz network.
    """

    def __init__(self, vertices, k, p, seed = int(time.time())):
        """
        Construct a Watts Strogatz Network. The vertices argument contains 
        objects of type Vertex representing the vertices of the network. Each 
        vertex is connected to k nearest neighbors in ring topology. The p 
        argument specifies the probability of rewiring each edge. The seed 
        argument is the seed for the random number generator used to build 
        the network. 
        """

        self.vertices = vertices
        self.k = k
        self.p = p
        self.seed = seed
        self.g = networkx.watts_strogatz_graph(len(vertices), k, p, 
                                               seed = seed)
        self.__adj_list_map__ = {}

    def __str__(self):
        """
        Return a string representation for this network.
        """

        return "Watts Strogatz Network (n = %d, k = %d, p = %f)" \
            %(len(self.vertices), self.k, self.p)

class GraphML_File_Network(Network):
    """
    Representation of a network constructed from a GraphML (Graph Markup 
    Language) file.
    """
    
    def __init__(self, vertices, fname):
        """
        Construct a network from the GraphML file specified by fname. The 
        vertices argument contains objects of type Vertex representing the 
        vertices of the network. The fname argument is the name of the 
        GraphML file.
        """

        self.vertices = vertices
        self.fname = fname
        self.g = networkx.convert_node_labels_to_integers(\
            networkx.read_graphml(fname))
        self.__adj_list_map__ = {}

    def __str__(self):
        """
        Return a string representation for this network.
        """

        return "GraphML File Network (n = %d, file = %s)" \
            %(len(self.vertices), self.fname)

def build_network(vertices, topology, params):
    """
    Build and return the requested network. The vertices argument contains 
    objects of type Vertex representing the vertices of the network. The 
    topology argument specifies the type of network that is to be built; 
    allowed topologies are:

      Complete, Erdos_Renyi, Barabasi_Albert, Powerlaw_Homophilic, 
      Powerlaw_Clustered, Random_Regular, Grid_2D, KNN, Watts_Strogatz, and 
      GraphML_File

    The params argument is a map specifying the parameters for the 
    requested network topology. Allowed parameters are:

      Complete: {}
      Erdos_Renyi: {"p": ..., "seed": ...}
      Barabasi_Albert: {"m": ..., "seed": ...}
      Powerlaw_Homophilic: {"m": ..., "h": ..., "maxiter": ..., 
                    "seed": ...}
      Powerlaw_Clustered: {"m": ..., "p": ..., "seed": ...}
      Random_Regular: {"degree": ..., "seed": ...}
      Grid_2D: {"m": ..., "n": ..., "periodic": ...}
      KNN: {"n": ..., "k": ..., "periodic": ...}
      Watts_Strogatz: {"k": ..., "p": ..., "seed": ...}
      GraphML_File: {"name": ...}
    """

    net = None
    if topology == "Complete":
        net = Complete_Network(vertices)
    elif topology == "Erdos_Renyi":
        net = Erdos_Renyi_Network(vertices, params["p"], params["seed"])
    elif topology == "Barabasi_Albert":
        net = Barabasi_Albert_Network(vertices, params["m"], params["seed"])
    elif topology == "Powerlaw_Homophilic":
        net = Powerlaw_Homophilic_Network(vertices, params["m"], 
                                           params["h"], params["maxiter"], 
                                           params["seed"])
    elif topology == "Powerlaw_Clustered":
        net = Powerlaw_Clustered_Network(vertices, params["m"], params["p"], 
                                       params["seed"])
    elif topology == "Random_Regular":
        net = Random_Regular_Network(vertices, params["degree"], 
                                     params["seed"])
    elif topology == "Grid_2D":
        net = Grid_2D_Network(vertices, params["m"], params["n"], 
                              params["periodic"])
    elif topology == "KNN":
        net = KNN_Network(vertices, params["n"], params["k"], 
                          params["periodic"])
    elif topology  == "Watts_Strogatz":
        net = Watts_Strogatz_Network(vertices, params["k"], params["p"], 
                                     params["seed"])
    elif topology == "GraphML_File":
        net = GraphML_File_Network(vertices, params["name"])
    return net

