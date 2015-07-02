// Disease.java

// Class that implements the agent-based virus dynamics model.
public class Disease
{
    // Every individual in the population is in one of these three states.
    private static int SUSCEPTIBLE = 0;
    private static int INFECTED = 1;
    private static int RECOVERED = 2;

    // Returns the virulence corresponding to the specified infectivity. 
    // The formula used is c = pow(infection / b, 1.0 / a).
    private static double fitness(double b, double a, double infection)
    {
	return Math.pow(infection / b, 1.0 / a);
    }
    
    // Returns a string representation of the items in the array a.
    private static String toString(int[] a)
    {
	String s = "";
	int i;
	for (i = 0; i < a.length - 1; i++) {
	    if (a[i] != 0) {
		s += i + ":" + a[i] + " "; 
	    }
	}
	s += (a[i] != 0) ? i + ":" + a[i] : "";
	return s;
    }

    // Entry point.
    public static void main(String[] args)
    {
	if (args.length != 1) {
	    System.out.println("Usage: java Disease <params file>");
	    System.exit(1);
	}
	Params params = Params.loadParams(args[0]);

	// States, infectivity values, and virulence values of individuals 
	// in the population.
	int[] states = new int[params.population];
	double[] infections = new double[params.population]; 
	double[] cs = new double[params.population];    

	// The underlying network structure.
	Network net;
	if (params.network_topology.equals("Complete")) {
	    net = new CompleteNetwork(params.population);
	}
	else {
	    net = new GraphMLNetwork(params.graphml_file);
	}

	// Infect an individual at random. This individual is the random 
	// neighbor of a randomly picked individual.
	int i, j;
	do {
	    i = net.getRandomVertex();
	}
	while ((j = net.getRandomNeighbor(i)) == -1);
	states[j] = INFECTED;
	infections[j] = params.infection;
	cs[j] = fitness(params.b, params.a, params.infection);

	System.out.printf("%d %d %d\n", params.population, params.generations, 
			  params.report_freq);

	// The dynamics.
	for (int t = 0; t <= params.generations; t++) {

	    // Gather and print some stats at t = report_freq. These include: 
	    // the number of susceptible and infected individuals, 
	    // the mean virulence level, and the distribution of virulence.
	    if (t % params.report_freq == 0) {
		int nS = 0;
		int nI = 0;
		int[] dist = new int[101];
		for (int p = 0; p < params.population; p++) {
		    if (states[p] == SUSCEPTIBLE) {
			nS += 1;
		    }
		    else if (states[p] == INFECTED) {
			nI += 1;
			dist[(int) Math.floor(cs[p] * 100)]++;
		    }
		}
		System.out.printf("%d %d %5.4f %s\n", nS, nI, 
				  Stats.sum(cs) / nI, 
				  toString(dist));
	    }

	    for (int p = 0; p < params.population; p++) {
		i = net.getRandomVertex();

		// i is susceptible.
		if (states[i] == SUSCEPTIBLE) {
		    if (Random.uniform() < params.death) {
			// i dies of natural cuases and is reborn as a 
			// susceptible individual.
			continue;
		    }
		    
		    // Find neighbor that could infect i. 
		    double probNot = 1.0;
		    double totalInfection = 0.0;
		    int[] neighbors = net.getNeighbors(i);
		    if (neighbors == null) {
			continue;
		    }
		    for (int k: neighbors) {
			if (states[k] == INFECTED) {
			    probNot *= (1.0 - infections[k]);
			    totalInfection += infections[k];
			}
		    }

		    if (Random.uniform() < 1.0 - probNot) {
			// i is infected by a neighbor, and takes on that 
			// individual's (possibly mutated) infectivity value.
			states[i] = INFECTED;
			double newInfection = 0.0;
			double r = Random.uniform() * totalInfection;
			double sum = 0.0;
			for (int k: neighbors) {
			    newInfection = infections[k];
			    sum += newInfection;
			    if (r < sum) {
				break;
			    }
			}
			if (Random.uniform() < params.mutation) {
			    newInfection = Math.max(0.0, 
					       Math.min(Random.
							gaussian(newInfection,
								 params.stddev),
							1.0));
			}
			infections[i] = newInfection;
			cs[i] = fitness(params.b, params.a, newInfection);
		    }
		}
		// i is infected.
		else if (states[i] == INFECTED) {
		    if (Random.uniform() < 
			1.0 - (1.0 - params.death) * (1.0 - cs[i])) {
			// i dies either of natural causes or due to the 
			// virus, and is reborn as a susceptible individual.
			states[i] = SUSCEPTIBLE;
			infections[i] = 0.0;
			cs[i] = 0.0;
		    }
		    else if (Random.uniform() < params.recovery) {
			// i recovers from the virus.
			states[i] = RECOVERED;
			infections[i] = 0.0;
			cs[i] = 0.0;
		    }
		}
		// i is recovered.
		else if (states[i] == RECOVERED) {
		    if (Random.uniform() < params.death) {
			// i dies of natural causes and is reborn as a 
			// susceptible individual.
			states[i] = SUSCEPTIBLE;
			infections[i] = 0.0;
			cs[i] = 0.0;
		    }
		}
	    }
	}
    }
}