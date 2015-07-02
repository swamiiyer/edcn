// Discrete.java

// Abstraction for a module for carrying out the game dynamics. 
abstract class Dynamics
{
    // The network on which the game is played.
    protected Network net;

    // Game parameters.
    protected Params params;

    // Traits of individuals in a population.
    protected int[] traits;

    // Copy of traits of individuals in a population. These values get 
    // commited to traits at the end of each round of updates.
    protected int[] traitsCopy;

    // Payoffs of individuals in the population.
    protected double[] payoffs;

    // Fitnesses of individuals in the population.
    protected double[] fitnesses;

    // Stores cumulative distribution of fitnesses, and is used in 
    // roulette selection.
    protected Pair[] fitnessCdf;

    // Total fitness the population.
    protected double totalFitness;

    // Min payoff in the population.
    protected double minPayoff;

    // Max payoff in the population.
    protected double maxPayoff;

    // Encapsulates the id of an individual in the population and the 
    // total fitness of the population upto and including the individual.
    protected class Pair
    {
	// Id of an individual.
	public int id;

	// Total fitness of the population upto and including this individual.
	public double sum;
	
	// Construct an instance of Pair given the id of an individual and 
	// the total fitness of the population upto and including the 
	// individual.
	public Pair(int id, double sum)
	{
	    this.id = id;
	    this.sum = sum;
	}
    }

    // Construct an instance of Dynamics given the network, game parameters, 
    // and the traits of the individuals in a population.
    protected Dynamics(Network net, Params params, int[] traits)
    {
	this.net = net;
	this.params = params;
	this.traits = traits;
	traitsCopy = this.traits.clone();
	payoffs = new double[params.population];
	fitnesses = new double[params.population];
    }

    // Invoked immediately before the round of interactions. Reinitialize 
    // the min and max payoff of the population.
    public void preInteraction()
    {
	minPayoff = Double.POSITIVE_INFINITY;
	maxPayoff = Double.NEGATIVE_INFINITY;
    }
    
    // Carry out an interaction.
    public abstract void interact();

    // Computations that need to be performed immediately after the round of 
    // interactions must be carried out here.
    public abstract void postInteraction();

    // Carry out an update.
    public abstract void update();

    // Called after the round of updates. Commit inheritance of trait for 
    // each individual in the population.
    public void postUpdate()
    {
	for (int i = 0; i < traitsCopy.length; i++) {
	    traits[i] = traitsCopy[i];
	}
    }
    
    // Return the fitness given the payoff p.
    protected double fitness(double p) 
    {
	//double fitness = Math.exp(params.selection_strength * p);
	double fitness = 
	    1.0 - params.selection_strength + params.selection_strength * p;
	assert(fitness >= 0.0);
	return fitness;
    }
    
    // Returns the probability that i will inherit j's trait given their 
    // payoffs pi and pj.
    protected double replicate(double pi, double pj)
    {
        if (pi >= pj) {
            return 0.0;
	}
        return params.selection_strength * 
            (pj - pi) / (maxPayoff - minPayoff);
    }

    // Returns the probability that i will inherit j's trait. pi and 
    // pj are the payoffs of i and j.
    protected double fermi(double pi, double pj)
    {
        return 1.0 / (1.0 + Math.exp(-params.selection_strength * 
				     (pj - pi)));
    }

    // Given the cumulative distribution l of fitness values and the total 
    // fitness of the population, using binary search pick and return an 
    // individual from l in proportion to fitness.
    protected int roulette(Pair[] l, double totalFitness)
    {
        double r = Random.uniform() * totalFitness;
        int lo = 0, hi = l.length - 1, mid = -1;
        while (lo <= hi) {
            mid = lo + (hi - lo) / 2;
            if (r > l[mid].sum) {
                lo = mid + 1;
	    }
            else if (r < l[mid].sum) {
                hi = mid - 1;
	    }
	}
        return l[mid].id;
    }

    // Rreturn the appropriate dynamics module given the network, game 
    // parameters, and the traits of the individuals in a population.
    public static Dynamics module(Network net, Params params, int[] traits)
    {
	Dynamics dynamics = null;
	if (params.network_topology.equals("Complete")) {
	    if (params.update_rule.equals("RE1")) {
		dynamics = new UnstructuredReplication1(net, params, traits); 
	    }
	    else if (params.update_rule.equals("RE2")) {
		dynamics = new UnstructuredReplication2(net, params, traits); 
	    }
	    else if (params.update_rule.equals("FE1")) {
		dynamics = new UnstructuredFermi1(net, params, traits); 
	    }
	    else if (params.update_rule.equals("FE2")) {
		dynamics = new UnstructuredFermi2(net, params, traits); 
	    }	
	    else if (params.update_rule.equals("IM")) {
		dynamics = new UnstructuredImitation(net, params, traits); 
	    }	
	    else if (params.update_rule.equals("BD")) {
		dynamics = new UnstructuredBirthDeath(net, params, traits); 
	    }	
	    else if (params.update_rule.equals("DB")) {
		dynamics = new UnstructuredDeathBirth(net, params, traits); 
	    }	
	}
	else {
	    if (params.update_rule.equals("RE1")) {
		dynamics = new StructuredReplication1(net, params, traits); 
	    }
	    else if (params.update_rule.equals("RE2")) {
		dynamics = new StructuredReplication2(net, params, traits); 
	    }
	    else if (params.update_rule.equals("FE1")) {
		dynamics = new StructuredFermi1(net, params, traits); 
	    }
	    else if (params.update_rule.equals("FE2")) {
		dynamics = new StructuredFermi2(net, params, traits); 
	    }	
	    else if (params.update_rule.equals("IM")) {
		dynamics = new StructuredImitation(net, params, traits); 
	    }	
	    else if (params.update_rule.equals("BD")) {
		dynamics = new StructuredBirthDeath(net, params, traits); 
	    }	
	    else if (params.update_rule.equals("DB")) {
		dynamics = new StructuredDeathBirth(net, params, traits); 
	    }	
	}
	return dynamics;
    }
}

// Replication dynamics on an unstructured network.
class UnstructuredReplication1 extends Dynamics
{
    // Build an instance of UnstructuredReplication1 given the network, game 
    // parameters, and the traits of the individuals in a population.
    public UnstructuredReplication1(Network net, Params params, int[] traits)
    {
	super(net, params, traits);
    }

    // Carry out an interaction.
    public void interact()
    {
        int i = net.getRandomVertex();
	double payoff = 0.0;
        if (Random.uniform() < params.assortativity) {
            payoff = params.payoff[traits[i]][traits[i]];
	}
	else {
            int j = net.getRandomNeighbor(i);
            payoff = params.payoff[traits[i]][traits[j]];
	}
        payoffs[i] = payoff;
    }

    // Calculate the min and max payoff of the population.
    public void postInteraction()
    {
	minPayoff = Stats.min(payoffs);
	maxPayoff = Stats.max(payoffs);
    }

    // Carry out an update.
    public void update()
    {
        int i = net.getRandomVertex();
        int j = net.getRandomNeighbor(i);
        double p = replicate(payoffs[i], payoffs[j]);
        if (Random.uniform() < p) {
	    traitsCopy[i] = traits[j];
	}
    }
}

// Replication dynamics on an unstructured network.
class UnstructuredReplication2 extends Dynamics
{
    // Build an instance of UnstructuredReplication2 given the network, game 
    // parameters, and the traits of the individuals in a population.
    public UnstructuredReplication2(Network net, Params params, int[] traits)
    {
	super(net, params, traits);
    }

    // Carry out an interaction followed by an update.
    public void interact()
    {
	// Interact.
        int i = net.getRandomVertex();
        int j = net.getRandomNeighbor(i);
	double pi = 0.0, pj = 0.0;
        if (Random.uniform() < params.assortativity) {
            pi = params.payoff[traits[i]][traits[i]];
	}
        else {
            int k = net.getRandomNeighbor(i);
            pi = params.payoff[traits[i]][traits[k]];
	}
        if (Random.uniform() < params.assortativity) {
            pj = params.payoff[traits[j]][traits[j]];
	}
        else {
            int l = net.getRandomNeighbor(j);
            pj = params.payoff[traits[j]][traits[l]];
	}
        payoffs[i] = pi;
        payoffs[j] = pj;
        minPayoff = Stats.min(new double[] {minPayoff, pi, pj});
        maxPayoff = Stats.max(new double[] {maxPayoff, pi, pj});

        // Update.
	double p = replicate(pi, pj);
        if (Random.uniform() < p) {
            traitsCopy[i] = traits[j];
	}
    }

    // Nothing to do here.
    public void postInteraction()
    {
	return;
    }

    // Nothing to do here.
    public void update()
    {
	return;
    }
}

// Fermi dynamics on an unstructured network.
class UnstructuredFermi1 extends Dynamics
{
    // Build an instance of UnstructuredFermi1 given the network, game 
    // parameters, and the traits of the individuals in a population.
    public UnstructuredFermi1(Network net, Params params, int[] traits)
    {
	super(net, params, traits);
    }

    // Carry out an interaction.
    public void interact()
    {
        int i = net.getRandomVertex();
	double payoff = 0.0;
        if (Random.uniform() < params.assortativity) {
            payoff = params.payoff[traits[i]][traits[i]];
	}
	else {
            int j = net.getRandomNeighbor(i);
            payoff = params.payoff[traits[i]][traits[j]];
	}
        payoffs[i] = payoff;
    }

    // Nothing to be do here.
    public void postInteraction()
    {
	return;
    }

    // Carry out an update.
    public void update()
    {
        int i = net.getRandomVertex();
        int j = net.getRandomNeighbor(i);
        double p = fermi(payoffs[i], payoffs[j]);
        if (Random.uniform() < p) {
	    traitsCopy[i] = traits[j];
	}
    }
}

// Fermi dynamics on an unstructured network.
class UnstructuredFermi2 extends Dynamics
{
    // Build an instance of UnstructuredFermi2 given the network, game 
    // parameters, and the traits of the individuals in a population.
    public UnstructuredFermi2(Network net, Params params, int[] traits)
    {
	super(net, params, traits);
    }

    // Carry out an interaction followed by an update.
    public void interact()
    {
	// Interact.
        int i = net.getRandomVertex();
        int j = net.getRandomNeighbor(i);
	double pi = 0.0, pj = 0.0;
        if (Random.uniform() < params.assortativity) {
            pi = params.payoff[traits[i]][traits[i]];
	}
        else {
            int k = net.getRandomNeighbor(i);
            pi = params.payoff[traits[i]][traits[k]];
	}
        if (Random.uniform() < params.assortativity) {
            pj = params.payoff[traits[j]][traits[j]];
	}
        else {
            int l = net.getRandomNeighbor(j);
            pj = params.payoff[traits[j]][traits[l]];
	}
        payoffs[i] = pi;
        payoffs[j] = pj;

        // Update.
	double p = fermi(pi, pj);
        if (Random.uniform() < p) {
            traitsCopy[i] = traits[j];
	}
    }

    // Nothing to do here.
    public void postInteraction()
    {
	return;
    }

    // Nothing to do here.
    public void update()
    {
	return;
    }
}

// Imitation dynamics on an unstructured network.
class UnstructuredImitation extends Dynamics
{
    // Build an instance of UnstructuredImitation given the network, game 
    // parameters, and the traits of the individuals in a population.
    public UnstructuredImitation(Network net, Params params, int[] traits)
    {
	super(net, params, traits);
    }

    // Carry out an interaction.
    public void interact()
    {
        int i = net.getRandomVertex();
	double payoff = 0.0;
        if (Random.uniform() < params.assortativity) {
            payoff = params.payoff[traits[i]][traits[i]];
	}
	else {
            int j = net.getRandomNeighbor(i);
            payoff = params.payoff[traits[i]][traits[j]];
	}
        payoffs[i] = payoff;
	fitnesses[i] = fitness(payoff);
    }

    // Calculate the fitness CDF and the total fitness of the population.
    public void postInteraction()
    {

        totalFitness = 0.0;
        fitnessCdf = new Pair[net.size()];
        for (int i = 0; i < net.size(); i++) {
            totalFitness += fitnesses[i];
	    fitnessCdf[i] = new Pair(i, totalFitness);
	}
    }

    // Carry out an update.
    public void update()
    {
        int i = net.getRandomVertex();
        int j = roulette(fitnessCdf, totalFitness);
	traitsCopy[i] = traits[j];
    }
}

// Birth-death dynamics on an unstructured network.
class UnstructuredBirthDeath extends Dynamics
{
    // Build an instance of UnstructuredBirthDeath given the network, game 
    // parameters, and the traits of the individuals in a population.
    public UnstructuredBirthDeath(Network net, Params params, int[] traits)
    {
	super(net, params, traits);
    }

    // Carry out an interaction.
    public void interact()
    {
        int i = net.getRandomVertex();
	double payoff = 0.0;
        if (Random.uniform() < params.assortativity) {
            payoff = params.payoff[traits[i]][traits[i]];
	}
	else {
            int j = net.getRandomNeighbor(i);
            payoff = params.payoff[traits[i]][traits[j]];
	}
        payoffs[i] = payoff;
	fitnesses[i] = fitness(payoff);
    }

    // Calculate the fitness CDF and total fitness of the population.
    public void postInteraction()
    {
        totalFitness = 0.0;
        fitnessCdf = new Pair[net.size()];
        for (int i = 0; i < net.size(); i++) {
            totalFitness += fitnesses[i];
	    fitnessCdf[i] = new Pair(i, totalFitness);
	}
    }

    // Carry out an update.
    public void update()
    {
        int j = roulette(fitnessCdf, totalFitness);
        int i = net.getRandomNeighbor(j);
	traitsCopy[i] = traits[j];
    }
}

// Death-birth dynamics on an unstructured network.
class UnstructuredDeathBirth extends Dynamics
{
    // Build an instance of UnstructuredDeathBirth given the network, game 
    // parameters, and the traits of the individuals in a population.
    public UnstructuredDeathBirth(Network net, Params params, int[] traits)
    {
	super(net, params, traits);
    }

    // Carry out an interaction.
    public void interact()
    {
        int i = net.getRandomVertex();
	double payoff = 0.0;
        if (Random.uniform() < params.assortativity) {
            payoff = params.payoff[traits[i]][traits[i]];
	}
	else {
            int j = net.getRandomNeighbor(i);
            payoff = params.payoff[traits[i]][traits[j]];
	}
        payoffs[i] = payoff;
	fitnesses[i] = fitness(payoff);
    }

    // Calculate the fitness CDF and the total fitness of the population.
    public void postInteraction()
    {

        totalFitness = 0.0;
        fitnessCdf = new Pair[net.size()];
        for (int i = 0; i < net.size(); i++) {
            totalFitness += fitnesses[i];
	    fitnessCdf[i] = new Pair(i, totalFitness);
	}
    }

    // Carry out an update.
    public void update()
    {
        int i = net.getRandomVertex();
	int j = roulette(fitnessCdf, totalFitness);
	traitsCopy[i] = traits[j];
    }
}

// Replication dynamics on a structured network.
class StructuredReplication1 extends Dynamics
{
    // Build an instance of StructuredReplication1 given the network, game 
    // parameters, and the traits of the individuals in a population.
    public StructuredReplication1(Network net, Params params, int[] traits)
    {
	super(net, params, traits);
    }

    // Carry out an interaction.
    public void interact()
    {
        int i = net.getRandomVertex();
        int j = net.getRandomNeighbor(i);
        if (j == -1) {
            return;
	}
        double payoff = params.payoff[traits[i]][traits[j]];
	payoffs[i] = payoff;
    }

    // Calculate the min and max payoff of the population.
    public void postInteraction()
    {
	minPayoff = Stats.min(payoffs);
	maxPayoff = Stats.max(payoffs);
    }

    // Carry out an interaction.
    public void update()
    {
        int i = net.getRandomVertex();
        int j = net.getRandomNeighbor(i);
        if (j == -1) {
            return;
	}
        double p = replicate(payoffs[i], payoffs[j]);
        if (Random.uniform() < p) { 
	    traitsCopy[i] = traits[j];
	}
    }
}

// Replication dynamics on a structured network.
class StructuredReplication2 extends Dynamics
{
    // Build an instance of StructuredReplication2 given the network, game 
    // parameters, and the traits of the individuals in a population.
    public StructuredReplication2(Network net, Params params, int[] traits)
    {
	super(net, params, traits);
    }

    // Carry out an interaction followed by an update.
    public void interact()
    {
        // Interact.
        int i = net.getRandomVertex();
        int j = net.getRandomNeighbor(i);
	if (j == -1) {
	    return;
	}
        int k = net.getRandomNeighbor(i);
        int l = net.getRandomNeighbor(j);
        if (k == -1 || l == -1 || i == l || j == k) {
            return;
	}
        double pi = params.payoff[traits[i]][traits[k]];
        double pj = params.payoff[traits[j]][traits[l]];
	payoffs[i] = pi;
	payoffs[j] = pj;
        minPayoff = Stats.min(new double[] {minPayoff, pi, pj});
        maxPayoff = Stats.max(new double[] {maxPayoff, pi, pj});

        // Update.
        double p = replicate(pi, pj);
        if (Random.uniform() < p) {
            traitsCopy[i] = traits[j];
	}
    }

    // Nothing to do here.
    public void postInteraction()
    {
	return;
    }

    // Nothing to do here.
    public void update()
    {
	return;
    }
}

// Fermi dynamics on a structured network.
class StructuredFermi1 extends Dynamics
{
    // Build an instance of StructuredFermi1 given the network, game 
    // parameters, and the traits of the individuals in a population.
    public StructuredFermi1(Network net, Params params, int[] traits)
    {
	super(net, params, traits);
    }

    // Carry out an interaction.
    public void interact()
    {
        int i = net.getRandomVertex();
        int j = net.getRandomNeighbor(i);
        if (j == -1) {
            return;
	}
        double payoff = params.payoff[traits[i]][traits[j]];
	payoffs[i] = payoff;
    }

    // Nothing to do here.
    public void postInteraction()
    {
	return;
    }

    // Carry out an update.
    public void update()
    {
        int i = net.getRandomVertex();
        int j = net.getRandomNeighbor(i);
        if (j == -1) {
            return;
	}
        double p = fermi(payoffs[i], payoffs[j]);
        if (Random.uniform() < p) { 
	    traitsCopy[i] = traits[j];
	}
    }
}

// Fermi dynamics on a structured network.
class StructuredFermi2 extends Dynamics
{
    // Build an instance of StructuredFermi2 given the network, game 
    // parameters, and the traits of the individuals in a population.
    public StructuredFermi2(Network net, Params params, int[] traits)
    {
	super(net, params, traits);
    }

    // Carry out an interaction followed by an update.
    public void interact()
    {
        // Interact.
        int i = net.getRandomVertex();
        int j = net.getRandomNeighbor(i);
	if (j == -1) {
	    return;
	}
        int k = net.getRandomNeighbor(i);
        int l = net.getRandomNeighbor(j);
        if (k == -1 || l == -1 || i == l || j == k) {
            return;
	}
        double pi = params.payoff[traits[i]][traits[k]];
        double pj = params.payoff[traits[j]][traits[l]];
	payoffs[i] = pi;
	payoffs[j] = pj;

        // Update.
        double p = fermi(pi, pj);
        if (Random.uniform() < p) {
            traitsCopy[i] = traits[j];
	}
    }

    // Nothing to do here.
    public void postInteraction()
    {
	return;
    }

    // Nothing to do here.
    public void update()
    {
	return;
    }
}

// Imitation dynamics on a structured network.
class StructuredImitation extends Dynamics
{
    // Build an instance of StructuredImitation given the network, game 
    // parameters, and the traits of the individuals in a population.
    public StructuredImitation(Network net, Params params, int[] traits)
    {
	super(net, params, traits);
    }

    // Carry out an interaction.
    public void interact()
    {
        int i = net.getRandomVertex();
        int j = net.getRandomNeighbor(i);
        if (j == -1) {
            return;
	}
        double payoff = params.payoff[traits[i]][traits[j]];
	payoffs[i] = payoff;
	fitnesses[i] = fitness(payoff);
    }

    // Nothing to do here.
    public void postInteraction()
    {
	return;
    }

    // Carry out an interaction.
    public void update()
    {
        int i = net.getRandomVertex();
        double totalFitness = 0.0;
	int[] neighbors = net.getNeighbors(i);
        Pair[] fitnessCdf = new Pair[neighbors.length + 1];
	for (int j = 0; j < neighbors.length; j++) {
	    totalFitness += fitnesses[neighbors[j]];
            fitnessCdf[j] = new Pair(neighbors[j], totalFitness);
	}
	totalFitness += fitnesses[i];
	fitnessCdf[neighbors.length] = new Pair(i, totalFitness);
        int j = roulette(fitnessCdf, totalFitness);
	traitsCopy[i] = traits[j];
    }
}

// Birth-death dynamics on a structured network.
class StructuredBirthDeath extends Dynamics
{
    // Build an instance of StructuredBirthDeath given the network, game 
    // parameters, and the traits of the individuals in a population.
    public StructuredBirthDeath(Network net, Params params, int[] traits)
    {
	super(net, params, traits);
    }

    // Carry out an interaction.
    public void interact()
    {
        int i = net.getRandomVertex();
        int j = net.getRandomNeighbor(i);
        if (j == -1) {
            return;
	}
        double payoff = params.payoff[traits[i]][traits[j]];
	payoffs[i] = payoff;
	fitnesses[i] = fitness(payoff);
    }

    // Calculate the fitness CDF and the total fitness of the population.
    public void postInteraction()
    {
        totalFitness = 0.0;
        fitnessCdf = new Pair[net.size()];
        for (int i = 0; i < net.size(); i++) {
            totalFitness += fitnesses[i];
	    fitnessCdf[i] = new Pair(i, totalFitness);
	}
    }

    // Carry out an update.
    public void update()
    {
        int j = roulette(fitnessCdf, totalFitness);
        int i = net.getRandomNeighbor(j);
        if (i == -1) {
            return;
	}
        traitsCopy[i] = traits[j];
    }
}

// Death-birth dynamics on a structured network.
class StructuredDeathBirth extends Dynamics
{
    // Build an instance of StructuredDeathBirth given the network, game 
    // parameters, and the traits of the individuals in a population.
    public StructuredDeathBirth(Network net, Params params, int[] traits)
    {
	super(net, params, traits);
    }

    // Carry out an interaction.
    public void interact()
    {
        int i = net.getRandomVertex();
        int j = net.getRandomNeighbor(i);
        if (j == -1) {
            return;
	}
        double payoff = params.payoff[traits[i]][traits[j]];
	payoffs[i] = payoff;
	fitnesses[i] = fitness(payoff);
    }

    // Nothing to do here.
    public void postInteraction()
    {
	return;
    }

    // Carry out an update.
    public void update()
    {
        int i = net.getRandomVertex();
        double totalFitness = 0.0;
	int[] neighbors = net.getNeighbors(i);
        Pair[] fitnessCdf = new Pair[neighbors.length];
	for (int j = 0; j < neighbors.length; j++) {
	    totalFitness += fitnesses[neighbors[j]];
            fitnessCdf[j] = new Pair(neighbors[j], totalFitness);
	}
        int j = roulette(fitnessCdf, totalFitness);
	traitsCopy[i] = traits[j];
    }
}

// Class that implements the agent-based model for simulating discrete 2x2 
// games.
public class Discrete
{
    // Every individual in the population possesses one of these two traits.
    private static int C = 0;
    private static int D = 1;

    // Carry out a single trial of the game given the parameters, and return 
    // an array containing the number of cooperators at each (saved) 
    // generation.
    private static int[] trial(Params params)
    {
	// Traits of individuals in the population.
	int[] traits = new int[params.population];
	
        // Start out with a population with the specified initial fraction 
	// of the population playing trait C and the remaining fraction 
        // playing trait D.
	for (int i = 0; i < params.population; i++) {
	    traits[i] = (Random.uniform() < params.init_cooperators)? 
		Discrete.C : Discrete.D;
	}

	// The underlying network structure.
	Network net;
	if (params.network_topology.equals("Complete")) {
	    net = new CompleteNetwork(params.population);
	}
	else {
	    net = new GraphMLNetwork(params.graphml_file);
	}

	// Create appropriate dynamics module.
	Dynamics dynamics = Dynamics.module(net, params, traits);

	// The dynamics.
	int[] nCs = new int[params.generations / params.report_freq + 1];
	int index = 0;
	for (int t = 0; t <= params.generations; t++) {

	    // Pre interaction.
	    dynamics.preInteraction();

            // Save the number of cooperators in this generation and 
	    // Plot results at report_freq and stop simulation if 
            // population has fixated to a trait.
	    if (t % params.report_freq == 0) {
		int nC = 0;
		for (int p = 0; p < params.population; p++) {
		    if (traits[p] == Discrete.C) {
			nC += 1;
		    }
		}
		nCs[index++] = nC;
		if (nC == params.population) {
		    for (int i = index; i < nCs.length; i++) {
			nCs[i] = params.population;
		    }
		    return nCs;
		}
		else if (nC == 0) {
		    return nCs;
		}
	    }

	    // Interact.
	    for (int count = 0; count < params.population; count++) {
		dynamics.interact();
	    }

	    // Post interaction.
	    dynamics.postInteraction();

	    // Update.
	    for (int count = 0; count < params.population; count++) {
		dynamics.update();
	    }

	    // Post update; commit trait inheritance.
	    dynamics.postUpdate();
	}   

	return nCs;
    }

    // Entry point.
    public static void main(String[] args)
    {
	// Parse command-line and load game parameters.
	if (args.length != 1) {
	    System.out.println("Usage: java Discrete <params file>");
	    System.exit(1);
	}
	Params params = Params.loadParams(args[0]);

	// Play the game the specified number of times and calculate 
	// the average fraction of cooperators in each (saved) generation.
	double[] aCs = new double[params.generations / params.report_freq + 1];
	for (int trial = 1; trial <= params.trials; trial++) {
	    int[] nCs = trial(params);
	    for (int i = 0; i < nCs.length; i++) {
		aCs[i] += nCs[i];
	    }
	}
	for (int i = 0; i < aCs.length; i++) {
	    aCs[i] /= (params.trials * params.population);
	}

	// Print to STDOUT the popultion size, the number of generations, 
	// report frequency, and the average fraction of cooperators 
	// over the last 10% of the (saved) generations. 
	int last = (int) (params.generations * 0.1);
	int b = params.generations / params.report_freq - 
	    last / params.report_freq;
	int e = params.generations / params.report_freq;
	System.out.printf("%d %d %d %f\n", 
			  params.population, params.generations, 
			  params.report_freq, Stats.mean(aCs, b, e));

	// Print to STDOUT the average fraction of cooperators in each of 
	// the (saved) generation.
	for (double aC : aCs) {
	    System.out.printf("%5.4f\n", aC);
	}
    }
}