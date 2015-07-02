// Continuous.java

// Abstraction for a module for carrying out the game dynamics. 
abstract class Dynamics
{
    // The network on which the game is played.
    protected Network net;

    // Game parameters.
    protected Params params;

    // Traits of individuals in a population.
    protected double[] traits;

    // Copy of traits of individuals in a population. These values get 
    // commited to traits at the end of each round of updates.
    protected double[] traitsCopy;

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
    
    // Payoff function object.
    protected Payoff payoffRep;

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
    protected Dynamics(Network net, Params params, double[] traits)
    {
	this.net = net;
	this.params = params;
	this.traits = traits;
	traitsCopy = this.traits.clone();
	payoffs = new double[params.population];
	fitnesses = new double[params.population];
	try {
	    payoffRep = 
		(Payoff) Class.forName(params.payoff).newInstance();
	}
	catch (Exception e) {
	    System.out.println("Payoff object could not be instantiated!");
	    System.exit(1);
	}
	payoffRep.params = params.payoff_params;
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
    
    // Return the payoff given the traits x and y of two individuals.
    protected double payoff(double x, double y)
    {
	return payoffRep.payoff(x, y);
    }
    
    // Return the fitness given the payoff p.
    protected double fitness(double p) 
    {
	double fitness = Math.exp(params.selection_strength * p);
	// double fitness = 
	//     1.0 - params.selection_strength + params.selection_strength * p;
	// assert(fitness >= 0.0);
	return fitness;
    }
    
    // Returns the probability that i will inherit j's trait given their 
    // payoffs pi and pj.
    protected double replicate(double pi, double pj)
    {
        if (pi >= pj) {
            return 0.0;
	}
	return params.selection_strength * (pj - pi) / (maxPayoff - minPayoff);
    }

    // Returns the probability that i will inherit j's trait. pi and 
    // pj are the payoffs of i and j.
    protected double fermi(double pi, double pj)
    {
        return 1.0 / (1.0 + Math.exp(-params.selection_strength * (pj - pi)));
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
    public static Dynamics module(Network net, Params params, double[] traits)
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
    public UnstructuredReplication1(Network net, Params params, double[] traits)
    {
	super(net, params, traits);
    }

    // Carry out an interaction.
    public void interact()
    {
        int i = net.getRandomVertex();
	double payoff = 0.0;
        if (Random.uniform() < params.assortativity) {
            payoff = payoff(traits[i], traits[i]);
	}
	else {
	    double groupTrait = 0.0;
	    if (params.group == 0) {
		for (int j : net.getNeighbors(i)) {
		    groupTrait += traits[j];
		}
		groupTrait -= traits[i];
		// Comment the following line if the sum of the traits of 
		// neighbors is desired instead of average.
		// groupTrait /= (params.population - 1);
	    }
	    else {
		for (int count = 0; count < params.group; count++) {
		    int j = net.getRandomNeighbor(i);
		    groupTrait += traits[j];
		}
		// Comment the following line if the sum of the traits of 
		// neighbors is desired instead of average.
		// groupTrait /= params.group;
	    }
            payoff = payoff(traits[i], groupTrait);
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
	double trait = traits[i];
        if (Random.uniform() < p) {
	    trait = traits[j];
	}
	if (Random.uniform() < params.mutation) {
	    trait = Math.max(0.0, Math.min(Random.gaussian(trait, 
							   params.stddev),
					   params.max_trait));
	}
	traitsCopy[i] = trait;
    }
}

// Replication dynamics on an unstructured network.
class UnstructuredReplication2 extends Dynamics
{
    // Build an instance of UnstructuredReplication2 given the network, game 
    // parameters, and the traits of the individuals in a population.
    public UnstructuredReplication2(Network net, Params params, double[] traits)
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
            pi = payoff(traits[i], traits[i]);
	}
        else {
            int k = net.getRandomNeighbor(i);
            pi = payoff(traits[i], traits[k]);
	}
        if (Random.uniform() < params.assortativity) {
            pj = payoff(traits[j], traits[j]);
	}
        else {
            int l = net.getRandomNeighbor(j);
            pj = payoff(traits[j], traits[l]);
	}
        payoffs[i] = pi;
        payoffs[j] = pj;
        minPayoff = Stats.min(new double[] {minPayoff, pi, pj});
        maxPayoff = Stats.max(new double[] {maxPayoff, pi, pj});

        // Update.
	double p = replicate(pi, pj);
	double trait = traits[i];
        if (Random.uniform() < p) {
	    trait = traits[j];
	}
	if (Random.uniform() < params.mutation) {
	    trait = Math.max(0.0, Math.min(Random.gaussian(trait, 
							   params.stddev),
					   params.max_trait));
	}
	traitsCopy[i] = trait;
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
    public UnstructuredFermi1(Network net, Params params, double[] traits)
    {
	super(net, params, traits);
    }

    // Carry out an interaction.
    public void interact()
    {
        int i = net.getRandomVertex();
	double payoff = 0.0;
        if (Random.uniform() < params.assortativity) {
            payoff = payoff(traits[i], traits[i]);
	}
	else {
	    double groupTrait = 0.0;
	    if (params.group == 0) {
		for (int j : net.getNeighbors(i)) {
		    groupTrait += traits[j];
		}
		groupTrait -= traits[i];
		// Comment the following line if the sum of the traits of 
		// neighbors is desired instead of average.
		// groupTrait /= (params.population - 1);
	    }
	    else {
		for (int count = 0; count < params.group; count++) {
		    int j = net.getRandomNeighbor(i);
		    groupTrait += traits[j];
		}
		// Comment the following line if the sum of the traits of 
		// neighbors is desired instead of average.
		// groupTrait /= params.group;
	    }
            payoff = payoff(traits[i], groupTrait);
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
	double trait = traits[i];
        if (Random.uniform() < p) {
	    trait = traits[j];
	    if (Random.uniform() < params.mutation) {
		trait = Math.max(0.0, Math.min(Random.gaussian(trait, 
							       params.stddev),
					       params.max_trait));
	    }
	}
	traitsCopy[i] = trait;
    }
}

// Fermi dynamics on an unstructured network.
class UnstructuredFermi2 extends Dynamics
{
    // Build an instance of UnstructuredFermi2 given the network, game 
    // parameters, and the traits of the individuals in a population.
    public UnstructuredFermi2(Network net, Params params, double[] traits)
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
            pi = payoff(traits[i], traits[i]);
	}
        else {
            int k = net.getRandomNeighbor(i);
            pi = payoff(traits[i], traits[k]);
	}
        if (Random.uniform() < params.assortativity) {
            pj = payoff(traits[j], traits[j]);
	}
        else {
            int l = net.getRandomNeighbor(j);
            pj = payoff(traits[j], traits[l]);
	}
        payoffs[i] = pi;
        payoffs[j] = pj;

        // Update.
	double p = fermi(pi, pj);
	double trait = traits[i];
        if (Random.uniform() < p) {
	    trait = traits[j];
	    if (Random.uniform() < params.mutation) {
		trait = Math.max(0.0, Math.min(Random.gaussian(trait, 
							       params.stddev),
					       params.max_trait));
	    }
	}
	traitsCopy[i] = trait;
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
    public UnstructuredImitation(Network net, Params params, double[] traits)
    {
	super(net, params, traits);
    }

    // Carry out an interaction.
    public void interact()
    {
        int i = net.getRandomVertex();
	double payoff = 0.0;
        if (Random.uniform() < params.assortativity) {
            payoff = payoff(traits[i], traits[i]);
	}
	else {
	    double groupTrait = 0.0;
	    if (params.group == 0) {
		for (int j : net.getNeighbors(i)) {
		    groupTrait += traits[j];
		}
		groupTrait -= traits[i];
		// Comment the following line if the sum of the traits of 
		// neighbors is desired instead of average.
		// groupTrait /= (params.population - 1);
	    }
	    else {
		for (int count = 0; count < params.group; count++) {
		    int j = net.getRandomNeighbor(i);
		    groupTrait += traits[j];
		}
		// Comment the following line if the sum of the traits of 
		// neighbors is desired instead of average.
		// groupTrait /= params.group;
	    }
            payoff = payoff(traits[i], groupTrait);
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
	double trait = traits[j];
	if (Random.uniform() < params.mutation) {
	    trait = Math.max(0.0, Math.min(Random.gaussian(trait, 
							   params.stddev),
					   params.max_trait));
	}
	traitsCopy[i] = trait;
    }
}

// Birth-death dynamics on an unstructured network.
class UnstructuredBirthDeath extends Dynamics
{
    // Build an instance of UnstructuredBirthDeath given the network, game 
    // parameters, and the traits of the individuals in a population.
    public UnstructuredBirthDeath(Network net, Params params, double[] traits)
    {
	super(net, params, traits);
    }

    // Carry out an interaction.
    public void interact()
    {
        int i = net.getRandomVertex();
	double payoff = 0.0;
        if (Random.uniform() < params.assortativity) {
            payoff = payoff(traits[i], traits[i]);
	}
	else {
	    double groupTrait = 0.0;
	    if (params.group == 0) {
		for (int j : net.getNeighbors(i)) {
		    groupTrait += traits[j];
		}
		groupTrait -= traits[i];
		// Comment the following line if the sum of the traits of 
		// neighbors is desired instead of average.
		// groupTrait /= (params.population - 1);
	    }
	    else {
		for (int count = 0; count < params.group; count++) {
		    int j = net.getRandomNeighbor(i);
		    groupTrait += traits[j];
		}
		// Comment the following line if the sum of the traits of 
		// neighbors is desired instead of average.
		// groupTrait /= params.group;
	    }
            payoff = payoff(traits[i], groupTrait);
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
	double trait = traits[j];
	if (Random.uniform() < params.mutation) {
	    trait = Math.max(0.0, Math.min(Random.gaussian(trait, 
							   params.stddev),
					   params.max_trait));
	}
	traitsCopy[i] = trait;
    }
}

// Death-birth dynamics on an unstructured network.
class UnstructuredDeathBirth extends Dynamics
{
    // Build an instance of UnstructuredDeathBirth given the network, game 
    // parameters, and the traits of the individuals in a population.
    public UnstructuredDeathBirth(Network net, Params params, double[] traits)
    {
	super(net, params, traits);
    }

    // Carry out an interaction.
    public void interact()
    {
        int i = net.getRandomVertex();
	double payoff = 0.0;
        if (Random.uniform() < params.assortativity) {
            payoff = payoff(traits[i], traits[i]);
	}
	else {
	    double groupTrait = 0.0;
	    if (params.group == 0) {
		for (int j : net.getNeighbors(i)) {
		    groupTrait += traits[j];
		}
		groupTrait -= traits[i];
		// Comment the following line if the sum of the traits of 
		// neighbors is desired instead of average.
		// groupTrait /= (params.population - 1);
	    }
	    else {
		for (int count = 0; count < params.group; count++) {
		    int j = net.getRandomNeighbor(i);
		    groupTrait += traits[j];
		}
		// Comment the following line if the sum of the traits of 
		// neighbors is desired instead of average.
		// groupTrait /= params.group;
	    }
            payoff = payoff(traits[i], groupTrait);
	}
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
	int j = roulette(fitnessCdf, totalFitness);
	double trait = traits[j];
	if (Random.uniform() < params.mutation) {
	    trait = Math.max(0.0, Math.min(Random.gaussian(trait, 
							   params.stddev),
					   params.max_trait));
	}
	traitsCopy[i] = trait;
    }
}

// Replication dynamics on a structured network.
class StructuredReplication1 extends Dynamics
{
    // Build an instance of StructuredReplication1 given the network, game 
    // parameters, and the traits of the individuals in a population.
    public StructuredReplication1(Network net, Params params, double[] traits)
    {
	super(net, params, traits);
    }

    // Carry out an interaction.
    public void interact()
    {
        int i = net.getRandomVertex();
	double groupTrait = 0.0;
	if (params.group == 0) {
	    for (int j : net.getNeighbors(i)) {
		groupTrait += traits[j];
	    }
	    // Comment the following line if the sum of the traits of 
	    // neighbors is desired instead of average.
	    // groupTrait /= net.getNeighbors(i).length;
	}
	else {
	    for (int count = 0; count < params.group; count++) {
		int j = net.getRandomNeighbor(i);
		if (j == -1) {
		    return;
		}
		groupTrait += traits[j];
		// Comment the following line if the sum of the traits of 
		// neighbors is desired instead of average.
		// groupTrait /= params.group;
	    }
	}
        double payoff = payoff(traits[i], groupTrait);
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
	double trait = traits[i];
        if (Random.uniform() < p) {
	    trait = traits[j];
	}
	if (Random.uniform() < params.mutation) {
	    trait = Math.max(0.0, Math.min(Random.gaussian(trait, 
							   params.stddev),
					   params.max_trait));
	}
	traitsCopy[i] = trait;
    }
}

// Replication dynamics on a structured network.
class StructuredReplication2 extends Dynamics
{
    // Build an instance of StructuredReplication2 given the network, game 
    // parameters, and the traits of the individuals in a population.
    public StructuredReplication2(Network net, Params params, double[] traits)
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
        double pi = payoff(traits[i], traits[k]);
        double pj = payoff(traits[j], traits[l]);
	payoffs[i] = pi;
	payoffs[j] = pj;
        minPayoff = Stats.min(new double[] {minPayoff, pi, pj});
        maxPayoff = Stats.max(new double[] {maxPayoff, pi, pj});

        // Update.
        double p = replicate(pi, pj);
	double trait = traits[i];
        if (Random.uniform() < p) {
	    trait = traits[j];
	}
	if (Random.uniform() < params.mutation) {
	    trait = Math.max(0.0, Math.min(Random.gaussian(trait, 
							   params.stddev),
					   params.max_trait));
	}
	traitsCopy[i] = trait;
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
    public StructuredFermi1(Network net, Params params, double[] traits)
    {
	super(net, params, traits);
    }

    // Carry out an interaction.
    public void interact()
    {
        int i = net.getRandomVertex();
	double groupTrait = 0.0;
	if (params.group == 0) {
	    for (int j : net.getNeighbors(i)) {
		groupTrait += traits[j];
	    }
	    // Comment the following line if the sum of the traits of 
	    // neighbors is desired instead of average.
	    // groupTrait /= net.getNeighbors(i).length;
	}
	else {
	    for (int count = 0; count < params.group; count++) {
		int j = net.getRandomNeighbor(i);
		if (j == -1) {
		    return;
		}
		groupTrait += traits[j];
		// Comment the following line if the sum of the traits of 
		// neighbors is desired instead of average.
		// groupTrait /= params.group;
	    }
	}
        double payoff = payoff(traits[i], groupTrait);
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
	double trait = traits[i];
        if (Random.uniform() < p) {
	    trait = traits[j];
	    if (Random.uniform() < params.mutation) {
		trait = Math.max(0.0, Math.min(Random.gaussian(trait, 
							       params.stddev),
					       params.max_trait));
	    }
	}
	traitsCopy[i] = trait;
    }
}

// Fermi dynamics on a structured network.
class StructuredFermi2 extends Dynamics
{
    // Build an instance of StructuredFermi2 given the network, game 
    // parameters, and the traits of the individuals in a population.
    public StructuredFermi2(Network net, Params params, double[] traits)
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
        double pi = payoff(traits[i], traits[k]);
        double pj = payoff(traits[j], traits[l]);
	payoffs[i] = pi;
	payoffs[j] = pj;

        // Update.
        double p = fermi(pi, pj);
	double trait = traits[i];
        if (Random.uniform() < p) {
	    trait = traits[j];
	    if (Random.uniform() < params.mutation) {
		trait = Math.max(0.0, Math.min(Random.gaussian(trait, 
							       params.stddev),
					       params.max_trait));
	    }
	}
	traitsCopy[i] = trait;
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
    public StructuredImitation(Network net, Params params, double[] traits)
    {
	super(net, params, traits);
    }

    // Carry out an interaction.
    public void interact()
    {
        int i = net.getRandomVertex();
	double groupTrait = 0.0;
	if (params.group == 0) {
	    for (int j : net.getNeighbors(i)) {
		groupTrait += traits[j];
	    }
	    // Comment the following line if the sum of the traits of 
	    // neighbors is desired instead of average.
	    // groupTrait /= net.getNeighbors(i).length;
	}
	else {
	    for (int count = 0; count < params.group; count++) {
		int j = net.getRandomNeighbor(i);
		if (j == -1) {
		    return;
		}
		groupTrait += traits[j];
		// Comment the following line if the sum of the traits of 
		// neighbors is desired instead of average.
		// groupTrait /= params.group;
	    }
	}
        double payoff = payoff(traits[i], groupTrait);
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
	double trait = traits[j];
	if (Random.uniform() < params.mutation) {
	    trait = Math.max(0.0, Math.min(Random.gaussian(trait, 
							   params.stddev),
					   params.max_trait));
	}
	traitsCopy[i] = trait;
    }
}

// Birth-death dynamics on a structured network.
class StructuredBirthDeath extends Dynamics
{
    // Build an instance of StructuredBirthDeath given the network, game 
    // parameters, and the traits of the individuals in a population.
    public StructuredBirthDeath(Network net, Params params, double[] traits)
    {
	super(net, params, traits);
    }

    // Carry out an interaction.
    public void interact()
    {
        int i = net.getRandomVertex();
	double groupTrait = 0.0;
	if (params.group == 0) {
	    for (int j : net.getNeighbors(i)) {
		groupTrait += traits[j];
	    }
	    // Comment the following line if the sum of the traits of 
	    // neighbors is desired instead of average.
	    // groupTrait /= net.getNeighbors(i).length;
	}
	else {
	    for (int count = 0; count < params.group; count++) {
		int j = net.getRandomNeighbor(i);
		if (j == -1) {
		    return;
		}
		groupTrait += traits[j];
		// Comment the following line if the sum of the traits of 
		// neighbors is desired instead of average.
		// groupTrait /= params.group;
	    }
	}
        double payoff = payoff(traits[i], groupTrait);
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
	double trait = traits[j];
	if (Random.uniform() < params.mutation) {
	    trait = Math.max(0.0, Math.min(Random.gaussian(trait, 
							   params.stddev),
					   params.max_trait));
	}
	traitsCopy[i] = trait;
    }
}

// Death-birth dynamics on a structured network.
class StructuredDeathBirth extends Dynamics
{
    // Build an instance of StructuredDeathBirth given the network, game 
    // parameters, and the traits of the individuals in a population.
    public StructuredDeathBirth(Network net, Params params, double[] traits)
    {
	super(net, params, traits);
    }

    // Carry out an interaction.
    public void interact()
    {
        int i = net.getRandomVertex();
	double groupTrait = 0.0;
	if (params.group == 0) {
	    for (int j : net.getNeighbors(i)) {
		groupTrait += traits[j];
	    }
	    // Comment the following line if the sum of the traits of 
	    // neighbors is desired instead of average.
	    // groupTrait /= net.getNeighbors(i).length;
	}
	else {
	    for (int count = 0; count < params.group; count++) {
		int j = net.getRandomNeighbor(i);
		if (j == -1) {
		    return;
		}
		groupTrait += traits[j];
		// Comment the following line if the sum of the traits of 
		// neighbors is desired instead of average.
		// groupTrait /= params.group;
	    }
	}
        double payoff = payoff(traits[i], groupTrait);
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
	double trait = traits[j];
	if (Random.uniform() < params.mutation) {
	    trait = Math.max(0.0, Math.min(Random.gaussian(trait, 
							   params.stddev),
					   params.max_trait));
	}
	traitsCopy[i] = trait;
    }
}

// Every payoff function must inherit this class and provide an implementation 
// for the abstract payoff function.
abstract class Payoff
{
    // Parameter values for the payoff function.
    public double[] params;
    
    // Return the payoff given the traits x and y.
    public abstract double payoff(double x, double y);
}

// Class that implements the agent-based model for simulating continuous games.
public class Continuous
{
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
	// Parse command-line and load game parameters.
	if (args.length != 1) {
	    System.out.println("Usage: java Continuous <params file>");
	    System.exit(1);
	}
	Params params = Params.loadParams(args[0]);

	// Traits of individuals in the population.
	double[] traits = new double[params.population];
	
        // Start out with a population with the specified initial fraction 
	// of the population playing trait C and the remaining fraction 
        // playing trait D.
	for (int i = 0; i < params.population; i++) {
	    traits[i] = params.init_trait;
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

	System.out.printf("%d %d %d %5.4f\n", params.population, 
			  params.generations, params.report_freq, 
			  params.max_trait);

	// The dynamics.
	for (int t = 0; t <= params.generations; t++) {

	    // Pre interaction.
	    dynamics.preInteraction();

	    // Gather and print some stats at t = report_freq. These include: 
	    // min, mean and max trait values, and the distribution of 
	    // traits.
	    if (t % params.report_freq == 0) {
		int[] dist = new int[101];
		for (int p = 0; p < params.population; p++) {
		    dist[(int) Math.floor(traits[p] / 
					  params.max_trait * 100)]++;
		}
		System.out.printf("%5.4f %5.4f %5.4f %s\n", Stats.min(traits), 
				  Stats.mean(traits), Stats.max(traits), 
				  toString(dist));
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
    }
}