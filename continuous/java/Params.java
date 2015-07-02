// Params.java

import com.google.gson.Gson;
import java.io.FileNotFoundException;
import java.io.FileReader;

// Stores experimental parameters read from a json file.
public class Params
{
    // Name of the payoff function.
    public String payoff;

    // Parameter values for the payoff function.
    public double[] payoff_params;

    // Selection strength.
    public double selection_strength;

    // Total number of players. 
    public int population;  

    // Number of generations. 
    public int generations; 

    // Frequency (in generations) at which results will be saved. 
    public int report_freq;

    // Initial trait for each player.
    public double init_trait;

    // Maximum value that the trait of any player can attain.
    public double max_trait;
    
    // The probability with which a player interacts with another player 
    // having the same trait as itself.
    public double assortativity;

    // group specifies the number of neighbors each player 
    // interacts with during each round of interactions. If 0, each 
    // player will interact with its entire neighborhood. The average of 
    // the neighbors' traits is used in calculating the payoff. Note 
    // that this setting is not applicable in "FE2" and "RE2" dynamics,  
    // and must be set to 1 if "assortativity" is positive.
    public int group;  

    // Update rule ("RE1", "RE2", "FE1", "FE2", "IM", "BD", "DB").
    public String update_rule; 

    // Rate at which Gaussian mutations of traits are carried out.
    public double mutation;

    // Standard deviation for the Gaussian mutations.
    public double stddev;  

    // Network topology ("Complete" or "GraphML_File").
    public String network_topology;

    // Filename if network topology is "GraphML_File".
    public String graphml_file;

    // Return a string representation of the experimental parameters.
    public String toString()
    {
	String m = "null"; 
	if (payoff_params != null) {
	    m = "[";
	    for (int i = 0; i < payoff_params.length - 1; i++) {
		m += payoff_params[i] + ", ";
	    }
	    m += payoff_params[payoff_params.length - 1] + "]";
	}
	return "payoff = " + payoff + "\n" + 
	    "payoff_params = " + m + "\n" + 
	    "selection_strength = " + selection_strength + "\n" +  
	    "population = " + population + "\n" + 
	    "generations = " + generations + "\n" + 
	    "report_freq = " + report_freq + "\n" + 
	    "init_trait = " + init_trait + "\n" + 
	    "max_trait = " + max_trait + "\n" + 
	    "assortativity = " + assortativity + "\n" + 
	    "group = " + group + "\n" + 
	    "update_rule = " + update_rule + "\n" + 
	    "mutation = " + mutation + "\n" + 
	    "stddev = " + stddev + "\n" + 
	    "network_topology = " + network_topology + "\n" +
	    "graphml_file = " + graphml_file;
    }

    // Returns a Params object built from the experimental parameters 
    // stored in the specified json file.
    public static Params loadParams(String fileName) 
    {
	Params params = null;
	try {
	    Gson gson = new Gson();
	    FileReader json = new FileReader(fileName);
	    params = gson.fromJson(json, Params.class);
	}
	catch (FileNotFoundException e) {
	    System.out.println(e.getMessage());
	}
	return params;
    }

    // Unit test.
    public static void main(String[] args)
    {
	System.out.println(loadParams(args[0]));
    }
}
