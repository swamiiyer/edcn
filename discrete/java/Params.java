// Params.java

import com.google.gson.Gson;
import java.io.FileNotFoundException;
import java.io.FileReader;

// Stores experimental parameters read from a json file.
public class Params
{
    // A 2x2 payoff matrix.
    public double[][] payoff;

    // Selection strength.
    public double selection_strength;

    // Total number of players. 
    public int population;  

    // Initial fraction of population playing the cooperative trait C.
    public double init_cooperators;

    // Number of generations. 
    public int generations; 

    // Frequency (in generations) at which results will be saved. 
    public int report_freq;

    // Number of trials.
    public int trials;

    // The probability with which a player interacts with another player 
    // having the same trait as itself.
    public double assortativity;

    // Update rule ("RE1", "RE2", "FE1", "FE2", "IM", "BD", "DB").
    public String update_rule; 

    // Network topology ("Complete" or "GraphML_File").
    public String network_topology;

    // Filename if network topology is "GraphML_File".
    public String graphml_file;

    // Return a string representation of the experimental parameters.
    public String toString()
    {
	String m = (payoff == null) ? "null" : "[[" + 
	    payoff[0][0] + ", " + 
	    payoff[0][1] + "], [" + 
	    payoff[1][0] + ", " + 
	    payoff[1][1] + "]]";
	return "payoff = " + m + "\n" + 
	    "selection_strength = " + selection_strength + "\n" +  
	    "population = " + population + "\n" + 
	    "init_cooperators = " + init_cooperators + "\n" + 
	    "generations = " + generations + "\n" + 
	    "report_freq = " + report_freq + "\n" + 
	    "trials = " + trials + "\n" + 
	    "assortativity = " + assortativity + "\n" + 
	    "update_rule = " + update_rule + "\n" + 
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
