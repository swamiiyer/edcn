// Params.java

import com.google.gson.Gson;
import java.io.FileNotFoundException;
import java.io.FileReader;

// Stores experimental parameters read from a json file.
public class Params
{
    // Total number of players. 
    public int population;  

    // Number of generations. 
    public int generations; 

    // Frequency (in generations) at which results will be saved. 
    public int report_freq;

    // Infection probability for the starting infected individual. 
    // This may vary for different individuals in the infected state. 
    public double infection;

    // Virulence level c of a pathogen is calculated from its infection 
    // probability as c = pow(infection / b, 1.0 / a).
    public double b;
    public double a;

    // Recovery probability.
    public double recovery;
    
    // Natural death probability.
    public double death;

    // Rate at which Gaussian mutations of infection probability is 
    // carried out.
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
	return "population = " + population + "\n" + 
	    "generations = " + generations + "\n" + 
	    "report_freq = " + report_freq + "\n" + 
	    "infection = " + infection + "\n" + 
	    "b = " + b + "\n" + 
	    "a = " + a + "\n" + 
	    "recovery = " + recovery + "\n" + 
	    "death = " + death + "\n" + 
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
