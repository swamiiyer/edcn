// Network.java

import edu.uci.ics.jung.graph.*;
import edu.uci.ics.jung.io.GraphMLMetadata;
import edu.uci.ics.jung.io.GraphMLReader;
import java.io.IOException;
import javax.xml.parsers.ParserConfigurationException;
import org.xml.sax.SAXException;
import org.apache.commons.collections15.BidiMap;
import org.apache.commons.collections15.Factory;

// Every network must support this interface.
interface Network
{
    // Return the size of the network.
    public int size();

    // Return the id of a random vertex in the network.
    public int getRandomVertex();

    // Return ids of the neighbors of vertex i, or null.
    public int[] getNeighbors(int i);

    // Return the id of a random neighbor of vertex i, or -1.
    public int getRandomNeighbor(int i);
}

// An abstraction for a complete network.
class CompleteNetwork implements Network
{
    // Size of the network.
    private int n;

    // Array of vertex ids.
    private int[] vertices;

    // Build an instance of CompleteNetwork given its size.
    public CompleteNetwork(int n)
    {
	this.n = n;
	vertices = new int[n];
	for (int i = 0; i < n; i++) {
	    vertices[i] = i;
	}
    }

    // Return the size of the network.
    public int size()
    {
	return n;
    }

    // Return the id of a random vertex in the network.
    public int getRandomVertex()
    {
	return vertices[Random.uniform(n)];
    }

    // Return ids of the neighbors of vertex i, or null.
    public int[] getNeighbors(int i)
    {
	return vertices;
    }

    // Return the id of a random neighbor of vertex i, or -1.
    public int getRandomNeighbor(int i)
    {
	return getRandomVertex();
    }
}

// An abstraction for a network stored in GraphML format.
class GraphMLNetwork implements Network
{
    // Size of the network.
    private int n;

    // Adjacency list of the underlying network.
    private int adjList [][];

    // Abstraction of a vertex in a network.
    private class Vertex
    {
	// Vertex's auto-generated id.
	private int id;

	// Vertex's id as in the input GraphML file.
	private int realId;

	// Build an instance of Vertex given its id.
	public Vertex(int id)
	{
	    this.id = id;
	}

	// Set the real id of the vertex to the specified value.
	public void setRealId(int realId)
	{
	    this.realId = realId;
	}

	// Return the real id of the vertex.
	public int getId()
	{
	    return realId;
	}
    }

    // Abstraction of an edge in a network.
    private class Edge
    {
	// Edge's id.
	private int id;

	// Build an instance of Edge given its id.
	public Edge(int id)
	{
	    this.id = id;
	}

	// Return the id of the edge.
	public int getId()
	{
	    return id;
	}
    }

    // A factory class for building instances of Vertex.
    private class VertexFactory implements Factory
    {
	// Vertex id counter.
	private int id = 0;

	// Return an instance of Vertex.
	public Vertex create()
	{
	    return (new Vertex(id++));
	}
    }

    // A factory class for building instances of Edge.
    private class EdgeFactory implements Factory
    {
	// Edge id counter.
	private int id = 0;

	// Return an instance of Edge.
	public Edge create()
	{
	    return (new Edge(id++));
	}
    }

    // Build an instance of GraphMLNetwork given the name of a GraphML file.
    public GraphMLNetwork(String fileName)
    {
	try {
	    GraphMLReader<UndirectedGraph<Vertex, Edge>, Vertex, Edge> gmlr =
		new GraphMLReader<UndirectedGraph<Vertex, 
	    Edge>, Vertex, Edge>(new VertexFactory(), new EdgeFactory());
	    UndirectedGraph<Vertex, Edge> graph = 
		new UndirectedSparseMultigraph<Vertex, Edge>();
	    gmlr.load(fileName, graph);
	    BidiMap<Vertex, String> vertexIds = gmlr.getVertexIDs();
	    for (Vertex n : graph.getVertices()) {
		n.setRealId(Integer.parseInt(vertexIds.get(n)));
	    }
	    n = graph.getVertexCount();
	    adjList = new int[n][];
	    for (Vertex n : graph.getVertices()) {
		if (adjList[n.getId()] == null) {
		    adjList[n.getId()] = 
			new int[graph.getNeighborCount(n)];
		}
		int i = 0;
		for (Vertex m : graph.getNeighbors(n)) {
		    adjList[n.getId()][i++] = m.getId();
		}
	    }
        }
	catch (ParserConfigurationException e) {
	    System.err.println(e.getMessage());
	}
	catch (SAXException e) {
	    System.err.println(e.getMessage());
	}
	catch (IOException e) {
	    System.err.println(e.getMessage());
	}
    }

    // Return the size of the network.
    public int size()
    {
	return n;
    }

    // Return the id of a random vertex in the network.
    public int getRandomVertex()
    {
	return Random.uniform(n);
    }

    // Return ids of the neighbors of vertex i, or null.
    public int[] getNeighbors(int i)
    {
	return adjList[i];
    }

    // Return the id of a random neighbor of vertex i, or -1.
    public int getRandomNeighbor(int i)
    {
	int[] neighbors = getNeighbors(i);
	if (neighbors != null && neighbors.length > 0) {
	    return neighbors[Random.uniform(neighbors.length)];
	}
	return -1;
    }
}