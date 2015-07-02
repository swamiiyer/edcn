// Random.java

// A utility class for working with random numbers.
class Random
{
    // A singleton instance of the underlying random number generator (RNG).
    private static java.util.Random random;  

    // Seed for the RNG.
    private static long seed; 

    // Set the seed for the RNG to the current time and create the singleton 
    // RNG instance.
    static 
    {
        seed = System.currentTimeMillis();
        random = new java.util.Random(seed);
    }
    
    // Singleton class, so disallow instantiation.
    private Random() {}

    // Return real number uniformly in [0, 1).
    public static double uniform() 
    {
        return random.nextDouble();
    }

    // Return int uniformly in [0, N).
    public static int uniform(int N) 
    {
        return random.nextInt(N);
    }

    // Return int uniformly in [a, b).
    public static int uniform(int a, int b) 
    {
        return a + uniform(b - a);
    }

    // Return real number uniformly in [a, b).
    public static double uniform(double a, double b) 
    {
        return a + uniform() * (b-a);
    }

    // Return a boolean, which is true with probability p, and false otherwise.
    public static boolean bernoulli(double p) 
    {
        return uniform() < p;
    }

    // Return a boolean, which is true with probability .5, and false otherwise.
    public static boolean bernoulli() 
    {
        return bernoulli(0.5);
    }

    // Return a real number with a standard Gaussian distribution.
    public static double gaussian() 
    {
        double r, x, y;
        do {
            x = uniform(-1.0, 1.0);
            y = uniform(-1.0, 1.0);
            r = x*x + y*y;
        } while (r >= 1 || r == 0);
        return x * Math.sqrt(-2 * Math.log(r) / r);
    }

    // Return a real number from a gaussian distribution with given mean and 
    public static double gaussian(double mean, double stddev) 
    {
        return mean + stddev * gaussian();
    }

    // Return a number from a discrete distribution: i with probability a[i].
    public static int discrete(double[] a) 
    {
        double r = uniform();
        double sum = 0.0;
        for (int i = 0; i < a.length; i++) {
            sum = sum + a[i];
            if (sum >= r) return i;
        }
        assert false;
        return -1;
    }
}