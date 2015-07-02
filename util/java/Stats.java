// Stats.java

// A utility class for computing basic statistics.
class Stats
{
    // Singleton class, so disallow instantiation.
    private Stats() { }

    // Return maximum value in array a, or -infinity.
    public static double max(double[] a) 
    {
	return max(a, 0, a.length - 1);
    }

    // Return maximum value in subarray a[lo..hi], or -infinity.
    public static double max(double[] a, int lo, int hi) 
    {
        double max = Double.NEGATIVE_INFINITY;
        for (int i = lo; i <= hi; i++) {
            if (a[i] > max) max = a[i];
        }
        return max;
    }
    
    // Return minimum value in array a, or infinity.
    public static double min(double[] a) 
    {
        return min(a, 0, a.length - 1);
    }

    // Return minimum value in subarray a[lo..hi], or infinity.
    public static double min(double[] a, int lo, int hi) 
    {
        double min = Double.POSITIVE_INFINITY;
        for (int i = lo; i <= hi; i++) {
            if (a[i] < min) min = a[i];
        }
        return min;
    }

    // Return average value in array a, or NaN.
    public static double mean(double[] a) 
    {
	return mean(a, 0, a.length - 1);
    }

    // Return average value in subarray a[lo..hi], or NaN.
    public static double mean(double[] a, int lo, int hi) 
    {
        if (a.length == 0) return Double.NaN;
        double sum = 0.0;
        for (int i = 0; i < a.length; i++) {
            sum = sum + a[i];
        }
        return sum / a.length;
    }

    // Return sample variance of array a, or NaN.
    public static double var(double[] a) 
    {
	return var(a, 0, a.length - 1);
    }

    // Return sample variance of subarray a[lo..hi], or NaN.
    public static double var(double[] a, int lo, int hi) 
    {
        int length = hi - lo + 1;
        if (length == 0) return Double.NaN;
        double avg = mean(a, lo, hi);
        double sum = 0.0;
        for (int i = lo; i <= hi; i++) {
            sum += (a[i] - avg) * (a[i] - avg);
        }
        return sum / (length - 1);
    }

    // Return population variance of array a, or NaN.
    public static double varp(double[] a) 
    {
	return varp(a, 0, a.length - 1);
    }

    // Return population variance of subarray a[lo..hi],  or NaN.
    public static double varp(double[] a, int lo, int hi) 
    {
        int length = hi - lo + 1;
        if (length == 0) return Double.NaN;
        double avg = mean(a, lo, hi);
        double sum = 0.0;
        for (int i = lo; i <= hi; i++) {
            sum += (a[i] - avg) * (a[i] - avg);
        }
        return sum / length;
    }

    // Return sample standard deviation of array a, or NaN.
    public static double stddev(double[] a) 
    {
        return stddev(a, 0, a.length - 1);
    }

    // Return sample standard deviation of subarray a[lo..hi], or NaN.
    public static double stddev(double[] a, int lo, int hi) 
    {
        return Math.sqrt(var(a, lo, hi));
    }

    // Return population standard deviation of array a, or NaN.
    public static double stddevp(double[] a) 
    {
        return stddevp(a, 0, a.length - 1);
    }

    // Return population standard deviation of subarray a[lo..hi], or NaN.
    public static double stddevp(double[] a, int lo, int hi) 
    {
        return Math.sqrt(varp(a, lo, hi));
    }

    // Return sum of all values in array a.
    public static double sum(double[] a) 
    {
        return sum(a, 0, a.length - 1);
    }

    // Return sum of all values in subarray a[lo..hi].
    public static double sum(double[] a, int lo, int hi) 
    {
        double sum = 0.0;
        for (int i = lo; i <= hi; i++) {
            sum += a[i];
        }
        return sum;
    }
}
