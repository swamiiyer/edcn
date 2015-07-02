// Payoff3.java
// 
// Defines the payoff function for the continuous snowdrift  
// game with quadratic benefit and cost functions.

// b1 = params[0]; b2 = params[1]; c1 = params[2]; c2 = params[3]
class Payoff3 extends Payoff
{
    // Benefit function; -b2 * x^2 + b1 * X
    public double B(double x)
    {
	return -params[1] * x * x + params[0] * x;
    }

    // Cost function; -c2 * x^2 + c1 * x
    public double C(double x)
    {
	return -params[3] * x * x + params[2] * x;
    }

    // Payoff function.
    public double payoff(double x, double y)
    {
	return B(x + y) - C(x);
    }
}