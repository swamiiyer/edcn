// Payoff4.java
// 
// Defines the payoff function for the continuous tragedy of the 
// commons game with cubic benefit and quadratic cost functions.

// b1 = params[0]; b2 = params[1]; b3 = params[2]; c1 = params[3]
class Payoff4 extends Payoff
{
    // Benefit function; -b3 * x^3 + b2 * x^2 + b1 * X
    public double B(double x)
    {
	return -params[2] * x * x * x + params[1] * x * x + params[0] * x;
    }

    // Cost function; c1 * x^2
    public double C(double x)
    {
	return params[3] * x * x;
    }

    // Payoff function.
    public double payoff(double x, double y)
    {
	return B(x) - C(x + y);
    }
}