// Payoff1.java

// Defines the payoff function for the continuous prisoner's 
// dilemma game with linear benefit and cost functions.

// b1 = params[0]; c1 = params[1]
class Payoff1 extends Payoff
{
    // Benefit function; b1 * X
    public double B(double x)
    {
	return params[0] * x;
    }

    // Cost function; c1 * x
    public double C(double x)
    {
	return params[1] * x;
    }

    // Payoff function.
    public double payoff(double x, double y)
    {
	return B(y) - C(x);
    }
}