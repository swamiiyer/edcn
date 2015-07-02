// Payoff2.java
// 
// Defines the payoff function for the continuous prisoner's 
// dilemma game with quadratic benefit and cost functions.

// b1 = params[0]; b2 = params[1]; c1 = params[2]
class Payoff2 extends Payoff
{
    // Benefit function; -b2 * x^2 + b1 * X
    public double B(double x)
    {
	return -params[1] * x * x + params[0] * x;
    }

    // Cost function; c1 * x^2
    public double C(double x)
    {
	return params[2] * x * x;
    }

    // Payoff function.
    public double payoff(double x, double y)
    {
	return B(y) - C(x);
    }
}