# Plots the variation of basic reproductive ratio R0 with virulence c given 
# alpha, beta_0, gamma, and d.

from pylab import *
import math
import sys

def main(argv):
    """
    Entry point.
    """

    if len(argv) != 4:
        print "Usage: python r0_plot.py <a> <b> <recovery> <d>"
        sys.exit(0)
    a = float(argv[0])
    b = float(argv[1])
    recovery = float(argv[2])
    d = float(argv[3])
    x = []
    y = []
    c = 0.0
    c_star = c
    r_max = 0.0
    for i in range(0, 1000):
        r = b * math.pow(c, a) / (1 - (1 - d) * (1 - c) * (1 - recovery))
        if r > r_max:
            r_max = r
            c_star = c
        x.append(c)
        y.append(r)
        c = c + 0.001

    figure(1, figsize = (7, 4.5), dpi = 500)
    plot(x, y, "#000000", alpha = 0.6, linewidth = 2.0)
    axvline(x = c_star, linestyle = "-.", linewidth = 1.0, color = "black", 
            label = "$c^\star = %4.3f$" %(c_star))
    xlabel(r"$c$")
    ylabel(r"$R_0$")
    legend()
    savefig("r0.pdf", format = "pdf")
    close(1)

if __name__ == "__main__":
    main(sys.argv[1:])

