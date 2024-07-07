# type: ignore

from math import atan2, degrees
from itertools import cycle
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 15})


# read data from csv and plot it
def plot_from_csv(filename: str):
    xs: list[float] = []
    ys: list[float] = []
    with open(filename) as run_rs:
        for line in run_rs.readlines()[1:]:
            split = line.split(",")
            xs += [float(split[1])]  # e_total
            ys += [float(split[2])]  # temperature
    color = next(plt.gca()._get_lines.prop_cycler)['color']
    plt.plot(xs, ys, next(linestyles), color=color)


# cycle through styles to improve readability
linestyles = cycle(["-", "--", "-.", ":"])
markerstyles = cycle(["o", "v", "s", "D", "X", "^", "*"])

# plot the data
plt.title(r"Temperature $T$ over Total Energy $E_{total}$")
plt.figtext(0.99, 0.01, r"Gupta/Ducastelle EAM, N=923 Isocahedron, Au",
            horizontalalignment='right')
plt.xlabel(r"Hamiltonian $E_{total}$ (eV)")
plt.ylabel(r"Temperature $T$ (K)")
plt.tight_layout()
plt.gcf().set_size_inches(12., 8., forward=True)
plot_from_csv("../builddir/gold_t_e.csv")
plt.savefig("gold_t_e.png")
plt.clf()
