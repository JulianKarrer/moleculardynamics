# type: ignore

from math import atan2, degrees
from itertools import cycle
import numpy as np
import matplotlib.pyplot as plt
# plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 15})


def labelLine(line, x, label=None, align=True, **kwargs):
    # FUNCTION LABELLINE FROM:
    # https://stackoverflow.com/questions/16992038/how-to-place-inline-labels-in-a-line-plot
    # AS AVAILABLE AT:
    # https://github.com/cphyc/matplotlib-label-lines
    # UNDER THE MIT LICENSE

    ax = line.axes
    xdata = line.get_xdata()
    ydata = line.get_ydata()

    if (x < xdata[0]) or (x > xdata[-1]):
        print('x label location is outside data range!')
        return

    # Find corresponding y co-ordinate and angle of the line
    ip = 1
    for i in range(len(xdata)):
        if x < xdata[i]:
            ip = i
            break

    y = ydata[ip-1] + (ydata[ip]-ydata[ip-1]) * \
        (x-xdata[ip-1])/(xdata[ip]-xdata[ip-1])

    if not label:
        label = line.get_label()

    if align:
        # Compute the slope
        dx = xdata[ip] - xdata[ip-1]
        dy = ydata[ip] - ydata[ip-1]
        ang = degrees(atan2(dy, dx))

        # Transform to screen co-ordinates
        pt = np.array([x, y]).reshape((1, 2))
        trans_angle = ax.transData.transform_angles(np.array((ang,)), pt)[0]

    else:
        trans_angle = 0

    # Set a bunch of keyword arguments
    if 'color' not in kwargs:
        kwargs['color'] = line.get_color()

    if ('horizontalalignment' not in kwargs) and ('ha' not in kwargs):
        kwargs['ha'] = 'center'

    if ('verticalalignment' not in kwargs) and ('va' not in kwargs):
        kwargs['va'] = 'center'

    if 'backgroundcolor' not in kwargs:
        kwargs['backgroundcolor'] = ax.get_facecolor()

    if 'clip_on' not in kwargs:
        kwargs['clip_on'] = True

    if 'zorder' not in kwargs:
        kwargs['zorder'] = 2.5

    ax.text(x, y, label, rotation=trans_angle, **kwargs)


def labelLines(lines, align=True, xvals=None, **kwargs):

    ax = lines[0].axes
    labLines = []
    labels = []

    # Take only the lines which have labels other than the default ones
    for line in lines:
        label = line.get_label()
        if "_line" not in label:
            labLines.append(line)
            labels.append(label)

    if xvals is None:
        xmin, xmax = ax.get_xlim()
        xvals = np.linspace(xmin, xmax, len(labLines)+2)[1:-1]

    for line, x, label in zip(labLines, xvals, labels):
        labelLine(line, x, label, align, **kwargs)


# MY OWN FUNCTIONS FROM HERE ON


# read data from csv and plot it
def plot_from_csv(filename: str, firstcol: str, name: str):
    xs: list[float] = []
    ys: list[float] = []
    with open(filename) as run_rs:
        for line in run_rs.readlines()[1:]:
            split = line.split(",")
            if split[0] == firstcol:
                xs += [float(split[1])]
                ys += [float(split[2])]
    color = next(plt.gca()._get_lines.prop_cycler)['color']
    # 
    plt.plot(xs, ys, next(linestyles), label=name, color=color, marker=next(markerstyles))


# cycle through styles to improve readability
linestyles = cycle([".-", "--", "-.", ":"])
markerstyles = cycle(["o", "v", "s", "D", "X", "^", "*"])


# TODO TAU ALS ZEITEINHEIT SETZEN 

# plot the data
def plot_all_sets(legend: bool):
    # setup titles, legends etc
    plt.suptitle(r"Temperature over time for different \tau of a Berendsen thermostat")
    plt.title(r"dt = 10^{-3}, 100 atoms, \sigma = 0.8, LJTS at r_c = 5\sigma, T^{target}=300K")
    plt.xlabel(r"Time (tau)")
    plt.ylabel(r"Temperature (K)")
    plt.tight_layout()
    plt.gcf().set_size_inches(18., 9., forward=True)
    for i in reversed(list(range(1, 10+1))):
        relax_time = str(0.01*i)
        plot_from_csv("../builddir/thermostat_temps.csv", relax_time, relax_time)
    if legend:
        plt.legend(fontsize="12")
    else:
        labelLines(plt.gca().get_lines(), zorder=2.5, fontsize=8)
   
    plt.show()
    plt.savefig("temps_relax"
                + ("_legend" if legend else "")
                + ".png")
    plt.clf()


# plot all variations of the graph and save to file
for legend in [True, False]:
    plot_all_sets(legend)
