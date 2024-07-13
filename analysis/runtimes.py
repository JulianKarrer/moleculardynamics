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
def plot_from_csv(filename: str, firstcol: str, name: str, plot_stddev: bool):
    stddev_low: list[float] = []
    stddev_high: list[float] = []
    xs: list[float] = []
    ys: list[float] = []
    with open(filename) as run_rs:
        for line in run_rs.readlines()[1:]:
            split = line.split(",")
            if split[0] == firstcol:
                xs += [float(split[1])]
                ys += [float(split[2])]
                stddev_low += [float(split[2])-float(split[5])]
                stddev_high += [float(split[2])+float(split[5])]
    color = next(plt.gca()._get_lines.prop_cycler)['color']
    plt.plot(xs, ys, next(linestyles), label=name,
             marker=next(markerstyles), color=color)
    if plot_stddev:
        plt.fill_between(xs, stddev_low, stddev_high, color=color, alpha=0.2)


# cycle through styles to improve readability
linestyles = cycle([".-", "--", "-.", ":"])
markerstyles = cycle(["o", "v", "s", "D", "X", "^", "*"])


# plot the data
def plot_all_sets(log: bool, legend: bool, plot_stddev: bool):
    # setup titles, legends etc
    plt.suptitle("Average runtime across 10 simulations of 50 timesteps each for the LJDS (direct summation) \
                and LJTS (truncated and shifted) potential in C++, Rust and Rust+Rayon respectively")
    plt.xlabel(r"Number of Atoms")
    plt.ylabel(r"Average Runtime in Microseconds ($\mu s$)")
    # plt.tight_layout()
    plt.gcf().set_size_inches(18., 9., forward=True)
    plot_from_csv("../builddir/runtimes.csv",
                  "direct", "LJDS-CPP", plot_stddev)
    plot_from_csv("../builddir/runtimes.csv", "ljts", "LJTS-CPP", plot_stddev)
    plot_from_csv("../rust/runtimes.csv", "direct", "LJDS-RS", plot_stddev)
    plot_from_csv("../rust/runtimes.csv", "ljts", "LJDS-RS", plot_stddev)
    # plot_from_csv("../rust/runtimes_par.csv",
    #   "direct", "LJDS-RS+R", plot_stddev)
    # plot_from_csv("../rust/runtimes_ts.csv", "ljts", "LJTS-RS", plot_stddev)
    # plot_from_csv("../rust/runtimes_par.csv", "ljts", "LJTS-RS+R", plot_stddev)
    if legend:
        plt.legend(fontsize="12")
    else:
        labelLines(plt.gca().get_lines(), zorder=2.5, fontsize=8)
    if log:
        plt.yscale("log")
        plt.xscale("log")
    else:
        plt.yscale("linear")
    if plot_stddev:
        plt.figtext(0.5, 0.02, "(shaded regions show the corrected standard deviation)", ha="center",
                    fontsize=10)
    # plt.show()
    plt.savefig("runtimes"
                + ("_log" if log else "")
                + ("_legend" if legend else "")
                + ("_stddev" if plot_stddev else "")
                + ".png")
    plt.clf()


# plot all variations of the graph and save to file
for log in [True, False]:
    for legend in [True, False]:
        for plot_stddev in [True, False]:
            plot_all_sets(log, legend, plot_stddev)

# plot_from_csv("../builddir/runtimes_ts.csv", "ljts", "LJTS-CPP")
# plot_from_csv("../rust/runtimes_ts.csv", "ljts", "LJTS-RS")
# plot_from_csv("../rust/runtimes_ts_par.csv", "ljts", "LJTS-RS+R")
