# type: ignore
from itertools import cycle
import matplotlib.pyplot as plt

type Graph = tuple[list[float], list[float]]
type CsvDict = dict[str,  Graph]


# read data from csv and return it


def data_from_csv(filename: str, run_identifier: str, name_x: str, name_y: str) -> dict[str,  tuple[list[float], list[float]]]:
    """From a given `filename`, parse a csv file containing data 
    for different runs of a measurement by inspecting the column names in the first row of the file:

    - `run_identifier` is the name of the column that is different between runs
    - `name_x` is the name of the column used for x-values
    - `name_y` is the name of the column used for y-values
    - returns a dictionary mapping the `run_identifier` value of each run 
    to the respective x and y values (given as a 2-tuple of lists of floats)
    """
    with open(filename) as run_rs:
        lines: list[str] = run_rs.readlines()
        firstline: str = lines[0]
        res: dict[str,  tuple[list[float], list[float]]] = dict()
        cat_index = [i for (i, x) in enumerate(
            firstline.split(",")) if x == run_identifier][0]
        x_index = [i for (i, x) in enumerate(
            firstline.split(",")) if x == name_x][0]
        y_index = [i for (i, y) in enumerate(
            firstline.split(",")) if y == name_y][0]
        for line in lines[1:]:
            split = line.split(",")
            res.setdefault(split[cat_index], ([], []))
            curxs, curys = res[split[cat_index]]
            res[split[cat_index]] = (
                curxs+[float(split[x_index])], curys+[float(split[y_index])])
        return res


def d_dx(xs: list[float], ys: list[float]) -> Graph:
    """
    Given two lists of length `n` conatining x and y values, differentiate y
    with respect to x numerically using a central difference scheme with O(dx^2) error.

    Returns lists of `n-2 ` length since the first and last elements have no neighbours
    """
    return xs[1:-1], [(ys[i+1]-ys[i-1])/(xs[i+1]-xs[i-1]) for i in range(1, len(ys)-1)]


# cycle through styles to improve readability
linestyles = cycle(["-", "--", "-.", ":"])
markerstyles = cycle(["o", "v", "s", "D", "X", "^", "*"])


def std_plot(title: str, xlabel: str, ylabel: str, info: str):
    plt.rcParams['text.usetex'] = True
    plt.rcParams.update({'font.size': 15})
    fig, ax = plt.subplots()
    fig.suptitle(title)
    fig.text(0.99, 0.01, info, horizontalalignment='right', fontsize="10")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.set_size_inches(16., 9., forward=True)
    return fig, ax


def plot_line(ax, xs, ys, label: str):
    colour = next(ax._get_lines.prop_cycler)['color']
    ax.plot(xs, ys, next(linestyles), label=label, color=colour)


# def w(d: float, h_bar: float) -> float:
#     """Compute the 1D normalized M4 Schoenberg B-Spline or Cubic Spline Kernel"""
#     q = d/h_bar
#     return (2/3)*(max(0., 2.0-q)**3 - 4*max(0., 1.0-q)**3)


# def smoothen(xs, ys):
#     """Smooth a graph by applying a Gaussian convolution (SPH filter) with a kernel support of two times the average spacing between x values"""
#     h_avg = sum([xs[i+1]-xs[i] for i in range(len(xs)-1)])/(len(xs)-1)
#     h_bar = 2*h_avg
#     return [sum([ys[j] * h_avg * w(abs(xs[i]-xs[j]), h_bar) for j in range(len(xs))]) for i in range(len(xs))]
