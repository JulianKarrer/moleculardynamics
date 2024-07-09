# type: ignore

from lib import *
import numpy as np
from scipy.optimize import curve_fit


def model_latent_heat(x, slope1: float, slope2: float, temp_melt: float, energy_melt: float, latent_heat: float) -> float:
    """
    Model of the expected temperature of gold for a given energy `x`. 
    - `slope1` and `slope2` are the slopes of the function in (K/(eV/N)) left and right of the plateau (specific heat capacity)
    - `temp_melt` is the melting temperature, after which the function is constant for a while
    - `energy_melt` is the energy per atom in (eV/N) at which the melting point is reached
    - `latent_heat` is the heat per atom in (eV/N) for which the function is constant, used for the phase transition
    the function rises with `slope` until reaching the point (`energy_melt`|`temp_melt`), 
    then stays constant for `latent_heat` before rising with `slope` again
    """
    res = np.zeros(len(x))
    # create masks
    low = x < energy_melt
    mid = (energy_melt < x) & (x < energy_melt+latent_heat)
    high = x > energy_melt+latent_heat
    # calculate each area of the funciton
    res[low] = temp_melt + slope1*(x[low]-energy_melt)
    res[mid] = temp_melt
    res[high] = temp_melt + slope2*(x[high]-(energy_melt+latent_heat))
    # return the result
    return res


def plot_gold_temp_over_energy_fitted():
    data: CsvDict = data_from_csv(
        "gold_t_e.csv", "n", "e_total", "temperature")
    info = r"Gupta/Ducastelle EAM for Au, Mackay Icosahedra, $\Delta t$=1fs, $\tau_{relax}=1000\Delta t$, $\Delta Q=10K \cdot \frac{3}{2}N k_B$"

    # plot models of the temperature and average total energy
    fig, ax = std_plot(
        r"Fitted Model of Temperature $T$ over Average Total Energy $\frac{E_{total}}{N}$ for Different Numbers of Atoms $N$ in the Cluster",
        r"Average total Energy $\frac{E_{total}}{N}$ (eV)",
        r"Temperature $T$ (K)",
        info,
    )
    params = dict()
    for name, [xs, ys] in data.items():
        # normalize x by number of atoms for plotting
        xs = [x/float(name) for x in xs]
        N = float(name)
        # use non-linear least squares to model to the data
        x = np.array(xs, dtype=np.float64)
        y = np.array(ys, dtype=np.float64)
        popt, _pcov = curve_fit(model_latent_heat, x, y,
                                p0=np.array([3000, 3000, 700, -3.42, 0.05]))
        # plot the result
        plot_line(ax, x, model_latent_heat(x, *(popt.tolist())), r"N = "+name)
        ax.scatter(
            xs, ys, color=ax.lines[-1].get_color(), marker=next(markerstyles), s=2,)
        # save the parameters used, plot the fitted model
        params[name] = popt.tolist()
    ax.legend(fontsize="12")

    # plot model of heat capacity and cluster size
    fig_c, ax_c = std_plot(
        r"Heat Capacity $C$ over the Number of Atoms $N$ in an Au Mackay Isocahedron",
        r"Numbers of Atoms $N$",
        r"Heat Capacity $C$ (ev/K)",
        info,
    )
    cluster_sizes = [float(name) for name in data.keys()]
    ax_c.plot(cluster_sizes, [
              float(name)/((params[name][0]+params[name][1])/2.0) for name in data.keys()])

    # plot model of melting point and cluster size
    fig_melt, ax_melt = std_plot(
        r"Melting Point $T_{melt}$ over the Number of Atoms $N$ in an Au Mackay Isocahedron",
        r"Numbers of Atoms $N$",
        r"Melting Point $T_{melt}$ (K)",
        info,
    )
    ax_melt.plot(cluster_sizes, [params[name][2] for name in data.keys()])

    # plot model of latent heat and cluster size
    fig_latent, ax_latent = std_plot(
        r"Latent Heat $\Delta Q_{latent}$ over the Number of Atoms $N$ in an Au Mackay Isocahedron",
        r"Numbers of Atoms $N$",
        r"Latent Heat $\Delta Q_{lat}$ (eV)",
        info,
    )
    ax_latent.plot(cluster_sizes, [params[name][4] for name in data.keys()])

    # save figures ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fig.savefig("fit_mackay_t_e.png", dpi=400)
    fig_c.savefig("fit_mackay_c_n.png", dpi=400)
    fig_melt.savefig("fit_mackay_t_melt_n.png", dpi=400)
    fig_latent.savefig("fit_mackay_ql_n.png", dpi=400)


def plot_gold_temp_over_energy():
    data: CsvDict = data_from_csv(
        "../builddir/gold_t_e.csv", "n", "e_total", "temperature")
    info = r"Gupta/Ducastelle EAM for Au, Mackay Icosahedra, $\Delta t$=1fs, $\tau_{relax}=1000\Delta t$, $\Delta Q=10K \cdot \frac{3}{2}N k_B$"
    fig, ax = std_plot(
        r"Temperature $T$ over Average Total Energy $\frac{E_{total}}{N}$ for Different Numbers of Atoms $N$ in the Cluster",
        r"Average total Energy $\frac{E_{total}}{N}$ (eV)",
        r"Temperature $T$ (K)",
        info,
    )
    # SETTINGS
    discard_first_n = 0
    # temperature deviation from linear model in kelvin that marks the melting point
    threshold = 10
    # number of subsequent samples beyond threshold to ensure melting point is reached
    n_over_tresh = 5

    # PLOT the temperature-energy curves and determine the heat capacity
    total_x_min = 0
    C_x_max = -3.46
    total_x_max = -9999
    C_graph = []
    colours = dict()
    for name, [xs, ys] in data.items():
        N = float(name)
        xs = [x/N for x in xs][discard_first_n:]
        ys = ys[discard_first_n:]
        total_x_min = min(total_x_min, min(xs))
        total_x_max = max(total_x_max, max(xs))
        plot_line(ax, xs, ys, r"N = "+name)
        colours[name] = ax.lines[-1].get_color()
        # determine heat capacity
        lin_area = [(x, y) for (x, y) in zip(xs, ys) if x < C_x_max]
        C = ((lin_area[-1][0] - lin_area[0][0]) *
             N / (lin_area[-1][1] - lin_area[0][1]))
        C_graph += [(N, C)]
    ax.axvspan(total_x_min, C_x_max, alpha=0.1)
    ax.legend(fontsize="12")

    # PLOT heat capacity as a function of cluster size
    fig_c, ax_c = std_plot(
        r"Heat Capacity $C$ over the Number of Atoms $N$ in an Au Mackay Isocahedron",
        r"Numbers of Atoms $N$",
        r"Heat Capacity $C$ (ev/K)",
        info,
    )
    xs_c = [x for (x, _y) in C_graph]
    ys_c = [y for (_x, y) in C_graph]
    ax_c.plot(xs_c, ys_c, "-", marker="o")

    # determine the melting points
    first = True
    melting_points: dict[tuple[float, float]] = dict()
    for name, [xs, ys] in data.items():
        # find a function that extrapolates the heat curve in the shaded area linearly
        xs = [x/float(name) for x in xs][discard_first_n:]
        ys = ys[discard_first_n:]
        (n, c) = [(n, c) for (n, c) in C_graph if n == float(name)][0]
        slope = 1.0/(c/n)
        marker_x, marker_y = [(x, y)
                              for (x, y) in zip(xs, ys) if x < C_x_max][-1]
        end_x = (1000-marker_y)/slope+marker_x
        def f_extrap(x) -> float: return (x-marker_x)*slope+marker_y
        # draw the extrapolation for the first plot
        if first:
            ax.axline((marker_x, marker_y), (end_x,
                                             1000), alpha=0.2, color=colours[name])
            first = False

        # find the first point that is further than some threshold from the extrapolated function
        b_x, b_y = [(x, y) for i, (x, y) in enumerate(zip(xs, ys))
                    if (
                        x > C_x_max
                        and i < len(xs)-(n_over_tresh+1)
                        and all([abs(ys[i+j]-f_extrap(xs[i+j])) > threshold for j in range(n_over_tresh+1)])
        )][0]
        melting_points[name] = (b_x, b_y)
        ax.plot(b_x, b_y, marker="o", color=colours[name])

    # determine latent heat
    lat_x_min = -3.365
    lat_heat_points: dict[tuple[float, float]] = dict()
    for name, [xs, ys] in data.items():
        # find a function that extrapolates the heat curve in the shaded area linearly
        xs = [x/float(name) for x in xs][discard_first_n:]
        ys = ys[discard_first_n:]
        lower_x, lower_y = [(x, y)
                            for (x, y) in zip(xs, ys) if x > lat_x_min][0]
        higher_x, higher_y = xs[-1], ys[-1]
        slope = (higher_y-lower_y)/(higher_x - lower_x)
        def f_extrap(x) -> float: return (x-lower_x)*slope+lower_y
        b_x, b_y = [(x, y) for i, (x, y) in enumerate(zip(xs, ys))
                    if (
                        x < lat_x_min
                        and i > (n_over_tresh+1)
                        and all([abs(ys[i-j]-f_extrap(xs[i-j])) > threshold for j in range(n_over_tresh+1)])
        )][-1]
        lat_heat_points[name] = (b_x, b_y)
        ax.plot(b_x, b_y, marker="o", color=colours[name])
    ax.axvspan(lat_x_min, total_x_max, alpha=0.1, color="green")

    # PLOT melting point as funciton of number of atoms
    fig_b, ax_b = std_plot(
        r"Melting Point $T_{melt}$ over the Number of Atoms $N$ in an Au Mackay Isocahedron",
        r"Numbers of Atoms $N$",
        r"Melting Point $T_{melt}$ (K)",
        info,
    )
    xs_b = []
    ys_b = []
    for name, [xs, ys] in data.items():
        xs_b += [float(name)]
        ys_b += [(lat_heat_points[name][1] -
                 melting_points[name][1])/2+melting_points[name][1]]
    ax_b.plot(xs_b, ys_b, "-", marker="o")

    # PLOT latent heat as funciton of number of atoms
    fig_l, ax_l = std_plot(
        r"Latent Heat $\Delta Q_{latent}$ over the Number of Atoms $N$ in an Au Mackay Isocahedron",
        r"Numbers of Atoms $N$",
        r"Latent Heat $\Delta Q_{lat}$ (eV)",
        info,
    )
    xs_l = []
    ys_l = []
    for name, [xs, ys] in data.items():
        xs_l += [float(name)]
        ys_l += [(lat_heat_points[name][0] - melting_points[name][0])]
    ax_l.plot(xs_l, ys_l, "-", marker="o")

    # save the figures
    fig.savefig("mackay_t_e.png", dpi=400)
    fig_c.savefig("mackay_c_n.png", dpi=400)
    fig_b.savefig("mackay_t_melt_n.png", dpi=400)
    fig_l.savefig("mackay_ql_n.png", dpi=400)


def plots_optimal_eam_timestep():
    # what runs to plot for each figure
    times_plot = ["1.00", "5.00", "10.00", "20.00"]
    times_diff_plot = ["1.00", "5.00", "10.00",]

    fig, ax = std_plot(
        r"Time Evolution of the Hamiltonian for different $\Delta t$",
        r"Time $t$ (fs)",
        r"$\frac{d}{dt} E_{total}$ (eV/fs)",
        r"Gupta/Ducastelle EAM, N=923 Isocahedron, Au at 500K"
    )
    fig_d, ax_d = std_plot(
        r"Time Derivative (Central Differences) of the Hamiltonian for different $\Delta t$",
        r"Time $t$ (fs)",
        r"$E_{total}$ (eV)",
        r"Gupta/Ducastelle EAM, N=923 Isocahedron, Au at 500K"
    )
    data: CsvDict = data_from_csv(
        "dt_eam.csv", "dt", "t", "Hamiltonian")
    for name, [xs, ys] in data.items():
        if name in times_diff_plot:
            x_d, y_d = d_dx(xs, ys)
            plot_line(ax_d, x_d[5:], y_d[5:], r"$\Delta t$ = "+name+"fs")
        if name in times_plot:
            plot_line(ax, xs, ys, r"$\Delta t$ = "+name+"fs")
    ax_d.legend(fontsize="12")
    fig_d.savefig("dt_diff_eam.png", dpi=400)
    fig.savefig("dt_eam.png", dpi=400)


if __name__ == "__main__":
    plots_optimal_eam_timestep()
    plot_gold_temp_over_energy()
    plot_gold_temp_over_energy_fitted()
