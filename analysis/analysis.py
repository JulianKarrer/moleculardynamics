# type: ignore

from lib import *
import numpy as np
from scipy.optimize import curve_fit
import subprocess


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
    info = r"Ducastelle EAM for Au, Cleri \& Rosato Parameters, Mackay Icosahedra, $\Delta t$=1fs, $\tau_{relax}=1000\Delta t$, $\Delta Q=10K \cdot \frac{3}{2}N k_B$"

    # plot models of the temperature and average total energy
    fig, ax = std_plot(
        r"Fitted Model of Temperature $T$ over Average Total Energy $\frac{E_{total}}{N}$ for Different Numbers of Atoms $N$ in the Cluster",
        r"Average total Energy $\frac{\mathcal{H}}{N}$ (eV)",
        r"Temperature $T$ (K)",
        info,
        1.5
    )
    params = dict()
    for name, [xs, ys] in data.items():
        # normalize x by number of atoms for plotting
        N = float(name)
        xs = [x/N for x in xs]
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
        r"Heat Capacity $C$ over the Number of Atoms $N$ in an Au Mackay Icosahedron",
        r"Numbers of Atoms $N$",
        r"Heat Capacity $C$ (ev/K)",
        info,
    )
    xs = [float(name) for name in data.keys()]
    ys_1 = [1.0/(params[name][0]/float(name)) for name in data.keys()]
    ys_2 = [1.0/(params[name][1]/float(name)) for name in data.keys()]
    plot_line(ax_c, xs, ys_1, 'solid', use_marker=True)
    plot_line(ax_c, xs, ys_2, 'fluid', use_marker=True)
    ax_c.legend()

    # plot model of melting point and cluster size
    fig_melt, ax_melt = std_plot(
        r"Melting Point $T_{melt}$ over the Number of Atoms $N$ in an Au Mackay Icosahedron",
        r"Numbers of Atoms $N$",
        r"Melting Point $T_{melt}$ (K)",
        info,
    )
    ax_melt.plot(xs, [params[name][2] for name in data.keys()])

    # plot model of latent heat and cluster size
    fig_latent, ax_latent = std_plot(
        r"Latent Heat $Q_{latent}$ over the Number of Atoms $N$ in an Au Mackay Icosahedron",
        r"Numbers of Atoms $N$",
        r"Latent Heat $Q_{lat}$ (eV)",
        info,
    )
    ax_latent.plot(xs, [params[name][4] for name in data.keys()])

    # save figures ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fig.savefig("fit_mackay_t_e.png", dpi=400)
    fig_c.savefig("fit_mackay_c_n.png", dpi=400)
    fig_melt.savefig("fit_mackay_t_melt_n.png", dpi=400)
    fig_latent.savefig("fit_mackay_ql_n.png", dpi=400)


def plot_gold_temp_over_energy():
    data: CsvDict = data_from_csv(
        "gold_t_e.csv", "n", "e_total", "temperature")
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
        r"Heat Capacity $C$ over the Number of Atoms $N$ in an Au Mackay Icosahedron",
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
        r"Melting Point $T_{melt}$ over the Number of Atoms $N$ in an Au Mackay Icosahedron",
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
        r"Latent Heat $\Delta Q_{latent}$ over the Number of Atoms $N$ in an Au Mackay Icosahedron",
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
    info = r"Ducastelle EAM, Cleri \& Rosato Parameters, N=923 Mackay Icosahedron, Au at 500K"

    fig, ax = std_plot(
        r"",  # r"Time Evolution of the Hamiltonian for different $\Delta t$",
        r"Time $t$ (fs)",
        r"$E_{total}$ (eV)",
        info
    )
    fig_d, ax_d = std_plot(
        r"",  # r"Time Evolution of the Derivative of the Hamiltonian for different $\Delta t$",
        r"Time $t$ (fs)",
        r"$\frac{d}{dt} E_{total}$ (eV/fs)",
        "derivative taken with central differences, " + info
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
    ax.legend(fontsize="12")
    fig_d.savefig("dt_diff_eam.png", dpi=400)
    fig.savefig("dt_eam.png", dpi=400)


def plots_optimal_ljds_timestep():
    # what runs to plot for each figure
    times_plot = ["0.001", "0.002", "0.005", "0.010"]
    times_diff_plot = ["0.001", "0.002", "0.005", "0.010"]
    info = r"Lennard-Jones direct summation, $\epsilon=1, \sigma=1, N=1000$, initially zero velocities, $r_0$ grid spacing"

    fig, ax = std_plot(
        r"",  # r"Time Evolution of the Hamiltonian for different $\Delta t$",
        r"Time $t$",
        r"$E_{total}$",
        info
    )
    fig_d, ax_d = std_plot(
        r"",  # r"Time Evolution of the Derivative of the Hamiltonian for different $\Delta t$",
        r"Time $t$",
        r"$\frac{d}{dt} E_{total}$ ",
        "derivative taken with central differences, " + info
    )
    data: CsvDict = data_from_csv(
        "dt_ljds.csv", "dt", "t", "Hamiltonian")
    for name, [xs, ys] in data.items():
        if name in times_diff_plot:
            x_d, y_d = d_dx(xs, ys)
            plot_line(ax_d, x_d[5:], y_d[5:], r"$\Delta t$ = "+name)
        if name in times_plot:
            plot_line(ax, xs, ys, r"$\Delta t$ = "+name)
    ax_d.legend(fontsize="12")
    ax.legend(fontsize="12")
    fig_d.savefig("dt_diff_ljds.png", dpi=400)
    fig.savefig("dt_ljds.png", dpi=400)


def plots_runtime_ljds_ljts():
    info = r"regular grid with spacing $r_0=2^{\frac{1}{6}}\sigma$, every other $xy$-plane offset by $0.5r_0$, $T_0=100K, \Delta t=0.001, \sigma=1.44, r_c=5\sigma$, 100 time steps per run"

    data: CsvDict = data_from_csv(
        "../builddir/runtimes.csv", "direct summation or ljts", "number of atoms", "average runtime")
    fig, ax = std_plot(
        r"",  # r"Time Evolution of the Hamiltonian for different $\Delta t$",
        r"Number of Atoms",
        r"Average Runtime in Microseconds across 10 Runs ($\mu s$)",
        info
    )
    for name, [xs, ys] in data.items():
        if name == "direct":
            plot_line(ax, xs, ys, r"Lennard-Jones Direct Summation")
        elif name == "ljts":
            plot_line(ax, xs, ys, r"Lennard-Jones Truncated and Shifted")
    ax.legend(fontsize="12")
    fig.savefig("runtimes_ljds_ljts.png", dpi=400)


def plot_parallel_energy_conserved():
    data = data_from_csv("hamiltonian_par.csv",
                         "nb_processes", "n", "e_total")
    fig, ax = std_plot(
        r"Time Evolution of the Hamiltonian for different number of processes",
        r"Time $t$ (fs)",
        r"$E_{total}$ (eV)",
        r"Same scenario each time: Au Mackay Icosahedron N=923 with initial zero velocities"
    )
    ys_prev = []
    for name, (xs, ys) in data.items():
        xs = [i for i in range(len(xs))]
        plot_line(ax, xs, ys, name+" processes")
    ax.legend(fontsize="12")
    fig.savefig("par_energy_consv.png", dpi=400)


def stress_strain():
    # N=3050, 300K, 10**8 strainrate
    fig, ax = std_plot(
        r"Force $F_{zz}$ acting on the $xy$-Plane as a Function of Strain $\frac{\Delta V_z}{V_{z_0}}$",
        r"Strain $\frac{\Delta V_z}{V_{z_0}}$",
        r"Force $F_{zz}=\frac{V\cdot\sigma_{zz}}{V_z}$ $\left(\frac{eV}{Å}\right)$",
        r"N=3050, Ducastelle EAM Au whisker aligned with z-axis, periodic boundary in z-direction, $\Delta t=5fs, \tau=1ps, \epsilon=10^8, T=300K$"
    )
    for v0, (xs, ys) in data_from_csv("whisker_3050/whisker_small.csv",
                                      "v_z_0", "v_z", "f_zz").items():
        ax.plot([x/float(v0)-1 for x in xs], ys)
    fig.tight_layout()
    fig.savefig("stress-strain-3050.png", dpi=400)

    # N=51_500, 300K, 10**8 strainrate
    fig, ax = std_plot(
        r"Force $F_{zz}$ acting on the $xy$-Plane as a Function of Strain $\frac{\Delta V_z}{V_{z_0}}$",
        r"Strain $\frac{\Delta V_z}{V_{z_0}}$",
        r"Force $F_{zz}=\frac{V\cdot\sigma_{zz}}{V_z}$ $\left(\frac{eV}{Å}\right)$",
        r"Ducastelle EAM Au whisker aligned with z-axis, periodic boundary in z-direction, $\Delta t=5fs, \tau=1ps$",
        1.7
    )

    # small one
    for v0, (xs, ys) in data_from_csv("whisker_3050/whisker_small.csv",
                                      "v_z_0", "v_z", "f_zz").items():
        plot_line(ax, [x/float(v0)-1 for x in xs],
                  ys, r"$N=3050, \epsilon=\cdot 10^8, T=300K$")

    for v0, (xs, ys) in data_from_csv("whisker_51500/whisker_large_fast.csv",
                                      "v_z_0", "v_z", "f_zz").items():
        plot_line(ax, [x/float(v0)-1 for x in xs],
                  ys, r"$N=51500,\epsilon=2\cdot 10^8, T=300K$")

    for v0, (xs, ys) in data_from_csv("whisker_51500/whisker_large_cold.csv",
                                      "v_z_0", "v_z", "f_zz").items():
        plot_line(ax, [x/float(v0)-1 for x in xs],
                  ys, r"$N=51500,\epsilon=1\cdot 10^8, T=200K$")
    for v0, (xs, ys) in data_from_csv("whisker_51500/whisker_large.csv",
                                      "v_z_0", "v_z", "f_zz").items():
        plot_line(ax, [x/float(v0)-1 for x in xs],
                  ys, r"$N=51500,\epsilon=1\cdot 10^8, T=300K$")
    fig.tight_layout()
    ax.legend(fontsize="12")
    fig.savefig("stress-strain-51_500-strainrate.png", dpi=400)


def run(command: str):
    subprocess.run(command.split(" "))


def stress_strain_animation():
    run("rm ../videos/whisker_51500_graph/*.png")
    data = data_from_csv("whisker_51500/whisker_large.csv",
                         "v_z_0", "v_z", "f_zz")
    plt.rcParams['text.usetex'] = True
    plt.rcParams.update({'font.size': 20})
    xlabel = r"Strain $\frac{\Delta V_z}{V_{z_0}}$"
    ylabel = r"Force $F_{zz}$"
    info = r"N=51500, Ducastelle EAM Au whisker aligned with z-axis, periodic boundary in z-direction, $\Delta t=5fs, \tau=1ps, \epsilon=10^8, T=300K$"
    v0, (xs, ys) = list(data.items())[0]
    for i in range(len(xs)):
        dpi = 100
        fig, ax = plt.subplots(figsize=(1920/dpi, 500/dpi))
        fig.text(0.99, 0.01, info, horizontalalignment='right', fontsize="10")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        fig.tight_layout()

        for v0, (xs, ys) in data.items():
            xs = [x/float(v0) for x in xs]
            ax.plot(xs, ys)
            ax.axvline(xs[i], c="red", alpha=0.4, lw=3)
        fig.tight_layout()
        print(i, end="\r")
        fig.savefig("../videos/whisker_51500_graph/" +
                    "{:05d}".format(i)+".png", dpi=dpi)
        plt.close()
    print()
    run("ffmpeg -framerate 60 -pattern_type glob -i '../videos/whisker_51500_graph/*.png' -an -preset veryslow -y ../videos/graph.mp4")
    run("ffmpeg -i ../videos/whisker_51500_3views.mp4 -i ../videos/graph.mp4 -filter_complex vstack=inputs=2 -an -preset veryslow -y ../videos/whisker_51500_graph.mp4")


if __name__ == "__main__":
    plots_optimal_eam_timestep()
    plots_optimal_ljds_timestep()
    plots_runtime_ljds_ljts()
    # plot_gold_temp_over_energy()
    plot_gold_temp_over_energy_fitted()
    plot_parallel_energy_conserved()
    stress_strain()

    # stress_strain_animation()
