import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from . import pop # for pop.compute_population

def plot_scalar(
    filename,
    traj,
    key,
    ylabel=None,
    time_units='au',
    legend_loc=1,
    state_colors=None,
    plot_average=True,
    clf=True,
    ):

    """ Plot an AIMS scalar property (e.g., a "spaghetti plot" for a bond distance).

    Params:
        filename (str) - the output PDF file
        traj (Trajectory) - the Trajectory to plot properties for
        key (str) - the key of the property
        ylabel (str) - the label of the y-axis (defaults to key)
        time_units (str) - "au" or "fs"
        legend_loc (legend location) - location indicator for legend
        state_colors (list of colors) - list of colors to use for state
            plotting. None defaults to interpolation on jet colormap. For two
            or three states, many people prefer ['r', 'b', 'g'] or similar.
        plot_average (bool) - Plot averages for each state and total?
        clf (bool) - clear plot or not?
    Result/Returns:
        returns plt handle for further modification
        saves figure to filename
    """

    if time_units == 'au':
        time_scale = 1.0
    elif time_units == 'fs':
        time_scale = 1.0 / 41.3413745758  # TODO: standardize this
    else:
        raise ValueError('Unknown time units: %s' % time_units)

    if state_colors:
        colors=state_colors
    else:
        cmap = matplotlib.cm.get_cmap('jet')
        colors = [cmap(float(x) / (len(traj.Is) - 1)) for x in reversed(list(range(len(traj.Is))))]

    if clf: plt.clf()
    if plot_average:
        # Plot average
        plt.plot(time_scale * np.array(traj.ts), traj.extract_property(key), '-k', linewidth=3.0)
        # Plot state averages
        for Iind, I in enumerate(traj.Is):
            traj2 = traj.subset_by_I(I)
            color = colors[Iind]
            plt.plot(time_scale * np.array(traj2.ts), traj2.extract_property(key), '-', color=color, linewidth=2.0)
    # Plot individual frames
    for Iind, I in enumerate(traj.Is):
        traj2 = traj.subset_by_I(I)
        color = colors[Iind]
        for lind, label in enumerate(traj2.labels):
            traj3 = traj2.subset_by_label(label)
            plt.plot(time_scale * np.array(traj3.ts), traj3.extract_property(key), '-', color=color, linewidth=1.0, label=('State=%d' % I if lind == 0 else None))

    plt.xlabel('t [%s]' % time_units)
    plt.ylabel(ylabel if ylabel else key)
    plt.axis('tight') # TODO: Does not seem to respect this
    plt.legend(loc=legend_loc)
    plt.savefig(filename)
    return plt

def plot_vector(
    filename,
    traj,
    key,
    y,
    ylabel=None,
    time_units='au',
    diff=False,
    cmap=plt.cm.bwr, 
    levels=None, 
    nlevel=65, 
    twosided=True, 
    clf=True,
    ):

    """ Plot an AIMS vector property (e.g., a heatmap of a UED or PES signal).

    Params:
        filename (str) - the output PDF file
        traj (Trajectory) - the Trajectory to plot properties for
        key (str) - the key of the property
        y (np.ndarray) - the indices of the y axis (e.g., R or Q or something
            like it)
        ylabel (str) - the label of the y-axis (defaults to key)
        time_units (str) - "au" or "fs"
        diff (bool) - is this a difference property from t=0?
        cmap (colormap) - the desired colormap
        levels (np.ndarray) - the explicitly desired contour levels (1st
            priority).
        nlevel (np.ndarray) - number of evenly spaced contour levels to
            saturate data (2nd priority).
        twosided (bool) - is the colormap two-sided (used only with nlevel)
        clf (bool) - clear plot or not?
    Result/Returns:
        returns plt handle for further modification
        saves figure to filename
    """
    
    if time_units == 'au':
        time_scale = 1.0
    elif time_units == 'fs':
        time_scale = 1.0 / 41.3413745758  # TODO: standardize this
    else:
        raise ValueError('Unknown time units: %s' % time_units)

    # Extract data to plot
    ts = np.array(traj.ts)
    ys = y
    Vs = traj.extract_property(key)

    # Difference properties
    if diff:
        Vs -= np.outer(np.ones_like(ts), Vs[0,:])

    # Levels and color ticks
    if levels is None:
        vmax = np.max(np.abs(Vs))
        if twosided:
            levels = np.linspace(-vmax,+vmax,nlevel)
            cticks = [-int(vmax),0,+int(vmax)]
        else:
            levels = np.linspace(0,+vmax,nlevel)
            cticks = [0,+int(vmax)]
    else:
        vmax = np.max(levels)
        if twosided:
            cticks = [-int(vmax),0,+int(vmax)]
        else:
            cticks = [0,+int(vmax)]

    if clf: plt.clf()
    hs = plt.contourf(time_scale * ts,ys,Vs.T,levels=levels,cmap=cmap,extend='both')
    for h in hs.collections:
        h.set_edgecolor('face')
    plt.colorbar(ticks=cticks)
    plt.xlabel('t [%s]' % time_units)
    plt.ylabel(ylabel if ylabel else key)
    plt.axis('tight') # TODO: Does not seem to respect this
    plt.savefig(filename)
    return plt

def plot_population(
    filename,
    traj,
    trajs,
    time_units='au',
    legend_loc='right',
    state_colors=None,
    plot_total=True,
    clf=True,
    tmax=None,
    ):

    """ Plot the AIMS state populations.

    Params:
        filename (str) - the output PDF file
        traj (Trajectory) - the Trajectory to plot populations for
        trajs (list of Trajectory) - a list of Trajectory objects, one for each
            IC. This is used for the "spaghetti" plots.
        time_units (str) - "au" or "fs"
        legend_loc (legend location) - location indicator for legend
        state_colors (list of colors) - list of colors to use for state
            plotting. None defaults to interpolation on jet colormap. For two
            or three states, many people prefer ['r', 'b', 'g'] or similar.
        plot_total (bool) - Plot total population?
        clf (bool) - clear plot or not?
    Result/Returns:
        returns plt handle for further modification
        saves figure to filename
    """

    if time_units == 'au':
        time_scale = 1.0
    elif time_units == 'fs':
        time_scale = 1.0 / 41.3413745758  # TODO: standardize this
    else:
        raise ValueError('Unknown time units: %s' % time_units)

    if state_colors:
        colors=state_colors
    else:
        cmap = matplotlib.cm.get_cmap('jet')
        colors = [cmap(float(x) / (len(traj.Is) - 1)) for x in reversed(list(range(len(traj.Is))))]

    if clf: plt.clf()

    # Spaghetti plots from trajs
    Imap = { I : Iind for Iind, I in enumerate(traj.Is) }
    for traj2 in trajs:
        ts = np.array(traj2.ts)
        for I in traj2.Is:
            traj3 = traj2.subset_by_I(I)
            t3s = traj3.ts
            w3s = [sum(frame.w for frame in traj3.subset_by_t(t).frames) for t in t3s]
            plt.plot(time_scale * t3s, w3s, '-', color=colors[Imap[I]], linewidth=1.0)
    
    # State populations
    ts = np.array(traj.ts)
    populations = pop.compute_population(traj)
    for Iind, I in enumerate(sorted(populations.keys())):
        plt.plot(time_scale * ts, populations[I], '-', color=colors[Iind], linewidth=2.0, label='State=%d' % I)

    # Total population (to check norm, state cutoffs, etc)
    if plot_total:
        total = np.zeros_like(ts)
        for population in list(populations.values()): total += population
        plt.plot(time_scale * ts, total, '-', color='k', linewidth=2.0, label='Total')
   
    if tmax is None: tmax = max(ts) * time_scale

    plt.xlabel('t [%s]' % time_units)
    plt.ylabel('Population [-]')
    plt.axis([time_scale * min(ts), tmax, -0.1, 1.1])
    plt.legend(loc=legend_loc)
    plt.tight_layout()
    plt.savefig(filename)

    return plt
