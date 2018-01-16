import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot_scalar(
    filename,
    traj,
    key,
    ylabel=None,
    time_units='au',
    legend_loc=1,
    state_colors=None,
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
    Result/Returns:
        returns plt handle for further modification
        saves figure to filename
    """

    plt.clf()

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
        colors = [cmap(float(x) / (len(traj.Is) - 1)) for x in reversed(range(len(traj.Is)))]

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
        for label in traj2.labels:
            traj3 = traj2.subset_by_label(label)
            plt.plot(time_scale * np.array(traj3.ts), traj3.extract_property(key), '-', color=color, linewidth=1.0)

    plt.xlabel('t [%s]' % time_units)
    plt.ylabel(ylabel if ylabel else key)
    plt.axis('tight')
    plt.legend(['Average'] + ['State=%d' % (I) for I in traj.Is], loc=legend_loc)
    plt.savefig(filename)
    return plt

    


