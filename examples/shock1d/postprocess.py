import numpy as np

from matplotlib.lines import Line2D
from matplotlib import pyplot as plt


def plot_var(
    path,
    x_true,
    x_pred,
    y_true,
    y_pred,
    var_name,
    title=None,
    x_label=None,
    y_label=None,
    scales=['log', 'linear'],
    x_lim=[1.e-9, 1.e-2]
):
    """Variable plotting."""

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.grid(True)

    if title:
        ax1.set_title(title)

    x_scale, y_scale = scales

    # X axis
    if x_label is not None:
        ax1.set_xlabel(x_label)
    if x_scale is not None:
        ax1.set_xscale(x_scale)
    ax1.set_xlim(x_lim)

    # Plotting
    lin = ['-', '--']
    col = ['r', 'b', 'g', 'tab:purple', 'tab:brown', 'tab:pink', \
           'tab:blue', 'tab:red', 'tab:green', 'tab:orange']
    marker_style = dict(marker='^', fillstyle='none', markersize=5)

    # Solution
    ax = [ax1, ax2]
    curves = []
    for d in range(y_true.shape[1]):

        if y_label[d] is not None:
            ax[d].set_ylabel(y_label[d])
            # ax[d].set_ylabel(y_label[d], color=col[d])
        # ax[d].tick_params(axis='y', labelcolor=col[d])

        delta = np.amin(y_true[:, d])*0.1*(1+4*d)
        ax[d].set_ylim([np.amin(y_true[:, d])-delta, np.amax(y_true[:, d])+delta])

        name = 'True ' + var_name[d]
        curves += ax[d].plot(x_true, y_true[:, d], c=col[d], ls=lin[0], lw=1, label=name)

        name = 'Pred ' + var_name[d]
        curves += ax[d].plot(x_pred, y_pred[:, d], c=col[d], ls=lin[1], lw=1, label=name, **marker_style)

    # legend_elements = [
    #     Line2D([0], [0], color='k', ls=lin[0], lw=1, label='True'),
    #     Line2D([0], [0], color='k', ls=lin[1], lw=1, label='Pred', marker='^', fillstyle='none')
    # ]
    # ax1.legend(handles=legend_elements)

    labels = [c.get_label() for c in curves]
    ax1.legend(curves, labels)

    fig.savefig(path)
    plt.close()


def plot_Y(
    path,
    x_true,
    x_pred,
    y_true,
    y_pred,
    var_name,
    title=None,
    labels=[None, None],
    scales=['log', 'linear'],
    x_lim=[1.e-9, 1.e-2],
    y_lim=[1.e-8, 5.]
):
    """Variable plotting."""

    fig = plt.figure()
    plt.grid()

    if title:
        plt.title(title)

    x_label, y_label = labels
    x_scale, y_scale = scales

    # X axis
    if x_label is not None:
        plt.xlabel(x_label)
    if x_scale is not None:
        plt.xscale(x_scale)
    plt.xlim(x_lim)

    # Y axis
    if y_label is not None:
        plt.ylabel(y_label)
    if y_scale == 'log':
        plt.yscale(y_scale)
        plt.ylim([np.amin(y_true)*5.e-1, np.amax(y_true)*5.e+0])
    else:
        delta = np.amax(y_true)*0.1
        plt.ylim([np.amin(y_true)-delta, np.amax(y_true)+delta])
    # plt.ylim(y_lim)

    # Plotting
    lin = ['-', '--']
    col = ['k', 'r', 'g', 'b', 'tab:purple', 'tab:brown', 'tab:pink', \
           'tab:blue', 'tab:red', 'tab:green', 'tab:orange']
    marker_style = dict(marker='^', fillstyle='none', markersize=5)

    # Solution
    if len(y_true.shape) < 2:
        y_true = np.expand_dims(y_true, 1)

    for d in range(y_true.shape[1]):

        name = 'True ' + var_name[d]
        plt.plot(x_true, y_true[:, d], c=col[d], ls=lin[0], lw=1, label=name)

        name = 'Pred ' + var_name[d]
        plt.plot(x_pred, y_pred[:, d], c=col[d], ls=lin[1], lw=1, label=name, **marker_style)

    plt.legend(fontsize='x-small')
    fig.savefig(path)
    plt.close()


def plot_levels(
    path,
    n_lev,
    E_lev,
    g_lev,
    var_name,
    title=None,
    labels=[r'$\epsilon_i\quad[eV]$', r'$n_i\//\/g_i\quad[m^{-3}]$'],
    scales=['linear', 'log']
):
    """Variable plotting."""

    fig = plt.figure()
    plt.grid()

    if title:
        plt.title(title)

    x_label, y_label = labels
    x_scale, y_scale = scales

    # X axis
    if x_label is not None:
        plt.xlabel(x_label)
    if x_scale is not None:
        plt.xscale(x_scale)

    # Y axis
    if y_label is not None:
        plt.ylabel(y_label)
    if x_scale is not None:
        plt.yscale(y_scale)

    # Plotting
    col = ['k', 'r', 'g', 'b', 'tab:purple', 'tab:brown', 'tab:pink', \
           'tab:blue', 'tab:red', 'tab:green', 'tab:orange']
    marker_style = dict(marker='o', fillstyle='none', markersize=1, linestyle='none')

    for g in range(len(n_lev)):
        plt.plot(E_lev[g], n_lev[g]/g_lev[g], c=col[g], lw=1, label=var_name[g], **marker_style)

    plt.legend()
    fig.savefig(path)
    plt.close()
