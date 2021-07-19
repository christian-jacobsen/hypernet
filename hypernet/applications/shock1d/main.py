import os
import sys
import numpy as np
import pandas as pd
import inputs.chemistry as inp_chem
import inputs.case as inp_case
# Remove the following line if you install `hypernet` as a Python package -----
sys.path.append(inp_case.hypernet)
# import hypernet as hy

from matplotlib import pyplot as plt
from hypernet.src.general import const
from hypernet.src.general import utils
from hypernet.src.algorithms import root
from hypernet.src.thermophysicalModels.reactionThermo import mixture as mixture_module
from hypernet.src.thermophysicalModels.chemistryModel import surrogate as surrogate_module


# Solver methods
###############################################################################
def fun(y, *args):
    cons, mix, chem = args
    # Conserved quantities
    M, Q, H = cons
    # Variables
    T, u = y[0], y[1]
    return np.array([
        T/u * mix.R + u - Q/M,
        mix.h(T) + 0.5 * u**2 - H/M
    ])

def jac(y, *args):
    cons, mix, chem = args
    # Variables
    T, u = y[0], y[1]
    return np.array([
        [ mix.R/u, -T/u**2*mix.R+1. ],
        [ mix.cp(T), u ]
    ])

def update_args(x, *args):
    cons, mix, chem = args
    mix.update(chem.update(x))
    return cons, mix, chem

# Thermo methods
###############################################################################
def X_to_Y(mix, X):
    n = X.shape[0]    
    m_sp = []
    for specie, thermo in mix.mixture.items():
        if hasattr(thermo.specie, 'n_bins') and thermo.specie.n_bins > 1:
            m_sp.extend( [ thermo.specie.m for _ in range(thermo.specie.n_bins) ]
            )
        else:
            m_sp.append(thermo.specie.m)
    m_sp = np.tile(np.array(m_sp), (n, 1))
    m_mix = np.sum(X * m_sp, axis=1, keepdims=True)
    return X*m_sp/m_mix

def conserved(mix, p, T, u):
    rho = p/(mix.R*T)
    M = rho * u
    Q = p + M * u
    H = M * (mix.h(T) + 0.5 * u**2)
    return M, Q, H

# Plotting methods
###############################################################################
def plot_test(
    path,
    x_true,
    x_pred,
    y_true,
    y_pred,
    var_name,
    title=None,
    labels=[None, None],
    scales=['log', 'linear'],
    data_val=None
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
    if x_scale != 'log':
        plt.xlim([min(x_true), max(x_true)])

    # Y axis
    if y_label is not None:
        plt.ylabel(y_label)
    if y_scale == 'log':
        plt.yscale(y_scale)
        plt.ylim([np.amin(y_true)*5.e-1, np.amax(y_true)*5.e+0])
    else:
        delta = np.amax(y_true)*0.1
        plt.ylim([np.amin(y_true)-delta, np.amax(y_true)+delta])

    # Plotting
    lin = ['-', '--']
    col = ['k', 'r', 'g', 'b', 'tab:purple', 'tab:brown', 'tab:pink', \
           'tab:blue', 'tab:red', 'tab:green', 'tab:orange']
    marker_style = dict(marker='^', fillstyle='none', markersize=5)#, markevery=y_true.shape[0]//20)

    # Solution
    for d in range(y_true.shape[1]):

        name = 'True ' + var_name[d]
        plt.plot(x_true, y_true[:, d], c=col[d], ls=lin[0], lw=1, label=name)

        name = 'Pred ' + var_name[d]
        plt.plot(x_pred, y_pred[:, d], c=col[d+1], ls=lin[1], lw=1, label=name, **marker_style)

    plt.legend(fontsize='x-small')
    fig.savefig(path)
    plt.close()


# Main function
###############################################################################
@utils.main_decorator
def main(*args):

    # Loadind data ============================================================
    utils.print_main("Loading data")
    data = pd.read_csv(inp_chem.path_to_data, sep="  ", header=None, \
        skiprows=2, dtype=np.float64, engine="python")
    data.columns = inp_chem.columns

    # Initilizing thermophysical models =======================================
    # Thermodynamic quantities
    utils.print_main("Initializing thermo")
    mix = utils.get_class(mixture_module, inp_case.mixture['name'])(
        inp_case.mixture['species'], inp_case.thermo, **inp_case.specie
    )
    # Mixture mass/molar fractions
    mix.update(inp_case.mixture['species'], var=inp_case.mixture['var'])
    # >> Get data
    X = np.hstack(
        [ np.expand_dims(data[col], axis=-1) \
            for col in inp_chem.columns if col.startswith('X') ]
    )
    Y = X_to_Y(mix, X)

    # Computing conserved flow quantities
    utils.print_main("Computing conserved flow quantities")
    cons = conserved(mix, **inp_case.freestream)

    # Chemistry model
    utils.print_main("Initializing chemistry")
    # >> Define IC: [M, rho, T, Y]
    x = np.expand_dims(data['x'], axis=-1)
    # ic = np.tile(np.array(
    #     [data['rho'][0]*data['u'][0], data['rho'][0], data['T'][0], *Y[0]]), \
    #         (x.shape[0], 1)
    # )
    # train = [np.hstack((ic,x)), Y]
    train = [x, Y]
    chem = surrogate_module.DeepNet(
        mix.mixture,
        inp_chem,
        train=train,
        test=[train]
    )
    if inp_chem.train_flg:
        chem.fit()

    # Solving =================================================================
    # Set up solver
    utils.print_main("Setting up solver")
    solver = root.Root(inp_case.algorithm)
    solver.fun = fun
    solver.jac = jac
    solver.update_args = update_args

    # Space-marching solution
    utils.print_main("Solving")
    x_true = np.expand_dims(data['x'].values, axis=-1)
    y_true = data[['T','u']]
    y0 = y_true.values[0]
    x_pred, y_pred = solver.solve(y0, cons, mix, chem)
    y_pred = pd.DataFrame(data=y_pred, index=None, columns=['T','u'])

    # Postprocessing ==========================================================
    utils.print_main("Postprocessing solution")
    if not os.path.exists(inp_case.postprocess):
        os.makedirs(inp_case.postprocess)
    for var in ['T','u']:
        unit = 'K' if var == 'T' else 'm/s'
        plot_test(
            inp_case.postprocess+var+'.pdf',
            x_true,
            x_pred,
            np.expand_dims(y_true[var].values, axis=-1),
            np.expand_dims(y_pred[var].values, axis=-1),
            var_name=[var],
            labels=[r'$x\quad[m]$', r'$%s\quad[%s]$' % (var, unit)],
            scales=['log', 'linear']
        )


if __name__ == "__main__":
    main()
