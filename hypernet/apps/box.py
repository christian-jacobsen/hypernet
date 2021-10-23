# *************************************************************************** #
#                                    PrODE                                    #
# --------------------------------------------------------------------------- #
#                     Machine-Learning-Based Approximators                    #
#                  for Ordinary Differential Equations (ODEs)                 #
#                                                                             #
# *************************************************************************** #

# -------------------------------- EXEC FILE -------------------------------- #
# Description:
# >> Build a ML-based approximator for ODEs
# --------------------------------------------------------------------------- #


NAME = 'fitRates'

import os
import sys
import argparse
import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from dataGenerator import DataGenerator
from hypernet.src.general import utils
from hypernet.src.thermophysicalModels.specie.specie import Specie

import hypernet.database as db
kinetic_db = os.path.dirname(db.__file__) + '/air/kinetics/'


# Parse arguments
###############################################################################
def dir_path(path):
    if os.path.exists(path+'/inputs'):
        return path
    else:
        raise argparse.ArgumentTypeError(
            "'inputs' folder not found in '{}'".format(
                os.path.normpath(path)+'/'
            )
        )

def get_opts():
    parser = argparse.ArgumentParser(
        prog=NAME,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Fit reactions rates following the Arrhenius law.",
        epilog=utils.app_epilog(name=NAME)
    )
    parser.add_argument('-d', '--dir',
        type=dir_path,
        default="./",
        help='path to the working directory (with `inputs` folder)'
    )
    parser.add_argument('-p', '--plot',
        type=int,
        default=1,
        choices=[0,1],
        help='plot fitted rates.'
    )
    parser.add_argument('-v', '--verbose',
        type=int,
        default=1,
        choices=[0,1],
        help='verbose mode'
    )
    return parser.parse_args()

# Arrhenius Law
###############################################################################
def log_arrhenius_law(T, a_log, b, c):
    # a_log = np.log(a)
    return a_log + b*np.log(T) - c/T

def modified_arrhenius_inv(T_inv, a, b, c):
    # T_inv = 1/T
    return a * np.power(1/T_inv,b) * np.exp(-c*T_inv)

# Plot rates
###############################################################################
def plot_rates(
    fig_name,
    x_true,
    x_pred,
    y_true,
    y_pred,
    var_name,
    title,
    labels=[None, None],
    scales=['log', 'linear']):
    """Variable plotting."""

    fig = plt.figure()
    plt.grid()
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
    if y_scale == 'log':
        plt.yscale(y_scale)
        plt.ylim([np.amin(y_true[y_true!=0.])*1.e-1, np.amax(y_true)*1.e+1])
    else:
        delta = np.amax(y_true)*0.1
        plt.ylim([np.amin(y_true)-delta, np.amax(y_true)+delta])

    # Solution
    for d in range(y_true.shape[1]):

        name = 'True ' + var_name[d]
        marker_style = dict(marker='x', lw=0.5, \
                            linestyle='none', fillstyle='none', markersize=5)
        plt.plot(
            x_true,
            y_true[:,d],
            **marker_style,
            label=name
        )

        name = 'Pred ' + var_name[d]
        plt.plot(
            x_pred,
            y_pred[:,d],
            c=plt.gca().lines[-1].get_color(),
            ls='--',
            lw=1,
            label=name
        )

    plt.legend(fontsize='x-small')
    fig.savefig(fig_name)
    plt.close()

# Main
###############################################################################
@utils.app_decorator(name=NAME)
def main():

    # Inizialization ==========================================================
    # Parse arguments ---------------------------------------------------------
    opts = get_opts()

    # Input module ------------------------------------------------------------
    sys.path.append(opts.dir)
    from inputs import general as inp_gen
    from inputs import postprocessing as inp_post

    # Initialize species ------------------------------------------------------
    utils.print_main("Initializing species ...", verbose=opts.verbose)
    species = {
        sp: Specie(sp, **info) if info != None else Specie(sp) \
            for sp, info in inp_gen.species.items()
    }

    # Generate Data ===========================================================
    utils.print_main('Generating data ...', verbose=opts.verbose)
    dataGen = DataGenerator(
        species,
        inp_gen.T,
        inp_gen.reacReader,
        inp_gen.reacWriter,
        verbose=opts.verbose
    )
    T, K = dataGen.training_data()

    # Fit Data ================================================================
    utils.print_main('Fitting rates ...', verbose=opts.verbose)
    for j in range(K.shape[1]):
        K_log_j = np.log( K[:,j][ K[:,j] != 0. ] )
        T_j     = T[:,0][ K[:,j] != 0. ]

        param_j, _ = curve_fit(log_arrhenius_law, T_j, K_log_j, \
            p0=[1,1,1.e4], method='trf')
        param_j[0] = np.exp(param_j[0])
        if j == 0:
            param = np.array(param_j)
        else:
            param = np.vstack((param, np.array(param_j)))

    # Write coefficient =======================================================
    paramDB = pd.DataFrame(param, columns=['A', 'n', 'Ta'])
    reactions = pd.concat([dataGen.reacDB, paramDB], axis=1)
    path = kinetic_db + inp_gen.reacReader['path']+'/reactions.csv'
    reactions.to_csv(path, float_format='{:e}'.format)

    print(reactions)
    input('=============================')

    # Plot fitted rates
    ###########################################################################
    if opts.plot:
        utils.print_main('Plotting rates ...', verbose=opts.verbose)
        x_true = np.reciprocal(np.array(inp_gen.T))
        y_true = K

        n = 1000
        x_pred = np.reciprocal(
            np.linspace(1000.0e0, max(inp_gen.T), n, dtype=np.float64)
        )
        y_pred = np.zeros((n, dataGen.n_reac), dtype=np.float64)
        for i, param_i in enumerate(param):
            y_pred[:,i] = modified_arrhenius_inv(x_pred, *param_i)

        path = opts.dir + '/plots/'
        if not os.path.exists(path):
            os.makedirs(path)

        # Dissociation
        start = 0
        end = dataGen.n_diss
        if y_pred.shape[1] >= end:
            y_true_diss = y_true[:,start:end]
            y_pred_diss = y_pred[:,start:end]
            fig_name = path + 'dissociation.pdf'
            plot_rates(
                fig_name,
                x_true,
                x_pred,
                y_true_diss,
                y_pred_diss,
                inp_post.var_names['diss'],
                'Dissociation Rates',
                labels=inp_post.labels,
                scales=inp_post.scales
            )

        # Exchange
        start = dataGen.n_diss
        end = dataGen.n_diss + dataGen.n_excit
        if y_pred.shape[1] >= end:
            y_true_exch = y_true[:,start:end]
            y_pred_exch = y_pred[:,start:end]
            fig_name = path + 'exchange.pdf'
            plot_rates(
                fig_name,
                x_true,
                x_pred,
                y_true_exch,
                y_pred_exch,
                inp_post.var_names['excit'],
                'Exchange Rates',
                labels=inp_post.labels,
                scales=inp_post.scales
            )

        # Inelastic
        start = dataGen.n_diss + dataGen.n_excit
        end = dataGen.n_diss + dataGen.n_excit * 2
        if y_pred.shape[1] >= end:
            y_true_inel = y_true[:,start:end]
            y_pred_inel = y_pred[:,start:end]
            fig_name = path + 'inelastic.pdf'
            plot_rates(
                fig_name,
                x_true,
                x_pred,
                y_true_inel,
                y_pred_inel,
                inp_post.var_names['excit'],
                'Inelastic Rates',
                labels=inp_post.labels,
                scales=inp_post.scales
            )



if __name__ == "__main__":
    main()
