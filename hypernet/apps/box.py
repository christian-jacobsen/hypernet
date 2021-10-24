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


NAME = 'box'

import os
import sys
import argparse
import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from hypernet.src.general import utils

from hypernet.src.algorithms.ode import ODE
from hypernet.src.thermophysicalModels import specie as specieMdl
from hypernet.src.thermophysicalModels.reactionThermo import mixture as mixMdl
from hypernet.src.thermophysicalModels import chemistry as chemMdl

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
        description="Solve 0D chemical reactor process.",
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
        help='plot results.'
    )
    parser.add_argument('-v', '--verbose',
        type=int,
        default=1,
        choices=[0,1],
        help='verbose mode'
    )
    return parser.parse_args()

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
    # from inputs import postprocessing as inp_post

    # Initilizing thermophysical models =======================================
    # Species thermo ----------------------------------------------------------
    utils.print_main("Initializing species thermo ...", verbose=opts.verbose)
    thermo = inp_gen.thermo if hasattr(inp_gen, 'thermo') else 'CoupledEnergyModes'
    EOS = inp_gen.EOS if hasattr(inp_gen, 'EOS') else 'PerfectGas'
    spTh = {
        name: specieMdl.specieThermos.SpecieThermos(
            name,
            info,
            thermo=thermo,
            EOS=EOS,
        ) for name, info in inp_gen.specie.items()
    }

    # Mixture -----------------------------------------------------------------
    utils.print_main("Initializing mixture ...", verbose=opts.verbose)
    mix = utils.get_class(mixMdl, inp_gen.mixture['name'])(spTh)
    # Mixture mass/molar fractions
    utils.print_submain("Update initial composition ...", verbose=opts.verbose)
    mix.update(inp_gen.mixture['composition'], var=inp_gen.mixture['var'])
    # Mass
    utils.print_submain("Get mass of the mixture ...", verbose=opts.verbose)
    rho = mix.rho_(**inp_gen.ambient)
    mix.rho_i(rho=rho)

    # Chemistry ---------------------------------------------------------------
    utils.print_main("Initializing chemistry ...", verbose=opts.verbose)
    chem = utils.get_class(chemMdl.chemistrySolver, inp_gen.chemistry['solver'])(
        mixture=mix,
        specieThermos=spTh,
        chemistryModel=inp_gen.chemistry['model'],
        reactionsList=inp_gen.chemistry['reactionsList'],
        processFlags=inp_gen.chemistry['processFlags'],
        constPV=inp_gen.chemistry['constPV'],
    )

    # Initilizing solver ======================================================
    utils.print_main("Initializing ODE solver ...", verbose=opts.verbose)
    solver = ODE(
        inp_gen.ode,
        function=chem.function,
        jacobian=None#chem.jacobian
    )

    # Solving =================================================================
    # Time-marching solution
    utils.print_submain("Solving ...", verbose=opts.verbose)
    T = utils.convert_to_array(inp_gen.ambient['T'])
    y0 = np.concatenate(
        tuple([mix.spTh['O2'].specie.rho, mix.spTh['O'].specie.rho, T])
    )

    # chem.chemModel.update(T)
    # y0 = np.concatenate(
    #     tuple([mix.spTh['O2'].specie.rho, mix.spTh['O'].specie.rho])
    # )
    t, y = solver.solve(y0, args=(rho,))

    # print(Y)
    plt.semilogx(t, y[:,:-1]/rho)
    plt.show()

    plt.semilogx(t, y[:,-1])
    plt.show()


    # # Postprocessing ==========================================================
    # utils.print_main("Postprocessing solution ...")
    # if not os.path.exists(inp_case.postprocess):
    #     os.makedirs(inp_case.postprocess)

    # # Plot `T` and `u`
    # var = ['T','u']
    # plot_var(
    #     inp_case.postprocess+'T_u.pdf',
    #     x_true,
    #     x_pred,
    #     data[var].values,
    #     y_pred[var].values,
    #     var_name=var,
    #     x_label=r'$x\quad[m]$',
    #     y_label=[r'$T\quad[K]$', r'$u\quad[m/s]$'],
    #     scales=['log', 'linear']
    # )

    # # Plot `Y`
    # Y_pred = chem.net.predict(x_pred)
    # Y_pred[0], x_pred[0] = Y[1], x_true[1]
    # O2_names = [ r'$O_2^{({%s})}$' % (i+1) for i in range(3) ]
    # y_scale = 'linear'
    # plot_Y(
    #     inp_case.postprocess+'Y_'+y_scale+'.pdf',
    #     x_true,
    #     x_pred,
    #     Y,
    #     Y_pred,
    #     O2_names+[r'$O$'],
    #     labels=[r'$x\quad[m]$', r'$Y$'],
    #     scales=['log', y_scale]
    # )


