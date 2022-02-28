# *************************************************************************** #
#                                   HyperNet                                  #
# --------------------------------------------------------------------------- #
#                 Machine Learning-Based library for modeling                 #
#           multi-component non-equilibrium thermochemical processes          #
#                                                                             #
# *************************************************************************** #

# -------------------------------- EXEC FILE -------------------------------- #
# Description:
# >> Solve 0D chemical reactor 
# --------------------------------------------------------------------------- #


NAME = 'box'

import os
import sys
import argparse
import numpy as np
import pandas as pd

from hypernet.src.general import utils

from hypernet.src.solvers.ode import ODE
from hypernet.src.thermophysicalModels import specie as specieMdl
from hypernet.src.thermophysicalModels.reactionThermo import mixture as mixMdl
from hypernet.src.thermophysicalModels import chemistry as chemMdl


LABELS = {
    't': 't [s]', 'T': 'T [K]', 'p': 'p [Pa]', 'n': 'n [m^-3]', 'E': 'E [J/kg]'
}


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
        description="Solve 0D chemical reactor",
        epilog=utils.app_epilog(name=NAME)
    )
    parser.add_argument('-d', '--dir',
        type=dir_path,
        default="./",
        help='path to the working directory (with `inputs` folder)'
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
    import inputs

    # Initilizing thermophysical models =======================================
    # Species thermo ----------------------------------------------------------
    utils.print_main(
        "Initializing species thermodynamics ...",
        start='',
        verbose=opts.verbose
    )
    thermo = getattr(inputs.general, 'thermo', 'CoupledEnergyModes')
    EOS = getattr(inputs.general, 'EOS', 'PerfectGas')
    constVP = getattr(inputs.general, 'constVP', 'V')
    spTh = {
        name: specieMdl.specieThermos.SpecieThermos(
            name,
            info,
            thermo=thermo,
            constVP=constVP,
            EOS=EOS
        ) for name, info in inputs.general.specie.items()
    }

    # Mixture -----------------------------------------------------------------
    utils.print_main("Initializing mixture ...", verbose=opts.verbose)
    mix = utils.get_class(mixMdl, inputs.general.mixture['name'])(spTh)
    # Mixture mass/molar fractions
    utils.print_submain("Update initial composition", verbose=opts.verbose)
    mix.update(
        XY=inputs.general.mixture['composition'],
        var=inputs.general.mixture['var']
    )
    # Mass/
    utils.print_submain("Get mixture mass", verbose=opts.verbose)
    rho = mix.rho_(**inputs.general.ambient)
    mix.rhoi_(rho=rho)

    # Chemistry ---------------------------------------------------------------
    utils.print_main("Initializing chemistry ...", verbose=opts.verbose)
    chem = utils.get_class(
        chemMdl.chemistrySolver, inputs.general.chemistry['solver']
    )(
        mixture=mix,
        specieThermos=spTh,
        chemistryModel=inputs.general.chemistry['model'],
        reactionsList=inputs.general.chemistry['reactionsList'],
        processFlags=inputs.general.chemistry['processFlags'],
        heatBath=inputs.general.chemistry['heatBath'],
    )

    # Initilizing solver ======================================================
    utils.print_main("Initializing ODE solver ...", verbose=opts.verbose)
    solver = ODE(
        setup=inputs.general.setup,
        function=chem.function,
        jacobian=chem.jacobian,
        eval_extras=chem.eval_extra_vars
    )

    # Solving =================================================================
    # Time-marching solution
    utils.print_submain("Solving", verbose=opts.verbose)
    y0 = np.concatenate(tuple([
        mix.spTh[chem.chemModel.molecule].specie.rho,
        mix.spTh[chem.chemModel.atom].specie.rho,
        utils.convert_to_array(inputs.general.ambient['T'])
    ]))
    t, y = solver.solve(y0, args=(rho,))
    y[:,:-1] = y[:,:-1] / rho

    # Printing ================================================================
    utils.print_submain("Printing solution", verbose=opts.verbose)
    # Evaluate extra variables
    for yi in y:
        chem.mix_update(yi[:-1])
        chem.eval_extra_vars(yi[-1], rho)
    # Create solution array
    sol = np.hstack([t, y])
    for v in chem.extraVars.values():
        sol = np.hstack([sol, np.array(v).reshape(-1,1)])
    # Create dataframe
    columns = ['t'] + list(chem.varNames) + list(chem.extraVars.keys())
    for i, c in enumerate(columns):
        columns[i] = LABELS[c] if c in LABELS.keys() else c
    columns = tuple(columns)
    df = pd.DataFrame(sol, columns=columns)
    # Write object to a .csv file
    out = opts.dir + '/outputs/'
    if not os.path.exists(out):
        os.makedirs(out)
    df.to_csv(opts.dir + '/outputs/out.csv', index=False)
