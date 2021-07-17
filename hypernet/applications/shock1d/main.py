import sys
import numpy as np
import pandas as pd
import inputs as inp
# Remove the following line if you install `hypernet` as a Python package -----
sys.path.append(inp.hypernet)
# import hypernet as hy

from hypernet.src.general import const
from hypernet.src.general import utils
from hypernet.src.algorithms import root
from hypernet.src.thermophysicalModels.reactionThermo import mixture as mixture_module


# Solver methods ##############################################################
def fun(y, *args):
    # Conserved quantities
    M, Q, H = cons
    # Variables
    T, u = y[0], y[1]
    return np.array([
        T/u * mix.R + u - Q/M,
        mix.h(T) + 0.5 * u**2 - H/M
    ])

def jac(y, *args):
    # Variables
    T, u = y[0], y[1]
    return np.array([
        [ mix.R/u, -T/u**2*mix.R+1. ],
        [ mix.cp(T), u ]
    ])

def update_args(self, x, *args):
    mix.update(chem.update(x))
    return cons, mix, chem

# Thermo methods ##############################################################
def conserved(mix, p, T, u):
    rho = p/(mix.R*T)
    M = rho * u
    Q = p + M * u
    H = M * (mix.h(T) + 0.5 * u**2)
    return M, Q, H/M

# Main function ###############################################################
@utils.main_decorator
def main(*args):

    # Loadind data ============================================================
    utils.print_main("Loading data")
    data = pd.read_csv(inp.data, sep="  ", header=None, \
        skiprows=2, dtype=np.float64, engine="python")
    data.columns = inp.columns

    # Initilizing thermophysical models =======================================
    # Thermodynamic quantities
    utils.print_main("Initializing thermo")
    mix = utils.get_class(mixture_module, inp.mixture['name'])(
        inp.mixture['species'], inp.thermo, **inp.specie
    )
    # Mixture mass/molar fractions
    mix.update(inp.mixture['species'], var=inp.mixture['var'])

    # Chemistry model
    utils.print_main("Initializing chemistry")
    chem = None # chemistryModel

    # Solving =================================================================
    # Set up solver
    utils.print_main("Setting up solver")
    solver = root.Root(inp.algorithm)
    solver.fun = fun
    solver.jac = jac
    solver.update_args = update_args

    # Computing conserved flow quantities
    utils.print_main("Computing conserved flow quantities")
    cons = conserved(mix, **inp.freestream)

    # Space-marching solution
    utils.print_main("Solving")
    y0 = np.array(list(inp.freestream.values()))
    x, y = solver.solve(y0, cons, mix, chem)


if __name__ == "__main__":
    main()
