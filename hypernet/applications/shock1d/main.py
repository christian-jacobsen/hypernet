import sys
import numpy as np
import pandas as pd
import inputs.chemistryModel as inp_chem
import inputs.case as inp_case
# Remove the following line if you install `hypernet` as a Python package -----
sys.path.append(inp_case.hypernet)
# import hypernet as hy

from hypernet.src.general import const
from hypernet.src.general import utils
from hypernet.src.algorithms import root
from hypernet.src.thermophysicalModels.reactionThermo import mixture as mixture_module
from hypernet.src.thermophysicalModels.chemistryModel import surrogate as surrogate_module


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
    chem.update(x)
    return cons, mix, chem

# Thermo methods ##############################################################
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

# Main function ###############################################################
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
    ic = np.tile(np.array(
        [data['rho'][0]*data['u'][0], data['rho'][0], data['T'][0], *Y[0]]), \
            (x.shape[0], 1)
    )
    train = [np.hstack((ic,x)), Y]
    chem = surrogate_module.DeepNet(
        mix,
        inp_chem,
        train=train,
        test=train
    )
    
    chem.fit()

    input('===============================')

    # Solving =================================================================
    # Set up solver
    utils.print_main("Setting up solver")
    solver = root.Root(inp_case.algorithm)
    solver.fun = fun
    solver.jac = jac
    solver.update_args = update_args

    # Space-marching solution
    utils.print_main("Solving")
    y0 = np.array(list(inp_case.freestream.values()))
    x, y = solver.solve(y0, cons, mix, chem)


if __name__ == "__main__":
    main()
