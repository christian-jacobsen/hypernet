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


def fun(y, *args):
    # Conserved quantities
    M, Q, E = cons
    # Variables
    p, T, u = y[0], y[1], y[2]
    psi = 1./(mix.R()*T)
    return np.array([
        p * psi * u - M,
        p * (u**2*psi + 1.) - Q,
        mix.h(T) + 1./2. * u**2 - E
    ])

def jac(y, *args):
    # Variables
    p, T, u = y[0], y[1], y[2]
    psi = 1./(mix.R()*T)
    return np.array([
        [ u*psi, -p*psi*u/T, p*psi ],
        [ u**2*psi+1., -p*psi*u**2/T, 2*p*psi*u ],
        [ 0., mix.cp(T), u ]
    ])

def conserved(mix, p, T, u):
    rho = p/(mix.R()*T)
    M = rho * u
    Q = p + M * u
    H = M * (mix.h(T) + 1./2. * u**2)
    return M, Q, H


@utils.timing
def main(*args):

    utils.print_main("Loading data")
    data = pd.read_csv(inp.data, sep="  ", header=None, \
        skiprows=2, dtype=np.float64, engine="python")
    data.columns = inp.columns
    data['X_O2'] = 1. - data['X_O']

    # for col in inp.columns:
    #     if col.startswith('X'):
    #         name = col.split('_')[-1]
    #         inp.mixture['species'][name] = data[col].to_numpy()[0]

    # Read input file
    # if inp not in args:
    #     if len(sys.argv) != 2:
    #         utils.raise_value_err(
    #             "Define the input file. Usage:\n"
    #             "$ python3 <application> <path_to_input>"
    #         )
    #     else:
    #         utils.print_main("Reading input file")
    #         inp = sys.argv[1]

    # Initialize specie mixture
    utils.print_main("Initializing thermo")
    mix = utils.get_class(mixture_module, inp.mixture['name'])(
        inp.mixture['species'], inp.thermo, **inp.specie
    )
    mix.update(inp.mixture['species'], var=inp.mixture['var'])

    # Conserved quantities
    utils.print_main("Initializing conserved flow quantities")
    # cons = conserved(mix, **inp.freestream)
    cons = conserved(mix, *data[['p', 'T', 'u']].values[0])
    print(cons)
    input('==============================')

    # Set up solver
    utils.print_main("Setting up solver")
    solver = root.Root(inp.algorithm)
    solver.fun = fun
    solver.jac = jac

    # Solve jump relations
    utils.print_main("Solving")
    y0 = np.array(list(inp.freestream.values()))
    y = solver.step(y0, cons, mix)
    y = np.stack((y, mix.rho(p=y[0], T=y[1])))

    return y


if __name__ == "__main__":
    main()




    if inp is None:
        if len(sys.argv) != 2:
            raise ValueError
        else:
            utils.print_main("Reading input file")
            inp = sys.argv[1]

    # initilize specie mixture
    utils.print_main("Initilizing thermo and chemistry models")
    mix = MultiComponent(inp.mixture, **inp.specie)

    # Conservative quantities
    M, Q, E = jumpRelations(inp.freestream)

    # Rankine-hugoiot
    y0 = jumpRelations.main(inp=inp)

    # fit Y spline
    utils.print_main("")

    # Algorithm
    utils.print_main("Setting up solver")
    algo = Root(inp.algorithm)
    algo.fun = fun()
    algo.jac = jac()
    algo.get_args = get_args()

    # solve
    utils.print_main("Solving")
    algo.solve()

if __name__ == "__main__":
    main()
