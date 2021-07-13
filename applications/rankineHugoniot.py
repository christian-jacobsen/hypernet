import numpy as np
import hypernet.src.general.const as const
from hypernet.src.thermophysicalModels.specie import Specie
from hypernet.src.thermophysicalModels.specie import Thermo

from scipy import optimize


import os
import sys
import tensorflow as tf
import numpy as np
from shutil import copy

def fun():
    return
def jac():
    return

def main(*args):

    if inp not in args:
        if len(sys.argv) != 2:
            raise ValueError
        else:
            utils.print_main("Reading input file")
            inp = sys.argv[1]

    # initilize specie mixture
    if mix not in args:
        utils.print_main("Initilizing thermo and chemistry models")
        mix = MultiComponent(
            inp.freestream['mixture'],
            **inp.specie
        )

    # Conservative quantities
    M, Q, E = jumpRelations(inp.freestream)

    # Rankine-hugoiot
    y0 = jumpRelations(inp.freestream)

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


