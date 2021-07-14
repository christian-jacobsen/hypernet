import numpy as np
import inputs as inp

from hypernet.src.thermophysicalModels.reactionThermo import mixture as mixture_module
from hypernet.src.algorithms import root


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

def conserved(p, T, u, mix):
    psi = 1./(mix.R()*T)
    M = p * psi * u
    Q = p * (u**2*psi + 1.)
    E = mix.h(T) + 1./2. * u**2
    return M, Q, E


@utils.timing()
def main(*args):

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
    if mix not in args:
        utils.print_main("Initializing Thermo model")
        mix = utils.get_class(mixture_module, inp.mixture['name'])(
            inp.mixture['species'], inp.thermo, **inp.specie
        )

    # Conserved quantities
    utils.print_main("Initializing conserved flow quantities")
    cons = conserved(**inp.freestream, mix)

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


