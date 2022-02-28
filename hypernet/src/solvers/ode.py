import numpy as np

from scipy.integrate import solve_ivp
from hypernet.src.general import utils
from hypernet.src.solvers import Basic


class ODE(Basic):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        setup,
        function=None,
        jacobian=None,
        *args,
        **kwargs
    ):
        super(ODE, self).__init__(
            setup,
            function=function,
            jacobian=jacobian,
            *args,
            **kwargs
        )
        # Integration method
        self.method = setup['method'].upper() if 'method' in setup else 'BDF'

    # Methods
    ###########################################################################
    @utils.timing
    def solve(self, y0, args=()):
        '''Solving the ODE.'''

        sol = solve_ivp(
            fun=self.fun,
            t_span=(self.t_st,self.t_end),
            y0=y0,
            method=self.method,
            args=args,
            first_step=self.dt_min,
            max_step=self.dt_max,
            atol=1.e-20,
            jac=self.jac
        )
        assert (sol.y > 0.).all()

        t = sol.t.reshape(-1,1)
        y = sol.y.T

        return [t, y]