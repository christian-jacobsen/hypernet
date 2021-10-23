import numpy as np

from scipy import optimize
# from profilehooks import profile
from hypernet.src.general import utils
from hypernet.src.algorithms import Basic


class Root(Basic):

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
        super(Root, self).__init__(
            setup,
            function=function,
            jacobian=jacobian,
            *args,
            **kwargs
        )
        # Root finding method
        self.method = setup['method'] if 'method' in setup else 'hybr'
        
    # Methods
    ###########################################################################
    # @profile
    @utils.timing
    def solve(self, y0, *args):
        # Initilize Data
        T = np.array([[self.t_st]])
        Y = np.expand_dims(y0,0)

        dt = self.dt_min
        t = self.t_st + dt
        while t <= self.t_end:
            # Update arguments
            args = self.update_args(t, *args)
            # Find the solution
            y = self.step(y0, *args)
            # Collect data
            T = np.vstack((X, np.expand_dims(t,0)))
            Y = np.vstack((Y, np.expand_dims(y,0)))
            # Update step
            dt = min(dt*self.dt_st, self.dt_max)
            t = t + dt
            y0 = y

        return T, Y

    def step(self, y0, *args):
        if not isinstance(args, tuple):
            args = tuple(args)
        sol = optimize.root(
            self.fun, y0, args=args, method=self.method, jac=self.jac
        )
        return sol.x

    def update_args(self, t, *args):
        pass
