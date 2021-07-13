import numpy as np
import hypernet.src.general.const as const
from hypernet.src.thermophysicalModels.specie import Specie
from hypernet.src.thermophysicalModels.specie import Thermo

from scipy import optimize

class Root(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        setup
    ):
        # Start & end points
        self.x_st = setup['x']['start']
        self.x_end = setup['x']['end']

        # Marching setup
        self.dx_min = setup['dx']['min']
        self.dx_max = setup['dx']['max']
        self.dx_str = setup['dx']['str']

        # Root finding method
        self.method = setup['method'] if 'method' in setup else 'hybr'
        
    # Methods
    ###########################################################################
    @utils.timing()
    def solve(self, y0):
        # Initilize Data
        X = np.array([[self.x_st]])
        Y = np.expand_dims(y0,0)

        dx = self.dx_min
        x = self.x_st + dx
        while x <= self.x_end:
            # Update arguments
            args = self.get_args(x)
            # Find the solution
            y = self.step(y0, args)
            # Collect data
            X = np.vstack((X, np.expand_dims(x,0)))
            Y = np.vstack((Y, np.expand_dims(y,0)))
            # Update step
            dx = min(dx*self.dx_str, self.dx_max)
            x = x + dx
            y0 = y

        return X, Y

    def step(self, y0, args):
        sol = optimize.root(self.fun, y0, args=tuple(args), \
            method=self.method, jac=self.jac)
        return sol.x

    def fun(self, y, args):
        return None

    def jac(self, y, args):
        return None

    def get_args(self, x):
        return None
