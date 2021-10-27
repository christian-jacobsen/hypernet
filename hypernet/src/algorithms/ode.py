import numpy as np

from scipy.integrate import ode
from profilehooks import profile
from hypernet.src.general import utils
from hypernet.src.algorithms import Basic


class ODE(Basic):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        setup,
        function=None,
        jacobian=None,
        update_args=None,
        *args,
        **kwargs
    ):
        super(ODE, self).__init__(
            setup,
            function=function,
            jacobian=jacobian,
            update_args=update_args,
            *args,
            **kwargs
        )
        # Integration method
        self.method = setup['method'] if 'method' in setup else 'bdf'

    # Methods
    ###########################################################################
    def jac_norm(self, t, y, *args):
        J = self.jac(t, y, *args)
        norm = np.linalg.norm(J)
        return norm

    @profile
    @utils.timing
    def solve(self, y0, args=(), jac_norm=False):
        '''Solving the ODE.'''

        r = ode(self.fun, self.jac)
        r.set_integrator(
            'vode',
            method=self.method,
            with_jacobian=True if self.jac is not None else False,
            atol=1.e-20
        )
        r.set_initial_value(y0, self.t_st)
        r.set_f_params(*args)
        r.set_jac_params(*args)

        # Initilize Data
        T = np.array([[self.t_st]])
        Y = np.expand_dims(y0,0)
        # J_ = self.jac(self.t_st, y0, args)
        if jac_norm:
            J = self.jac_norm(self.t_st, y0, *args)

        # Integrate
        dt = self.dt_min
        while r.successful() and r.t <= self.t_end:
            r.integrate(r.t+dt)
            if jac_norm:
                J = np.vstack((J,self.jac_norm(self.t_st, y0, *args)))
            # Collect data
            T  = np.vstack((T, np.expand_dims(r.t,0)))
            Y  = np.vstack((Y, np.expand_dims(r.y,0)))
            print(Y.shape[0])
            # Update step
            dt = min(dt*self.dt_str, self.dt_max)
        assert (Y > 0.).all()

        # Return outputs
        out = [T, Y]
        if jac_norm:
            out.append(J)

        return out
