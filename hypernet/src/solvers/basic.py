import abc


class Basic(object):

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
        # Start/end points
        self.t_st = setup['start']
        self.t_end = setup['end']

        # Marching setup
        self.dt_min = setup['delta']['min']
        self.dt_max = setup['delta']['max']
        self.dt_str = setup['delta']['str']

        # Integrand function/jacobian
        self.fun = function
        self.jac = jacobian

        # Function for updating function/jacobian arguments
        self.update_args = update_args
        
    # Methods
    ###########################################################################
    @abc.abstractmethod
    def solve(self, y0, *args):
        pass
