import abc


class Basic(object):

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
        
    # Methods
    ###########################################################################
    @abc.abstractmethod
    def solve(self, y0, *args):
        pass