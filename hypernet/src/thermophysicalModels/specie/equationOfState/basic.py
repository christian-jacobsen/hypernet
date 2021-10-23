import abc


class Basic(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        specie,
        *args,
        **kwargs
    ):
        # Specie Properties
        self.specie = specie

    # Methods
    ###########################################################################
    @abc.abstractmethod
    def p(self, T):
        pass

    @abc.abstractmethod
    def psi(self, T):
        pass
