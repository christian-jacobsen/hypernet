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
    # Update ------------------------------------------------------------------
    def update(self, T):
        self.p = self.p_(T)
        self.psi = self.psi_(T)

    @abc.abstractmethod
    def p_(self, T):
        pass

    @abc.abstractmethod
    def psi_(self, T):
        pass
