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

        self.Q = 1.
        self.dQdT = 0.
    
    # Methods
    ###########################################################################
    # Update ------------------------------------------------------------------
    def update(self, T):
        self.Q = self.Q_(T)
        self.dQdT = self.dQdT_(T)

    # Partition functions -----------------------------------------------------
    @abc.abstractmethod
    def Q_(self, T):
        pass

    @abc.abstractmethod
    def dQdT_(self, T):
        pass
