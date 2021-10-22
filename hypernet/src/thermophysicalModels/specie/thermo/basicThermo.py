import abc


class Thermo(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        specie,
        *args,
        **kwargs
    ):
        # Specie Properties ===================================================
        self.specie = specie

    # Methods
    ###########################################################################
    # Enthalpy ================================================================
    @abc.abstractmethod
    def cp(self, T):
        # [J/(mol K)]
        pass

    @abc.abstractmethod
    def h(self, T):
        # [J/mol]
        pass

    # Energy ==================================================================
    @abc.abstractmethod
    def cv(self, T):
        # [J/(mol K)]
        pass

    @abc.abstractmethod
    def e(self, T):
        # [J/mol]
        pass

    # Energy of formation -----------------------------------------------------
    def e_f(self):
        # [J/mol]
        return self.specie.Ef