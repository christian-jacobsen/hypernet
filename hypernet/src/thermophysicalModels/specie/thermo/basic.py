import abc

from hypernet.src.general import const


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
    # Enthalpy ----------------------------------------------------------------
    @abc.abstractmethod
    def cp(self, T):
        # [J/(mol K)]
        return self.cv(T) + const.URG

    @abc.abstractmethod
    def h(self, T):
        # [J/mol]
        return self.e(T) + const.URG*T

    # Energy ------------------------------------------------------------------
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