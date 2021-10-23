import abc


class Basic(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        specieThermos,
        *args,
        **kwargs
    ):
        self.spTh = specieThermos

    # Methods
    ###########################################################################
    # Mixture mass fractions --------------------------------------------------
    @abc.abstractmethod
    def update(self, mixture, var='Y', *args, **kwargs):
        pass

    # Enthalpy ----------------------------------------------------------------
    @abc.abstractmethod
    def cp(self, T, *args, **kwargs):
        # [J/(kg K)]
        pass

    @abc.abstractmethod
    def h(self, T, *args, **kwargs):
        # [J/kg]
        pass

    @abc.abstractmethod
    def dhdY(self, T, dY, *args, **kwargs):
        # [J/kg]
        pass

    # Energy ------------------------------------------------------------------
    @abc.abstractmethod
    def cv(self, T, *args, **kwargs):
        # [J/(kg K)]
        pass

    @abc.abstractmethod
    def e(self, T, *args, **kwargs):
        # [J/kg]
        pass

    @abc.abstractmethod
    def dedY(self, T, dY, *args, **kwargs):
        # [J/kg]
        pass
