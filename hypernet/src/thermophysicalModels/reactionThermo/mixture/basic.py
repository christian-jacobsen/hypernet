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
        # Thermodynamic specie properties
        self.spTh = specieThermos

    # Properties
    ###########################################################################
    @property
    def M(self):
        return self._M
    @M.setter
    def M(self, value):
        self._M = value

    @property
    def R(self):
        return self._R
    @R.setter
    def R(self, value):
        self._R = value

    # Methods
    ###########################################################################
    # Update ------------------------------------------------------------------
    @abc.abstractmethod
    def update(self, *args, **kwargs):
        pass

    # Enthalpy/Energy ---------------------------------------------------------
    @abc.abstractmethod
    def he_(self, *args, **kwargs):
        # [J/kg]
        pass

    @abc.abstractmethod
    def cpv_(self, *args, **kwargs):
        # [J/(kg K)]
        pass

    @abc.abstractmethod
    def dhedY_(self, *args, **kwargs):
        # [J/kg]
        pass

    # Mixture properties ------------------------------------------------------
    @abc.abstractmethod
    def M_(self, *args, **kwargs):
        # [kg/mol]
        pass

    @abc.abstractmethod
    def R_(self, *args, **kwargs):
        # [J/(kg K)]
        pass
