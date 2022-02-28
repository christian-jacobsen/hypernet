import abc


class Basic(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        specie,
        constVP='V',
        *args,
        **kwargs
    ):
        # Specie Properties
        self.specie = specie

        # Constant pressure/volume process
        self.constVP = constVP
        if self.constVP == 'P':
            self.he_ = self.h_
            self.cpv_ = self.cp_
        elif self.constVP == 'V':
            self.cpv_ = self.cv_
            self.he_ = self.e_
        self.dcpvdT_ = self.dcvdT_

    # Methods
    ###########################################################################
    # Update ------------------------------------------------------------------
    def update(self, T):
        self.he = self.he_(T)
        self.cpv = self.cpv_(T)
        self.dcpvdT = self.dcpvdT_(T)

    # Enthalpy ----------------------------------------------------------------
    def h_(self, T):
        # [J/kg]
        return self.e_(T) + self.specie.R * T

    def cp_(self):
        # [J/(kg K)]
        return self.cv_(T) + self.specie.R

    # Energy ------------------------------------------------------------------
    @abc.abstractmethod
    def e_(self, T):
        # [J/kg]
        pass

    @abc.abstractmethod
    def cv_(self, T):
        # [J/(kg K)]
        pass

    @abc.abstractmethod
    def dcvdT_(self, T):
        # [J/(kg K^2)]
        pass

    # Energy of formation -----------------------------------------------------
    def e_f_(self):
        # [J/kg]
        return self.specie.Ef / self.specie.M
