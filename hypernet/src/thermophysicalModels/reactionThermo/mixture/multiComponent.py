import numpy as np

from hypernet.src.general import const
from hypernet.src.general import utils
from hypernet.src.thermophysicalModels.reactionThermo.mixture import Basic


class MultiComponent(Basic):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        specieThermos,
        *args,
        **kwargs
    ):
        super(MultiComponent, self).__init__(specieThermos)

    # Methods
    ###########################################################################
    # Mixture properties ------------------------------------------------------
    def update(self, XY, var='Y'):
        # Update mass/molar fractions
        for name, value in XY.items():
            value = utils.check_XY(utils.convert_to_array(value))
            setattr(self.spTh[name].specie, var, value)
        # Update mixture/species properties
        self.M = self.M_(var=var)
        if var == 'Y':
            self.Xi_()
        elif var == 'X':
            self.Yi_()
        self.R = self.R_()

    # Mixture properties ------------------------------------------------------
    # Mass
    def M_(self, var='Y'):
        # [kg/mol]
        if var == 'Y':
            M = [spTh.specie.Y / spTh.specie.M for spTh in self.spTh.values()]
            return 1./np.sum(np.concatenate(M))
        elif var == 'X':
            M = [spTh.specie.X * spTh.specie.M for spTh in self.spTh.values()]
            return np.sum(np.concatenate(M))

    # Specific gas constant
    def R_(self):
        R = [spTh.specie.Y * spTh.specie.R for spTh in self.spTh.values()]
        return np.sum(np.concatenate(R))

    # Pressure
    def p_(self, rho, T):
        return rho*self.R*T

    # Density
    def rho_(self, p, T):
        return p/(self.R*T)

    # Number density
    def n_(self, rho):
        self.ni_(rho=rho, var='Y')
        n = [spTh.specie.n for spTh in self.spTh.values()]
        return np.sum(np.concatenate(n))

    # Enthalpy/Energy
    def he_(self):
        # [J/kg]
        he = [spTh.specie.Y * spTh.thermo.he for spTh in self.spTh.values()]
        return np.sum(np.concatenate(he))

    def cpv_(self):
        # [J/(kg K)]
        cpv = [spTh.specie.Y * spTh.thermo.cpv for spTh in self.spTh.values()]
        return np.sum(np.concatenate(cpv))

    def dcpvdT_(self):
        # [J/kg]
        dcpvdT = [
            spTh.specie.Y * spTh.thermo.dcpvdT for spTh in self.spTh.values()
        ]
        return np.sum(np.concatenate(dcpvdT))

    def dhedY_(self, dY):
        # [J/kg]
        dhedY = [
            np.sum(dY[name] * spTh.thermo.he) \
                for name, spTh in self.spTh.items()
        ]
        return np.sum(dhedY)

    # Species properties ------------------------------------------------------
    def Yi_(self):
        for spTh_ in self.spTh.values():
            sp = spTh_.specie
            sp.Y = sp.X * sp.M / self.M

    def Xi_(self):
        for spTh_ in self.spTh.values():
            sp = spTh_.specie
            sp.X = sp.Y * self.M / sp.M

    def ni_(self, rho=None, n=None, var='Y'):
        for spTh_ in self.spTh.values():
            sp = spTh_.specie
            if var == 'Y':
                sp.n = sp.Y * rho / sp.M * const.UNA
            elif var == 'X':
                sp.n = sp.X * n

    def rhoi_(self, rho=None, n=None, var='Y'):
        for spTh_ in self.spTh.values():
            sp = spTh_.specie
            if var == 'Y':
                sp.rho = sp.Y * rho
            elif var == 'X':
                sp.rho = sp.X * n * sp.M / const.UNA
