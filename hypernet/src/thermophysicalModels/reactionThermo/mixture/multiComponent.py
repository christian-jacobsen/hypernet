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
    # Mixture mass fractions --------------------------------------------------
    def update(self, mixture, var='Y'):
        for name, value in mixture.items():
            setattr(self.spTh[name].specie, var, utils.convert_to_array(value))
        self.M_(var=var)
        if var == 'Y':
            self.X_i()
        elif var == 'X':
            self.Y_i()
        self.R_()

    # Enthalpy ----------------------------------------------------------------
    def cp(self, T):
        # [J/(kg K)]
        return self.mix_quantity('cp', T)

    def h(self, T):
        # [J/kg]
        return self.mix_quantity('h', T)

    def dhdY(self, T, dY):
        # [J/kg]
        dhdY_ = []
        for name, spTh_ in self.spTh.values():
            sp, th = spTh_.specie, spTh_.thermo
            dhdY_.append(dY[name] / sp.M * th.h(T))
        return np.sum(np.concatenate(tuple(dhdY_)))

    # Energy ------------------------------------------------------------------
    def cv(self, T):
        # [J/(kg K)]
        return self.mix_quantity('cv', T)

    def e(self, T):
        # [J/kg]
        return self.mix_quantity('e', T)

    def dedY(self, T, dY):
        # [J/kg]
        dedY_ = []
        for name, spTh_ in self.spTh.values():
            sp, th = spTh_.specie, spTh_.thermo
            dedY_.append(dY[name] / sp.M * th.e(T))
        return np.sum(np.concatenate(tuple(dedY_)))

    # Energy of formation
    def e_f(self):
        # [J/kg]
        return self.mix_quantity('e_f')

    # Translational Energy
    def cv_tr(self):
        # [J/(kg K)]
        return self.mix_quantity('cv_tr')

    def e_tr(self, T):
        # [J/kg]
        return self.mix_quantity('e_tr', T)

    # Ro-vibrational Energy
    def cv_int(self, T):
        # [J/(kg K)]
        return self.mix_quantity('cv_int', T)

    def e_int(self, T):
        # [J/kg]
        return self.mix_quantity('e_int', T)

    # Mixture properties ------------------------------------------------------
    def M_(self, var='Y'):
        # [kg/mol]
        M_ = []
        if var == 'Y':
            for spTh_ in self.spTh.values():
                M_.append(spTh_.specie.Y / spTh_.specie.M)
            self.M = 1./np.sum(np.concatenate(tuple(M_)))
        elif var == 'X':
            for spTh_ in self.spTh.values():
                M_.append(spTh_.specie.X * spTh_.specie.M)
            self.M = np.sum(np.concatenate(tuple(M_)))

    def R_(self):
        # [J/(kg K)]
        self.R = self.mix_quantity(const.URG)

    def p_(self, rho, T):
        return rho*self.R*T

    def rho_(self, p, T):
        return p/(self.R*T)

    # Species properties ------------------------------------------------------
    def Y_i(self):
        for spTh_ in self.spTh.values():
            sp = spTh_.specie
            sp.Y = sp.X * sp.M / self.M

    def X_i(self):
        for spTh_ in self.spTh.values():
            sp = spTh_.specie
            sp.X = sp.Y * self.M / sp.M

    def n_i(self, rho=None, n=None, var='Y'):
        for spTh_ in self.spTh.values():
            sp = spTh_.specie
            if var == 'Y':
                sp.n = sp.Y * rho / sp.M * const.UNA
            elif var == 'X':
                sp.n = sp.X * n

    def rho_i(self, rho=None, n=None, var='Y'):
        for spTh_ in self.spTh.values():
            sp = spTh_.specie
            if var == 'Y':
                sp.rho = sp.Y * rho
            elif var == 'X':
                sp.rho = sp.X * n * sp.M / const.UNA

    # Utils -------------------------------------------------------------------
    def mix_quantity(self, quantity, *args):
        var = []
        for spTh_ in self.spTh.values():
            sp, th = spTh_.specie, spTh_.thermo
            if isinstance(quantity, str):
                var_ = getattr(th, quantity)(*args)
            else:
                var_ = quantity
            var.append(sp.Y / sp.M * var_)
        return np.sum(np.concatenate(tuple(var)))
