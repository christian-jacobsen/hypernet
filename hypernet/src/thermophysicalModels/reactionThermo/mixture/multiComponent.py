import numpy as np

from hypernet.src.general import const
from hypernet.src.general import utils


class MultiComponent(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        thermoSpecies,
        *args,
        **kwargs
    ):
        self.thSp = thermoSpecies

    # Methods
    ###########################################################################
    # Mixture mass fractions ==================================================
    def update(self, mixture, var='Y'):
        for name, value in mixture.items():
            setattr(self.thSp[name].specie, var, utils.convert_to_array(value))
        self.m_(var=var)
        if var == 'Y':
            self.X_i()
        elif var == 'X':
            self.Y_i()
        self.R_()

    # Enthalpy ================================================================
    def cp(self, T):
        # [J/(kg K)]
        return self.mix_quantity('cp', T)

    def h(self, T):
        # [J/kg]
        return self.mix_quantity('h', T)

    # Energy ==================================================================
    def cv(self, T):
        # [J/(kg K)]
        return self.mix_quantity('cv', T)

    def e(self, T):
        # [J/kg]
        return self.mix_quantity('e', T)

    # Energy of formation -----------------------------------------------------
    def e_f(self):
        # [J/kg]
        return self.mix_quantity('e_f')

    # Translational Energy ----------------------------------------------------
    def cv_tr(self):
        # [J/(kg K)]
        return self.mix_quantity('cv_tr')

    def e_tr(self, T):
        # [J/kg]
        return self.mix_quantity('e_tr', T)

    # Ro-vibrational Energy ---------------------------------------------------
    def cv_int(self, T):
        # [J/(kg K)]
        return self.mix_quantity('cv_int', T)

    def e_int(self, T):
        # [J/kg]
        return self.mix_quantity('e_int', T)

    # Mixture properties ======================================================
    def m_(self, var='Y'):
        # [kg/mol]
        M_ = []
        if var == 'Y':
            for thSp_ in self.thSp.values():
                M_.append(thSp_.specie.Y / thSp_.specie.M)
            self.M = 1./np.sum(np.concatenate(tuple(M_)))
        elif var == 'X':
            for thSp_ in self.thSp.values():
                M_.append(thSp_.specie.X * thSp_.specie.M)
            self.M = np.sum(np.concatenate(tuple(M_)))

    def R_(self):
        # [J/(kg K)]
        self.R = self.mix_quantity(const.URG)

    def p_(self, rho, T):
        return rho*self.R()*T

    def rho_(self, p, T):
        return p/(self.R()*T)

    # Species properties ======================================================
    def Y_i(self):
        for thSp_ in self.thSp.values():
            sp = thSp_.specie
            sp.Y = sp.X * sp.M / self.M

    def X_i(self):
        for thSp_ in self.thSp.values():
            sp = thSp_.specie
            sp.X = sp.Y * self.M / sp.M

    def n_i(self, rho=None, n=None, var='Y'):
        for thSp_ in self.thSp.values():
            sp = thSp_.specie
            if var == 'Y':
                sp.n = sp.Y * rho / sp.M * const.UNA
            elif var == 'X':
                sp.n = sp.X * n

    def rho_i(self, rho=None, n=None, var='Y'):
        for thSp_ in self.thSp.values():
            sp = thSp_.specie
            if var == 'Y':
                sp.rho = sp.Y * rho
            elif var == 'X':
                sp.rho = sp.X * n * sp.M / const.UNA

    # Utils ===================================================================
    def mix_quantity(self, quantity, *args):
        var = []
        for thSp_ in self.thSp.values():
            sp = thSp_.specie
            th = thSp_.thermo
            if isinstance(quantity, str):
                var_ = getattr(th, quantity)(*args)
            else:
                var_ = quantity
            var.append(sp.Y / sp.M * var_)
        return tf.reduce_sum(tf.concat(var, 1), axis=1, keepdims=True)
