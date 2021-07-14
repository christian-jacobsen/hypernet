import numpy as np

from hypernet.src.general import const
from hypernet.src.general import utils
from hypernet.src.thermophysicalModels.specie import specie as specie_module
from hypernet.src.thermophysicalModels.specie import thermo as thermo_module


class MultiComponent(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        mixture,
        thermo,
        *args,
        **kwargs
    ):
        self.mixture = {}
        for name in mixture:
            self.mixture[name] = utils.get_class(thermo_module, thermo)(
                specie_module.Specie(name, *args, **kwargs)
            )

        
    # Methods
    ###########################################################################
    # Mixture mass fractions ==================================================
    def update(self, mixture, var='Y'):
        for name, value in mixture.items():
            setattr(self.mixture[name].specie, var, utils.convert_to_array(value))
        if var == 'Y':
            self.X_i()
        elif var == 'X':
            self.Y_i()

    # Enthalpy ================================================================
    def cp(self,T):
        # [J/(kg K)]
        return self.mix_quantity('cp',T)

    def h(self,T):
        # [J/kg]
        return self.mix_quantity('h',T)

    # Internal Energy =========================================================
    def cv(self,T):
        # [J/(kg K)]
        return self.mix_quantity('cv',T)

    def e(self,T):
        # [J/kg]
        return self.mix_quantity('e',T)

    # Energy of formation -----------------------------------------------------
    def e_f(self):
        # [J/kg]
        return self.mix_quantity('e_f',T)

    # Translational Internal Energy -------------------------------------------
    def cv_tr(self):
        # [J/(kg K)]
        return self.mix_quantity('cv_tr',T)

    def e_tr(self,T):
        # [J/kg]
        return self.mix_quantity('e_tr',T)

    # Ro-Vibrational Internal Energy ------------------------------------------
    def cv_rv(self,T):
        # [J/(kg K)]
        return self.mix_quantity('cv_rv',T)

    def e_rv(self,T):
        # [J/kg]
        return self.mix_quantity('e_rv',T)

    # Mixture properties ======================================================
    def m(self, var='Y'):
        # [kg/mol]
        m_ = []
        if var == 'Y':
            for specie, thermo in self.mixture.items():
                m_.append(thermo.specie.Y / thermo.specie.m)
            return 1./np.sum(np.concatenate(tuple(m_)))
        elif var == 'X':
            for specie, thermo in self.mixture.items():
                m_.append(thermo.specie.X * thermo.specie.m)
            return np.sum(np.concatenate(tuple(m_)))

    def R(self):
        # [J/(kg K)]
        return self.mix_quantity(const.URG)

    def p(self, rho, T):
        return rho*self.R()*T

    def rho(self, p, T):
        return p/(self.R()*T)

    # Species properties ======================================================
    def Y_i(self):
        for specie, thermo in self.mixture.items():
            thermo.specie.Y = np.clip(
                thermo.specie.X * thermo.specie.m / self.m(var='X'), a_min=0., a_max=1.
            )

    def X_i(self):
        for specie, thermo in self.mixture.items():
            thermo.specie.X = np.clip(
                thermo.specie.Y * self.m(var='Y') / thermo.specie.m, a_min=0., a_max=1.
            )

    def n_i(self, rho=None, n=None, var='Y'):
        for specie, thermo in self.mixture.items():
            if var == 'Y':
                thermo.specie.n = np.clip(
                    thermo.specie.Y * rho / thermo.specie.m * const.UNA, a_min=0., a_max=None
                )
            elif var == 'X':
                thermo.specie.n = np.clip(
                    thermo.specie.X * n, a_min=0.
                )

    def rho_i(self, rho=None, n=None, var='Y'):
        for specie, thermo in self.mixture.items():
            if var == 'Y':
                thermo.specie.rho = np.clip(
                    thermo.specie.Y * rho, a_min=0., a_max=None
                )
            elif var == 'X':
                thermo.specie.rho = np.clip(
                    thermo.specie.X * n * thermo.specie.m / const.UNA, a_min=0., a_max=None
                )

    # Utils ===================================================================
    def mix_quantity(self, quantity, *args):
        var = []
        for specie, thermo in self.mixture.items():
            if isinstance(quantity, str):
                var_ = getattr(thermo, quantity)(*args)
            else:
                var_ = quantity
            var_ = utils.convert_to_array(var_)
            var.append(thermo.specie.Y / thermo.specie.m * var_)
        return np.sum(np.concatenate(tuple(var)))
