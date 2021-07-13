import numpy as np
import hypernet.src.general.const as const
from hypernet.src.thermophysicalModels.specie import Specie
from hypernet.src.thermophysicalModels.specie import Thermo

class MultiComponent(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        input_file,
        mixture,
        *args,
        **kwargs
    ):
        self.mix = {}
        for specie in mixture:
            s = Specie(specie, args, kwargs)
            self.mix[specie] = Thermo(s)
        self.update(mixture)
        
    # Methods
    ###########################################################################

    def update(self, mixture):
        for specie, Y in mixture.items():
            self.mix[specie].sp.Y_(Y)

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

    def m(self):
        # [kg/mol]
        return 1./self.mix_quantity(1.)

    def R(self):
        # [J/(kg K)]
        return self.mix_quantity(const.URG)

    def p(self, rho, T):
        return self.mix_quantity(rho*const.URG*T)

    def X_i(self):
        for specie, thermo in self.mix.items():
            specie.X_(np.clip(thermo.sp.Y * self.m() / thermo.sp.m, a_min=0.))

    def n_i(self, rho):
        for specie, thermo in self.mix.items():
            specie.n_(np.clip(thermo.sp.Y * const.UNA * rho / thermo.sp.m, a_min=0.))

    def rho_i(self, rho):
        for specie, thermo in self.mix.items():
            specie.rho_(np.clip(thermo.sp.Y * rho, a_min=0.))

    def mix_quantity(self, quantity, args):
        var = []
        for specie, thermo in self.mix.items():
            if isinstance(quantity, str):
                var_ = getattr(thermo, quantity)(*args)
            else:
                var_ = quantity
            var_ = utils.convert_to_array(var_)
            var.append(thermo.sp.Y / thermo.sp.m * var_)
        return np.sum(np.concatenate(tuple(var)))
