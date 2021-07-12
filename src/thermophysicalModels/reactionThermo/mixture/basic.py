import numpy as np
import hypernet.src.general.const as const
from hypernet.src.thermophysicalModels.specie import Specie
from hypernet.src.thermophysicalModels.specie import Specie

class Basic(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        input_file,
        mixture
    ):
        self.mix = {}
        for s, Y in mixture.items():
            sp = Specie(input_file, s, Y)
            self.mix[s] = Thermo(sp)
        
    # Methods
    ###########################################################################
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

    def to_array(x):        
        if not isinstance(x, np.ndarray):
            x = np.array(x)
        return x

    def mix_quantity(self, quantity, args):
        var = []
        for s, th in self.mix.items():
            var.append(th.sp.Y * utils.to_array(getattr(th, quantity)(*args)))
        return np.sum(np.concatenate(tuple(var)))