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

    def cp_i(self, T):
        # [J/(kg K)]
        _cp_i = self.quantity_i('cp', T)
        _cp_i = np.concatenate(tuple(_cp_i))
        return _cp_i

    def dcp_idT(self, T):
        # [J/(kg K)]
        _dcp_idT = self.quantity_i('dcpdT', T)
        _dcp_idT = np.concatenate(tuple(_dcp_idT))
        return _dcp_idT

    def h(self, T):
        # [J/kg]
        return self.mix_quantity('h', T)
        
    def h_i(self, T):
        # [J/(kg K)]
        _h_i = self.quantity_i('h', T)
        _h_i = np.concatenate(tuple(_h_i))
        return _h_i

    def dhdY(self, T, dY):
        # [J/kg]
        dhdY_ = 0.
        for name, spTh_ in self.spTh.items():
            dhdY_ = dhdY_ + np.sum(
                dY[name] / spTh_.specie.M * spTh_.thermo.h(T)
            )
        return dhdY_

    # Energy ------------------------------------------------------------------
    def cv(self, T):
        # [J/(kg K)]
        return self.mix_quantity('cv', T)

    def cv_i(self, T):
        # [J/(kg K)]
        _cv_i = self.quantity_i('cv', T)
        _cv_i = np.concatenate(tuple(_cv_i))
        return _cv_i

    def dcv_idT(self, T):
        # [J/(kg K)]
        _dcv_idT = self.quantity_i('dcvdT', T)
        _dcv_idT = np.concatenate(tuple(_dcv_idT))
        return _dcv_idT

    def e(self, T):
        # [J/kg]
        return self.mix_quantity('e', T)
        
    def e_i(self, T):
        # [J/(kg K)]
        _e_i = self.quantity_i('e', T)
        _e_i = np.concatenate(tuple(_e_i))
        return _e_i

    def dedY(self, T, dY):
        # [J/kg]
        dedY_ = 0.
        for name, spTh_ in self.spTh.items():
            dedY_ = dedY_ + np.sum(
                dY[name] / spTh_.specie.M * spTh_.thermo.e(T)
            )
        return dedY_

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
        _R = 0.
        for spTh_ in self.spTh.values():
            _R = _R + np.sum(
                spTh_.specie.Y / spTh_.specie.M * const.URG
            )
        self.R = _R

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
        var_i = self.quantity_i(quantity, *args)
        for i, spTh_ in enumerate(self.spTh.values()):
            var.append(spTh_.specie.Y * var_i[i])
        return np.sum(np.concatenate(tuple(var)))

    def quantity_i(self, quantity, *args):
        var = []
        for spTh_ in self.spTh.values():
            var_ = getattr(spTh_.thermo, quantity)(*args)
            var_ = utils.convert_to_array(var_)
            var.append(var_ / spTh_.specie.M)
        return var
