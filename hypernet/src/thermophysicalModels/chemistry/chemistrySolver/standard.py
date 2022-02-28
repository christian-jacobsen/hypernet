import numpy as np
import tensorflow as tf

from hypernet.src.general import utils
from hypernet.src.thermophysicalModels.chemistry.chemistrySolver import Basic


class Standard(Basic):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        mixture,
        specieThermos,
        chemistryModel,
        processFlags,
        reactionsList=None,
        heatBath='isothermal',
        *args,
        **kwargs
    ):
        super(Standard, self).__init__(
            mixture,
            specieThermos,
            chemistryModel,
            reactionsList=reactionsList,
            processFlags=processFlags,
            heatBath=heatBath,
            *args,
            **kwargs
        )

        self.update_chem_thermo = True

        self.function = self.function_
        if self.heatBath == 'isothermal':
            self.jacobian = self.jacobian_
        else:
            self.jacobian = None


    # Methods
    ###########################################################################
    # Update method -----------------------------------------------------------
    def update(self, r, T, rho):
        if self.update_chem_thermo:
            # Update thermodynamics
            for spTh in self.spTh.values():
                spTh.update(T)
            # Update chemistry model
            self.chemModel.update(T)
            if self.heatBath == 'isothermal':
                self.update_chem_thermo = False
        # Update mixture
        self.mix_update(r/rho)

    def mix_update(self, Y):
        # Update mixture
        self.mixture.update(
            {name: np.take(Y, idx) \
                for name, idx in self.chemModel.specieIndices.items()}
        )

    def w_(self, r, T, rho):
        # Get first order source terms
        self.wr = self.wr_(r)
        if self.heatBath == 'isothermal':
            self.wT = np.zeros(1)
        else:
            self.wT = self.wT_(self.wr, T, rho)
        self.w = np.concatenate(tuple([self.wr, self.wT]))

    def dw_(self, r, T, rho):
        # Get second order source terms
        self.dwrdr = self.dwrdr_(r, self.wr)
        if self.heatBath == 'isothermal':
            self.dwrdT = np.zeros((self.chemModel.nSpecies,1))
            self.dwTdr = np.zeros((1,self.chemModel.nSpecies))
            self.dwTdT = np.zeros((1,1))
        else:
            self.dwrdT = self.dwrdT_(r)
            self.dwTdr = self.dwTdr_(r, self.dwrdr, T, self.wT, rho)
            self.dwTdT = self.dwTdT_(r, self.wr, T, self.wT, rho)
        # Evaluate contributions from reactions
        dwr = np.concatenate([self.dwrdr, self.dwrdT], axis=1)
        # Evaluate the effect on the thermodynamic system
        dwT = np.concatenate([self.dwTdr, self.dwTdT], axis=1)
        self.dw = np.concatenate([dwr, dwT])

    # Function ----------------------------------------------------------------
    def function_(self, t, y, rho):
        # Split variables into `rho_i` vector and `T`
        r, T = np.split(y, [self.chemModel.nSpecies])
        # Call update method
        self.update(r, T, rho)
        self.w_(r, T, rho)
        return self.w

    # Source terms from reactions
    def wr_(self, r):
        # Get Master Equation matrices
        Ke, Kd, Kr = self.chemModel.K
        # Evaluate contributions from each process
        wr = np.matmul(Ke, r) * r[-1]
        wr = wr + np.matmul(Kd, r) * r[-1]
        wr = wr + np.squeeze(Kr) * r[-1]**3
        return wr

    # Source terms from thermodynamics
    def wT_(self, wr, T, rho):
        # Get rho fractions derivative
        dY = {
            name: np.take(wr, idx)/rho \
                for name, idx in self.chemModel.specieIndices.items()
        }
        # Evaluate source term
        wT = - self.mixture.dhedY_(dY) / self.mixture.cpv_()
        return utils.convert_to_array(wT)

    # Jacobian ----------------------------------------------------------------
    def jacobian_(self, t, y, rho):
        # Split variables into `rho_i` vector and `T`
        r, T = np.split(y, [self.chemModel.nSpecies])
        # Call update method
        self.update(r, T, rho)
        self.w_(r, T, rho)
        self.dw_(r, T, rho)
        return self.dw

    # Source terms from reactions
    def dwrdr_(self, r, wr):
        # Get Master Equation matrices
        Ke, Kd, Kr = self.chemModel.K
        # Evaluate sources
        I = np.array([[0]*(self.chemModel.nSpecies-1)+[1]])
        wr = wr.reshape(-1,1) if len(wr.shape) == 1 else wr
        dwrdr = r[-1] * (Ke + Kd)
        dwrdr = dwrdr + np.matmul(wr + 2*r[-1]**3*Kr, I) / r[-1]
        return dwrdr

    def dwrdT_(self, r):
        # Get Master Equation matrices
        dKedT, dKddT, dKrdT = self.chemModel.dKdT
        # Evaluate contributions from each process
        dwrdT = np.matmul(dKedT, r) * r[-1]
        dwrdT = dwrdT + np.matmul(dKddT, r) * r[-1]
        dwrdT = dwrdT + np.squeeze(dKrdT) * r[-1]**3
        return dwrdT.reshape(-1,1)

    # Source terms from thermodynamics
    def dwTdr_(self, r, dwrdr, T, wT, rho):
        # Get needed thermo quantities
        he_i = self.thermo_quantity('he')
        cpv_i = self.thermo_quantity('cpv')
        # Evaluate source term
        dwTdr = - np.matmul(he_i, dwrdr)
        dwTdr = dwTdr + wT * cpv_i
        dwTdr = dwTdr / ( self.mixture.cpv_() * rho )
        return dwTdr.reshape(1,-1)

    def dwTdT_(self, r, wr, T, wT, rho):
        # Get needed thermo quantities
        cpv_i = self.thermo_quantity('cpv')
        dcpv_idT = self.thermo_quantity('dcpvdT')
        # Evaluate source term
        dwTdT = - np.sum(cpv_i * wr, keepdims=True)
        dwTdT = dwTdT + wT * np.sum(dcpv_idT * r, keepdims=True)
        dwTdT = dwTdT / ( self.mixture.cpv_() * rho )
        return dwTdT.reshape(1,-1)

    # Utils -------------------------------------------------------------------
    def thermo_quantity(self, quantity, *args):
        var = []
        for spTh in self.spTh.values():
            var.append(
                utils.convert_to_array(getattr(spTh.thermo, quantity))
            )
        return np.concatenate(tuple(var))
