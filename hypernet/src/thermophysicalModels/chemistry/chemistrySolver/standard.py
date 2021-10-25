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
        constPV='V',
        *args,
        **kwargs
    ):
        super(Standard, self).__init__(
            mixture,
            specieThermos,
            chemistryModel,
            reactionsList=reactionsList,
            processFlags=processFlags,
            constPV=constPV,
            *args,
            **kwargs
        )

    # Methods
    ###########################################################################
    # Update method -----------------------------------------------------------
    def update(self, r, T, mass):
        # Update chemistry model
        self.chemModel.update(T)
        # Update mixture
        self.mixture.update(
            {name: np.take(r, idx)/mass \
                for name, idx in self.chemModel.specieIndices.items()}
        )

    def w(self, r, T, mass):
        # Get first order source terms
        self.wr = self.wr_(r)
        self.wT = self.wT_(self.wr, T, mass)

    def dw(self, r, T, mass):
        # Get second order source terms
        self.dwrdr = self.dwrdr_(r, self.wr)
        self.dwrdT = self.dwrdT_(r)
        self.dwTdr = self.dwTdr_(r, self.dwrdr, T, self.wT, mass)
        self.dwTdT = self.dwTdT_(r, self.wr, T, self.wT, mass)

    # Function ----------------------------------------------------------------
    # def function(self, t, y, mass, T):
    #     # Call update method
    #     # self.update(y, T, mass)
    #     # Evaluate contributions from reactions
    #     drdt = self.wr_(y)
    #     return drdt

    def function(self, t, y, mass):
        # Split variables into `rho` vector and `T`
        r, T = np.split(y, self.varIndices)
        # Call update method
        self.update(r, T, mass)
        self.w(r, T, mass)
        return np.concatenate(tuple([self.wr, self.wT]))

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
    def wT_(self, wr, T, mass):
        # Get mass fractions derivative
        dYdt = {
            name: np.take(wr, idx)/mass \
                for name, idx in self.chemModel.specieIndices.items()
        }
        # Evaluate source term
        wT = - self.dehdY(T, dYdt) / self.cvp(T)
        return utils.convert_to_array(wT)

    # Jacobian ----------------------------------------------------------------
    def jacobian(self, t, y, mass):
        # Split variables into `rho` vector and `T`
        r, T = np.split(y, self.varIndices)
        # Call update method
        self.update(r, T, mass)
        self.w(r, T, mass)
        self.dw(r, T, mass)
        # Evaluate contributions from reactions
        dwr = np.concatenate([self.dwrdr, self.dwrdT], axis=1)
        # Evaluate the effect on the thermodynamic system
        dwT = np.concatenate([self.dwTdr, self.dwTdT])
        dwT = np.expand_dims(dwT, 0)
        return np.concatenate([dwr, dwT])

    # Source terms from reactions
    def dwrdr_(self, r, wr):
        # Get Master Equation matrices
        Ke, Kd, Kr = self.chemModel.K
        # Evaluate sources
        I = np.array([[0]*(self.chemModel.nSpecies-1)+[1]])
        wr = np.expand_dims(wr, 1) if len(wr.shape) == 1 else wr
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
        return np.expand_dims(dwrdT, 1)

    # Source terms from thermodynamics
    def dwTdr_(self, r, dwrdr, T, wT, mass):
        # Get needed thermo quantities
        eh_i_ = self.eh_i(T)
        cvp_i_ = self.cvp_i(T)
        cvp_ = self.cvp(T)
        # Evaluate source term
        dwTdr = - np.matmul(eh_i_, dwrdr)
        dwTdr = dwTdr + wT * cvp_i_
        dwTdr = dwTdr / ( cvp_ * mass )
        return dwTdr

    def dwTdT_(self, r, wr, T, wT, mass):
        # Get needed thermo quantities
        cvp_i_ = self.cvp_i(T)
        dcvp_idT_ = self.dcvp_idT(T)
        cvp_ = self.cvp(T)
        # Evaluate source term
        dwTdT = - np.sum(cvp_i_ * wr, keepdims=True)
        dwTdT = dwTdT + wT * np.sum(dcvp_idT_ * r, keepdims=True)
        dwTdT = dwTdT / ( cvp_ * mass )
        return dwTdT

    # @tf.function
    # def jacobian(self, t, y, mass):
    #     y = tf.constant(y)
    #     with tf.GradientTape() as g:
    #         g.watch(y)
    #         w = self.function(t, y, mass)
    #     dwdy = g.jacobian(w, y)
    #     return dwdy.numpy()