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

    # Function ----------------------------------------------------------------
    # def function(self, t, y, mass, T):
    #     # Call update method
    #     # self.update(y, T, mass)
    #     # Evaluate contributions from reactions
    #     drdt = self.wr(y)
    #     return drdt

    def function(self, t, y, mass):
        # Split variables into `rho` vector and `T`
        r, T = np.split(y, self.varIndices)
        # Call update method
        self.update(r, T, mass)
        # Evaluate contributions from reactions
        drdt = self.wr(r)
        # Evaluate the effect on the thermodynamic system
        dTdt = self.wT(drdt, T, mass)
        return np.concatenate(tuple([drdt, dTdt]))

    def wr(self, r):
        # Get Master Equation matrices
        Ke, Kd, Kr = self.chemModel.K
        # Evaluate contributions from each process
        wr_ = np.matmul(Ke, r) * r[-1]
        wr_ = wr_ + np.matmul(Kd, r) * r[-1]
        wr_ = wr_ + np.squeeze(Kr) * r[-1]**3
        return wr_

    def wT(self, drdt, T, mass):
        # Get mass fractions derivative
        dYdt = {
            name: np.take(drdt, idx)/mass \
                for name, idx in self.chemModel.specieIndices.items()
        }
        # Evaluate source term
        wT_ = - self.dehdY(T, dYdt) / self.cvp(T)
        return utils.convert_to_array(wT_)

    # Jacobian ----------------------------------------------------------------
    # @tf.function
    # def jacobian(self, t, y, mass):
    #     y = tf.constant(y)
    #     with tf.GradientTape() as g:
    #         g.watch(y)
    #         w = self.function(t, y, mass)
    #     dwdy = g.jacobian(w, y)
    #     return dwdy.numpy()


    # def jacobian(self, t, y, mass, T):
    #     return self.dwrdr(y)


    def dwrdr(self, r):
        # Get Master Equation matrices
        Ke, Kd, Kr = self.chemModel.K
        # Evaluate sources
        I = np.array([[0]*self.chemModel.spTh['O2'].specie.n_bins+[1]])
        wr_ = np.expand_dims(self.wr(r), 1)
        dwrdr_ = r[-1] * (Ke + Kd)
        dwrdr_ = dwrdr_ + np.matmul(wr_ + 2*r[-1]**3*Kr, I) / r[-1]
        return dwrdr_

    def dwrdT(self, r):
        # Get Master Equation matrices
        dKedT, dKddT, dKrdT = self.chemModel.dKdT
        # Evaluate contributions from each process
        dwrdT_ = np.matmul(dKedT, r) * r[-1]
        dwrdT_ = dwrdT_ + np.matmul(dKddT, r) * r[-1]
        dwrdT_ = dwrdT_ + np.squeeze(dKrdT) * r[-1]**3
        return dwrdT_

    def dwTdr(self, r, drdt, T, mass):
        # Get needed thermo quantities
        eh_i_ = np.expand_dims(self.eh_i(T), 0)
        cvp_i_ = np.expand_dims(self.cvp_i(T), 0)
        # Get needed source terms
        dwrdr_ = self.dwrdr(r)
        wT_ = self.wT(drdt, T, mass)
        # Evaluate source term
        dwTdr_ = - np.matmul(eh_i_, dwrdr_)
        dwTdr_ = dwTdr_ + wT_ * cvp_i_
        dwTdr_ = dwTdr_ / ( self.cvp(T) * mass )
        return dwTdr_

    def dwTdT(self, drdt, T, mass):
        dwTdr_ = 0.
        return dwTdT_
