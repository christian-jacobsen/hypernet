import numpy as np
import tensorflow as tf

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
    # def update(self, t, y, mass):
    def update(self, r, T, mass):
        # r, T = np.split(y, self.varIndeces)
        # Update chemistry model
        self.chemModel.update(T)
        # Update mixture
        self.mixture.update(
            {name: np.take(r, idx)/mass for idx in self.chemModel.specieIndeces}
        )

    # Function ----------------------------------------------------------------
    def function(self, t, y, mass):
        # Split variables into `rho` vector and `T`
        r, T = np.split(y, self.varIndeces)
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
        wr_ += np.matmul(Kd, r) * r[-1]
        wr_ += np.squeeze(Kr) * r[-1]**3
        return wr_

    def wT(self, drdt, T, mass):
        # Get mass fractions derivative
        dYdt = {name: np.take(drdt, idx)/mass for idx in self.specieIndeces}
        # Evaluate source term
        wT_ = - self.dehdY(T, dYdt) / self.cvp(T)
        return wT_

    # Jacobian ----------------------------------------------------------------
    # @tf.function
    def jacobian(self, t, y, mass):
        with tf.GradientTape() as g:
            g.watch(y)
            w = self.function(y, mass)
        dwdy = g.jacobian(w, y)
        return dwdy.numpy()


    # def jacobian(self, r, T, mass):
    #     Y, T = y[:-1], y[-1]

    #     # Evaluate contributions from reactions
    #     dYdt = self.omega(Y, T)

    #     # Evaluate the effect on the thermodynamic system
    #     # >> Update mixture
    #     self.mixture.update(
    #         {name: np.take(Y, idx) for idx in self.specieIndeces}
    #     )
    #     # >> Evaluate mixture Cv
    #     cv_ = self.cvp(T)
    #     # >> Evaluate dTdt
    #     dYdt_ = {name: np.take(omega_, idx) for idx in self.specieIndeces}
    #     dedY_ = self.mixture.dedY(T, dYdt_)
    #     dTdt = - dedY_ / self.cvp(T)
    #     return excit + diss + recom


    # def dwrdr(self, r):
    #     # Get Master Equation matrices
    #     Ke, Kd, Kr = self.chemModel.K

    #     # Get Master Equation matrices
    #     I = np.array([[0]*self.chemModel.spTh['O2'].n_bins+[1]])
    #     wr_ = np.expand_dims(self.wr(r), 1)

    #     dwrdr_ = r[-1] * (Ke + Kd)
    #     dwrdr_ += np.matmul(wr_ + 2*r[-1]**3*Kr, I) / r[-1]

    #     return dwrdr_

    # def dwrdT(self, r):

    #     # Get Master Equation matrices
    #     dKedT, dKddT, dKrdT = self.chemModel.dKdT

    #     # Evaluate contributions from each process
    #     dwrdT_ = np.matmul(dKedT, r) * r[-1]
    #     dwrdT_ += np.matmul(dKddT, r) * r[-1]
    #     dwrdT_ += np.squeeze(dKrdT) * r[-1]**3

    #     return dwrdT_



    # def domegaYdY(self, Y, T):
    #     # Get Master Equation matrices
    #     K_e, K_d, K_r = self.chemModel.K

    # def jac(self, t, y, arg):
    #     '''Jacobian calculation: jac[i,j] = df[i] / dy[j].'''
    #     d = np.shape(y)[0]
    #     J = np.zeros((d,d), dtype=np.float64)
    #     J[:,:-1] = (arg[0][:,:-1] + arg[1][:,:-1]) * y[-1]
    #     J[:, -1] = np.matmul(arg[0][:,:-1], y[:-1]) \
    #         + np.matmul(arg[1][:,:-1], y[:-1]) \
    #         + 3*np.squeeze(arg[2])*np.power(y[-1], 2)
    #     return J



    #     # Evaluate contributions from each process
    #     excit = np.matmul(K_e, Y) * Y[-1]
    #     diss = np.matmul(K_d, Y) * Y[-1]
    #     recom = np.squeeze(K_r) * np.power(Y[-1], 3)

    #     return excit + diss + recom

    # def domegaYdT(self, Y, dYdt, T):
    #     dYdt = {name: np.take(dYdt, idx) for idx in self.specieIndeces}
    #     return - self.dehdY(T, dYdt) / self.cvp(T)

    # def domegaTdY(self, Y, T):
    #     # Get Master Equation matrices
    #     K_e, K_d, K_r = self.matrices(T)

    #     # Evaluate contributions from each process
    #     excit = np.matmul(K_e, Y) * Y[-1]
    #     diss = np.matmul(K_d, Y) * Y[-1]
    #     recom = np.squeeze(K_r) * np.power(Y[-1], 3)

    #     return excit + diss + recom

    # def domegaTdT(self, Y, dYdt, T):
    #     dYdt = {name: np.take(dYdt, idx) for idx in self.specieIndeces}
    #     return - self.dehdY(T, dYdt) / self.cvp(T)
