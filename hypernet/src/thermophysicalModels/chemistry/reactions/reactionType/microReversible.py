import numpy as np

from hypernet.src.thermophysicalModels.chemistry.reactions.reactionType import Basic


class MicroReversible(Basic):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        specieThermos,
        reactionRate,
        processIndices,
        *args,
        **kwargs
    ):
        super(MicroReversible, self).__init__(
            specieThermos,
            reactionRate,
            processIndices,
            *args,
            **kwargs
        )

    # Properties
    ###########################################################################
    # Dissociation equilibrium constants --------------------------------------
    @property
    def Keq(self):
        return self._Keq
    @Keq.setter
    def Keq(self, value):
        self._Keq = value

    @property
    def dKeqdT(self):
        return self._dKeqdT
    @dKeqdT.setter
    def dKeqdT(self, value):
        self._dKeqdT = value

    # Partition functions products --------------------------------------------
    @property
    def PFdot(self):
        return self._PFdot
    @PFdot.setter
    def PFdot(self, value):
        self._PFdot = value

    @property
    def dPFdotdT(self):
        return self._dPFdotdT
    @dPFdotdT.setter
    def dPFdotdT(self, value):
        self._dPFdotdT = value

    # Methods
    ###########################################################################
    # Update method -----------------------------------------------------------
    def update(self, T):
        super(MicroReversible, self).update(T)
        # Evaluate partition functions products
        self.PFdot = self.PFdot_()
        self.dPFdotdT = self.dPFdotdT_()
        # Evaluate dissociation equilibrium constants
        self.Keq = self.Keq_()
        self.dKeqdT = self.dKeqdT_()

    # Partition functions products --------------------------------------------
    def PFdot_(self):
        return {
            name: spTh.intPF.Q * spTh.transPF.Q \
                for name, spTh in self.spTh.items()
        }

    def dPFdotdT_(self):
        return {
            name: np.sum(spTh.intPF.dQdT * spTh.transPF.Q \
                + spTh.intPF.Q * spTh.transPF.dQdT) \
                for name, spTh in self.spTh.items()
        }

    # Dissociation equilibrium constants --------------------------------------
    def Keq_(self):
        return self.PFdot[self.atom]**2 / self.PFdot[self.molecule]

    def dKeqdT_(self):
        dKeqdT = 2 * self.PFdot[self.atom] * self.dPFdotdT[self.atom]
        dKeqdT = dKeqdT - self.Keq * self.dPFdotdT[self.molecule]
        dKeqdT = dKeqdT / self.dPFdotdT[self.molecule]
        return dKeqdT

    # Reverse reaction rates --------------------------------------------------
    def kr_(self, kf, reacIndex, indices=None):
        if reacIndex in self.processIndices['diss']:
            i = indices[0]
            Keq = self.Keq[i]
        else:
            l, r = indices
            Q = self.spTh[self.molecule].intPF.Q
            Keq = Q[r] / Q[l]
        return kf / Keq

    def dkrdT_(self, kr, dkfdT, reacIndex, indices=None):
        if reacIndex in self.processIndices['diss']:
            i = indices[0]
            dkrdT = dkfdT / self.Keq[i]
            dkrdT = dkrdT - dkfdT / self.Keq[i]**2 * self.dKeqdT[i]
            return dkrdT
        else:
            l, r = indices
            Q = self.spTh[self.molecule].intPF.Q
            dQdT = self.spTh[self.molecule].intPF.dQdT
            dkrdT = dkfdT * Q[l] / Q[r]
            dkrdT = dkrdT + kr * ( dQdT[l] / Q[l] - dQdT[r] / Q[r] )
            return dkrdT
