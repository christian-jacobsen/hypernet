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
            name: spTh_.intPF.Q * spTh_.transPF.Q \
                for name, spTh_ in self.spTh.items()
        }

    def dPFdotdT_(self):
        return {
            name: np.sum(spTh_.intPF.dQdT * spTh_.transPF.Q \
                + spTh_.intPF.Q * spTh_.transPF.dQdT) \
                for name, spTh_ in self.spTh.items()
        }

    # Dissociation equilibrium constants --------------------------------------
    def Keq_(self):
        return self.PFdot['O']**2 / self.PFdot['O2']

    def dKeqdT_(self):
        dKeqdT = 2 * self.PFdot['O'] * self.dPFdotdT['O']
        dKeqdT = dKeqdT - self.Keq * self.dPFdotdT['O2']
        dKeqdT = dKeqdT / self.dPFdotdT['O2']
        return dKeqdT

    # Reverse reaction rates --------------------------------------------------
    def kr_(self, kf, reacIndex, indices=None):
        if reacIndex in self.processIndices['diss']:
            i = indices[0]
            Keq_ = self.Keq[i]
        else:
            l, r = indices
            Q_ = self.spTh['O2'].intPF.Q
            Keq_ = Q_[r] / Q_[l]
        return kf / Keq_

    def dkrdT_(self, kr, dkfdT, reacIndex, indices=None):
        if reacIndex in self.processIndices['diss']:
            i = indices[0]
            dkrdT_ = dkfdT / self.Keq[i]
            dkrdT_ = dkrdT_ - dkfdT / self.Keq[i]**2 * self.dKeqdT[i]
            return dkrdT_
        else:
            l, r = indices
            Q_ = self.spTh['O2'].intPF.Q
            dQdT_ = self.spTh['O2'].intPF.dQdT
            dkrdT_ = dkfdT * Q_[l] / Q_[r]
            dkrdT_ = dkrdT_ + kr * ( dQdT_[l] / Q_[l] - dQdT_[r] / Q_[r] )
            return dkrdT_
