import numpy as np

from hypernet.src.thermophysicalModels.chemistryModel.reactions.reactionType import BasicReactionType


class MicroReversibleReaction(BasicReactionType):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        specieThermos,
        reactionRate,
        processIndeces,
        *args,
        **kwargs
    ):
        super(MicroReversibleReaction, self).update(
            specieThermos,
            reactionRate,
            processIndeces,
            *args,
            **kwargs
        )

    # Update method -----------------------------------------------------------
    def update(self, T):
        super(MicroReversibleReaction, self).update(T)
        self.Keq = self.Keq_()
        self.dKeqdT = self.dKeqdT_()

    # Properties --------------------------------------------------------------
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

    # Dissociation equilibrium constant ---------------------------------------
    def Keq_(self):
        '''Diss.-Recomb. equilibrium constants matrix'''
        PFdot = {
            name: np.sum(t.intPF.Q * t.transPF.Q) \
                for name, t in self.thSp.items()
            }
        K = PFdot['O']**2 / PFdot['O2']
        return K

    def dKeqdT_(self):
        '''Diss.-Recomb. equilibrium constants matrix'''
        PFdot = {
            name: np.sum(t.intPF.Q * t.transPF.Q) \
                for name, t in self.thSp.items()
            }
        dPFdotdT = {
            name: np.sum(t.intPF.dQdT * t.transPF.Q \
                + t.intPF.Q * t.transPF.dQdT) \
                for name, t in self.thSp.items()
            }
        K = 2*PFdot['O']*dPFdotdT['O']
        K -= self.Keq(T)*dPFdotdT['O2']
        return K/dPFdotdT['O2']

    # Reverse reaction rates --------------------------------------------------
    def kr(self, kf, reacIndex, indeces=None):
        if reacIndex == self.processIndeces['diss']:
            Keq_ = self.Keq
        else:
            l, r = indeces
            Keq_ = self.thSp['O2'].Q[r] / self.thSp['O2'].Q[l]
        return kf / Keq_

    def dkrdT(self, kr, dkfdT, reacIndex, indeces=None):
        if reacIndex == self.processIndeces['diss']:
            dkrdT_ = dkfdT / self.Keq
            dkrdT_ -= dkfdT / self.Keq**2 * self.dKeqdT
            return dkrdT_
        else:
            Q_, dQdT_ = self.thSp['O2'].Q, self.thSp['O2'].dQdT
            dkrdT_ = dkfdT * Q_[l] / Q_[r]
            dkrdT_ += kr * ( dQdT_[l] / Q_[l] - dQdT_[r] / Q_[r] )
            return dkrdT_
