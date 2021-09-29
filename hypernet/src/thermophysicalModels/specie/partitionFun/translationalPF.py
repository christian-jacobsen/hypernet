import numpy as np

from hypernet.src.general import const
from hypernet.src.specie.partitionFun.basicPF import PartitionFun


class Translational(PartitionFun):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        specie,
        *args,
        **kwargs
    ):
        super(Translational, self).__init__(specie)
        self.basis = 2.*np.pi*const.UKB*self.specie.m/const.UH**2

    # Translational partition functions ---------------------------------------
    def Q(self, T):
        return np.power(self.basis*T, 3./2.)

    def dQ_dT(self, T):
        return 3./2.*self.basis*np.power(self.basis*T, 1./2.)