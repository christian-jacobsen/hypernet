import numpy as np

from hypernet.src.general import const
from hypernet.src.thermophysicalModels.specie.partitionFun import Basic


class Translational(Basic):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        specie,
        *args,
        **kwargs
    ):
        super(Translational, self).__init__(specie)
        # Constant basis
        self.basis = 2.*np.pi*const.UKB*self.specie.m/const.UH**2

    # Methods
    ###########################################################################
    def Q_(self, T):
        return np.ones(1)*np.power(self.basis*T, 3./2.)

    def dQ_dT_(self, T):
        return np.ones(1)*3./2.*self.basis*np.power(self.basis*T, 1./2.)
