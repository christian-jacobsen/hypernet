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
        # Constant base
        self.base = 2.*np.pi*self.specie.m*const.UKB/(const.UH**2)

    # Methods
    ###########################################################################
    def Q_(self, T):
        return np.power(self.base*T, 3./2.)

    def dQdT_(self, T):
        return 3./2.*self.base*np.power(self.base*T, 1./2.)
