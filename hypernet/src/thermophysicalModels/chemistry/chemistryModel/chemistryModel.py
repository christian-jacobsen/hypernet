import numpy as np
import pandas as pd

from hypernet.src.general import const
from hypernet.src.general import utils

from hypernet.src.thermophysicalModels.specie import specie as specieModule
from hypernet.src.thermophysicalModels.specie import partitionFun as PFModule
from hypernet.src.thermophysicalModels.specie import thermo as thermoModule
from hypernet.src.thermophysicalModels.specie import equationOfState as EOSModule


class ThermoSpecie(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        name,
        rovib=None,
        thermo='CoupledEnergyModes',
        EOS='PerfectGas',
        *args,
        **kwargs
    ):
        self.specie = specieModule.Specie(name, rovib)

        self.intPF = PFModule.internalPF.Internal(self.specie)
        self.transPF = PFModule.translationalPF.Translational(self.specie)

        self.thermo = utils.get_class(thermoModule, thermo)(
            self.specie, self.intPF
        )

        self.EOS = utils.get_class(EOSModule, EOS)(self.specie)
