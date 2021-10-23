import numpy as np

from hypernet.src.general import const
from hypernet.src.general import utils
from hypernet.src.thermophysicalModels.specie import specie as specieMdl
from hypernet.src.thermophysicalModels.specie import partitionFun as PFMdl
from hypernet.src.thermophysicalModels.specie import thermo as thermoMdl
from hypernet.src.thermophysicalModels.specie import equationOfState as EOSMdl


class SpecieThermos(object):

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
        self.specie = specieMdl.Specie(name, rovib)

        self.intPF = PFMdl.internalPF.Internal(self.specie)
        self.transPF = PFMdl.translationalPF.Translational(self.specie)

        self.thermo = utils.get_class(thermoMdl, thermo)(
            self.specie, self.intPF
        )

        self.EOS = utils.get_class(EOSMdl, EOS)(self.specie)
