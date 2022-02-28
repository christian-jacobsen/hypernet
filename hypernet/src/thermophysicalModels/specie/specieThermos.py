from hypernet.src.thermophysicalModels.specie import specie as specieMdl
from hypernet.src.thermophysicalModels.specie import partitionFun as PFMdl
from hypernet.src.thermophysicalModels.specie import thermo as thermoMdl
from hypernet.src.thermophysicalModels.specie import equationOfState as EOSMdl
from hypernet.src.general import utils


class SpecieThermos(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        name,
        rovib=None,
        thermo='CoupledEnergyModes',
        constVP='V',
        EOS='PerfectGas',
        *args,
        **kwargs
    ):
        self.specie = specieMdl.Specie(name, **rovib) \
            if rovib != None else specieMdl.Specie(name)

        self.intPF = PFMdl.internal.Internal(self.specie)
        self.transPF = PFMdl.translational.Translational(self.specie)

        self.thermo = utils.get_class(thermoMdl, thermo)(
            self.specie, self.intPF, constVP=constVP
        )

        self.EOS = utils.get_class(EOSMdl, EOS)(self.specie)

    # Methods
    ###########################################################################
    # Update ------------------------------------------------------------------
    def update(self, T):
        self.intPF.update(T)
        self.transPF.update(T)
        self.thermo.update(T)
        self.EOS.update(T)
