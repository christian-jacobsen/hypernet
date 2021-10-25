import abc

from hypernet.src.general import utils
from hypernet.src.thermophysicalModels.chemistry import chemistryModel as chemMdl


class Basic(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        mixture,
        specieThermos,
        chemistryModel,
        reactionsList=None,
        processFlags=None,
        constPV='V',
        *args,
        **kwargs
    ):
        # Mixture
        self.mixture = mixture

        # Thermodynamic specie properties
        self.spTh = specieThermos

        # Constant pressure/volume process
        if constPV == 'P':
            self.cvp = self.mixture.cp
            self.cvp_i = self.mixture.cp_i
            self.dcvp_idT = self.mixture.dcp_idT
            self.eh = self.mixture.h
            self.dehdY = self.mixture.dhdY
            self.eh_i = self.mixture.h_i
        else:
            self.cvp = self.mixture.cv
            self.cvp_i = self.mixture.cv_i
            self.dcvp_idT = self.mixture.dcv_idT
            self.eh = self.mixture.e
            self.dehdY = self.mixture.dedY
            self.eh_i = self.mixture.e_i

        # Chemistry model
        self.chemModel = utils.get_class(chemMdl, chemistryModel)(
            self.spTh,
            processFlags,
            reactionsList=reactionsList,
            *args,
            **kwargs
        )

        # Variables
        self.n_var = self.chemModel.nSpecies + 1
        self.varIndices = [self.n_var-1]
        self.varNames = self.get_names()

    # Methods
    ###########################################################################
    # Update method -----------------------------------------------------------
    @abc.abstractmethod
    def update(self, *args, **kwargs):
        pass

    # Variables ---------------------------------------------------------------
    def get_names(self):
        varNames = []
        for name, spTh_ in self.spTh.items():
            if spTh_.specie.n_at > 1:
                varNames.extend(
                    [ name+'('+str(b+1)+')' for b in range(spTh_.specie.n_bins) ]
                )
            else:
                varNames.append(name)
        varNames.append('T')
        return varNames
