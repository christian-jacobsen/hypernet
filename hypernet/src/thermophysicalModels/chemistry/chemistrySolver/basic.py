import abc

from hypernet.src.general import utils
from hypernet.src.thermophysicalModels.chemistry import chemistryModel as chemMdl


class Basic(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        mixture,
        chemistryModel,
        specieThermos,
        reactionsList=None,
        processFlags=None,
        constPV='V',
        *args,
        **kwargs
    ):
        # Mixture
        self.mixture = mixture

        # Constant pressure/volume process
        if constPV == 'P':
            self.cvp = self.mixture.cp
            self.dehdY = self.mixture.dhdY
        else:
            self.cvp = self.mixture.cv
            self.dehdY = self.mixture.dedY

        # Chemistry model
        self.chemModel = utils.get_class(chemMdl, chemistryModel)(
            specieThermos,
            reactionsList,
            processFlags,
            *args,
            **kwargs
        )

    # Methods
    ###########################################################################
    # Update method -----------------------------------------------------------
    @abc.abstractmethod
    def update(self, *args, **kwargs):
        pass