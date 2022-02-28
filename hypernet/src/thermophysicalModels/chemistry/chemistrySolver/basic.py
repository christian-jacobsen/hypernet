import abc
import numpy as np

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
        heatBath='isothermal',
        *args,
        **kwargs
    ):
        # Mixture
        self.mixture = mixture

        # Thermodynamic specie properties
        self.spTh = specieThermos

        # Isothermal/Adiabatic heat bath
        self.heatBath = heatBath

        # Chemistry model
        self.chemModel = utils.get_class(chemMdl, chemistryModel)(
            self.spTh,
            processFlags,
            reactionsList=reactionsList,
            *args,
            **kwargs
        )

        # Variables
        self.varNames = self.get_names()
        self.extraVars = dict(p=[], n=[], E=[])

    # Methods
    ###########################################################################
    # Update method -----------------------------------------------------------
    @abc.abstractmethod
    def update(self, *args, **kwargs):
        pass

    # Variables ---------------------------------------------------------------
    def get_names(self):
        varNames = []
        for name, spTh in self.spTh.items():
            if spTh.specie.n_at > 1:
                varNames.extend([
                    'Y_'+name+'('+str(b+1)+')' \
                        for b in range(spTh.specie.n_bins)
                ])
            else:
                varNames.append('Y_'+name)
        varNames.append('T')
        return tuple(varNames)

    # Extra physical quantities -----------------------------------------------
    def eval_extra_vars(self, T, rho):
        self.extraVars['p'].append(float(self.mixture.p_(rho, T)))
        self.extraVars['n'].append(float(self.mixture.n_(rho)))
        self.extraVars['E'].append(float(self.mixture.he_()))