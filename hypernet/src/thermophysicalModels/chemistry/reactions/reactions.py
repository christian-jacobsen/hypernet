import re
import numpy as np
import pandas as pd

from hypernet.src.general import const
from hypernet.src.general import utils
from hypernet.src.thermophysicalModels.chemistry.reactions import reactionRate as reacRateMdl
from hypernet.src.thermophysicalModels.chemistry.reactions import reactionType as reacTypeMdl


class Reactions(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        specieThermos,
        reactList,
        processIndices,
        *args,
        **kwargs
    ):
        # Processes indices
        self.processIndices = processIndices

        # Species thermos
        self.spTh = specieThermos

        # Reactions database
        self.reacDB = self.get_reacDB(pd.read_csv(reactList))

        # Reaction rates
        reactRate = 'Arrhenius' #if 'Arrhenius' in self.reacDB['description']
        self.reacRate = utils.get_class(reacRateMdl, reactRate)(self.reacDB)

        # Reaction types
        reacType = 'MicroReversible' #if 'MicroReversible' in self.reacDB['description']
        self.reacType = utils.get_class(reacTypeMdl, reacType)(
            self.spTh, self.reacRate, self.processIndices
        )

    # Methods
    ###########################################################################
    # Update method -----------------------------------------------------------
    def update(self, T):
        self.reacType.update(T)
        K = dict(
            kf=self.reacType.kf, dkfdT=self.reacType.dkfdT, kr=[], dkrdT=[]
        )
        for i, row in self.reacDB.iterrows():
            K['kr'].append(
                self.reacType.kr_(
                    K['kf'][i], row['reacIndex'], row['indices']
                )
            )
            K['dkrdT'].append(
                self.reacType.dkrdT_(
                    K['kr'][i], K['dkfdT'][i], row['reacIndex'], row['indices']
                )
            )
        K = {k: np.array(v) for k, v in K}
        return pd.concat([self.reacDB, pd.DataFrame.from_dict(K)], axis=1)

    # Manipulate reactions database -------------------------------------------
    def get_reacDB(self, reacDB):
        indices = []
        for i, row in reacDB.iterrows():
            indices.append(
                self.read_indices(row['reactStr'], row['reacIndex'])
            )
        reacDB['indices'] = indices
        return reacDB

    def read_indices(self, reacStr, reacIndex):
        indices = re.findall(r"\((.*?)\)",reacStr)
        indices = [int(i)-1 for i in indices]
        if reacIndex == self.processIndices['diss']:
            return indices[0]
        else:
            return indices
