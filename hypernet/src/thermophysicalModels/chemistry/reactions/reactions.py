import re
import numpy as np
import pandas as pd
import tensorflow as tf

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
        reactRateList = ('Arrhenius',)
        for rr in reactRateList:
            if rr in self.reacDB['description'][0]:
                self.reacRate = utils.get_class(reacRateMdl, rr)(self.reacDB)
                reactRate = rr
                break

        # Reaction types
        reacTypeList = ('MicroReversible',)
        for rt in reacTypeList:
            if rt in self.reacDB['description'][0]:
                self.reacType = utils.get_class(reacTypeMdl, rt)(
                    self.spTh, self.reacRate, self.processIndices
                )
                break
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
        K = tf.nest.map_structure(lambda x: np.array(x), K)
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
