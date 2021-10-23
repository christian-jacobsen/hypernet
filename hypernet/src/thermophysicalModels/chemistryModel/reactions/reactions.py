import re
import numpy as np
import pandas as pd

from hypernet.src.general import const
from hypernet.src.general import utils
from hypernet.src.thermophysicalModels.chemistryModel.reactions import reactionRate as reacRateModule
from hypernet.src.thermophysicalModels.chemistryModel.reactions import reactionType as reactionTypeModule


class Reactions(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        specieThermos,
        reactList,
        processIndeces,
        *args,
        **kwargs
    ):
        self.processIndeces = processIndeces
        self.spTh = specieThermos
        self.reacDB = self.get_reacDB(
            pd.read_csv(reactList, index_col=0)
        )

        # Reaction rates ------------------------------------------------------
        reactRate = 'ArrheniusReactionRate' #if 'Arrhenius' in self.reacDB['description']
        self.reacRate = utils.get_class(reacRateModule, reactionRate)(
            self.reacDB
        )

        # Reaction types ------------------------------------------------------
        reactType = 'MicroReversibleReaction' #if 'microReversible' in self.reacDB['description']
        self.reacType = utils.get_class(reacTypeModule, reactionType)(
            self.spTh, self.reacRate, self.processIndeces
        )

    # Methods
    ###########################################################################
    def update(self, T):
        self.reacType.update(T)
        K = dict(
            kf=self.reacType.kf, dkfdT=self.reacType.dkfdT, kr=[], dkrdT=[]
        )
        for i, row in self.reacDB.iterrows():
            K['kr'].append(
                self.reacType.kr(
                    K['kf'][i], row['reacIndex'], row['indeces']
                )
            )
            K['dkrdT'].append(
                self.reacType.kr(
                    K['kr'][i], K['dkfdT'][i], row['reacIndex'], row['indeces']
                )
            )
        K = {k: np.array(v) for k, v in K}
        return pd.concat([self.reacDB, pd.DataFrame.from_dict(K)], axis=1)

    def get_reacDB(self, reacDB):
        indeces = []
        for i, row in reacDB.iterrows():
            indeces.append(
                self.read_indeces(row['reactStr'], row['reacIndex'])
            )
        reacDB['indeces'] = indeces
        return reacDB

    def read_indeces(self, reacStr, reacIndex):
        indeces = re.findall(r"\((.*?)\)",reacStr)
        indeces = [int(i)-1 for i in indeces]
        if reacIndex == self.processIndeces['diss']:
            return *indeces
        else:
            return indeces
