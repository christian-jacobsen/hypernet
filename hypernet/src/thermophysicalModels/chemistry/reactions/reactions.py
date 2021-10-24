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
        processIndeces,
        *args,
        **kwargs
    ):
        # Processes indeces
        self.processIndeces = processIndeces

        # Species thermos
        self.spTh = specieThermos

        # Reactions database
        self.reacDB = self.get_reacDB(pd.read_csv(reactList, index_col=0))

        print(self.reacDB)

        # Reaction rates
        reactRate = 'Arrhenius' #if 'Arrhenius' in self.reacDB['description']
        self.reacRate = utils.get_class(reacRateMdl, reactRate)(self.reacDB)

        # Reaction types
        reacType = 'MicroReversible' #if 'MicroReversible' in self.reacDB['description']
        self.reacType = utils.get_class(reacTypeMdl, reacType)(
            self.spTh, self.reacRate, self.processIndeces
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

    # Manipulate reactions database -------------------------------------------
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
            return indeces[0]
        else:
            return indeces
