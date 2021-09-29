import abc


class PartitionFun(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        specie,
        *args,
        **kwargs
    ):
        # Specie Properties ===================================================
        self.specie = specie

    # Partition functions -----------------------------------------------------
    @abc.abstractmethod
    def Q(self, T):
        pass

    @abc.abstractmethod
    def dQ_dT(self, T):
        pass
