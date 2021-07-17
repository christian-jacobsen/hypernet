class Surrogate(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        data_x,
        data_Y,
        *args,
        **kwargs
    ):
        # Data ================================================================
        self.x = data_x
        self.Y = data_Y
        self.n_species = data_Y.shape[1]

    # Methods
    ###########################################################################
    @abc.abstractmethod
    def fit(self):
        pass

    @abc.abstractmethod
    def update(self, x):
        pass