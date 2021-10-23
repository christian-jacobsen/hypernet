import abc


class Basic(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        mixture
    ):
        # Mixture =============================================================
        self.species = list(mixture.keys())
        self.sp_to_Y = self.specie_to_Y(mixture)
        self.var_names = self.get_var_names(mixture)

    # Methods
    ###########################################################################
    @abc.abstractmethod
    def fit(self):
        pass

    @abc.abstractmethod
    def update(self, x):
        pass

    def specie_to_Y(self, mixture):
        _sp_to_Y = {}
        pos = 0
        for specie, thermo in mixture.items():
            _sp_to_Y[specie] = [ i+pos for i in range(len(thermo.specie.Y)) ]
            pos = _sp_to_Y[specie][-1] + 1 
        return _sp_to_Y

    def get_var_names(self, mixture):
        names = []
        for specie, thermo in mixture.items():
            if hasattr(thermo.specie, 'n_bins') and thermo.specie.n_bins > 1:
                names.extend(
                    [ r'$%s^{({%s})}$' % (specie, i+1) \
                        for i in range(thermo.specie.n_bins) ]
                )
            else:
                names.append(r'$%s$' % specie)
        return names

