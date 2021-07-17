import numpy as np

from hypernet.src.general import const
from hypernet.src.general import utils

from scipy import interpolate


class DeepNet(Surrogate):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        mixture,
        inputs,
        train=None,
        valid=None,
        test=None,
        verbose=0
    ):  

        self.mixture = mixture
        self.sp_to_Y = self.specie_to_Y(mixture)

        Y = []
        pos = 0
        for specie, thermo in mixture.items():
            Y_i = thermo.specie.Y
            for i in range(len(Y_i)):
                self.Y_to_sp[specie] = [ i+pos for i in range(len(Y_i)) ]
            pos = Y_to_sp[specie][-1]
            Y.append(thermo.specie.Y)
        Y = np.concatenate(tuple(Y))

        # Variables Name ----------------------------------------------------------
        var_names = [ r'$O_2^{({%s})}$' % (i+1) for i in range(dim['out']-1) ]
        var_names.append(r'$O$')


        self.inp = inputs

        # Data ================================================================
        self.data = pro_utils.get_class(pro.data, self.inp.data['type'])(
            dim=self.inp.net['args']['dim'],
            train_data=train,
            valid_data=valid,
            test_data=test,
            training=self.inp.train_flg,
            **self.inp.data['args']
        )

        # Network =============================================================
        self.net = pro_utils.get_class(pro.networks, self.inp.net['architecture'])(
            training=self.inp.train_flg,
            **self.inp.net['args']
        )

        # Apply I/O Transformation --------------------------------------------
        self.net.apply_input_transform(self.inp.input_fn)
        self.net.apply_output_transform(self.inp.output_fn)

        # Model ===============================================================
        self.model = pro.model.Model(
            self.data, self.net, self.inp.name, path='./chemSurrogateModel/'
        )

        # Build Model ---------------------------------------------------------
        self.model.build(
            self.inp.dim['inp'], verbose=verbose, load_net=self.inp.train_flg, \
            visualize_graph=False
        )

        # Compile Model -------------------------------------------------------
        self.model.compile(
            self.inp.optimizer, lr=self.inp.lr, loss=self.inp.loss, \
            metrics=self.inp.metrics, decay=self.inp.decay
        )

    # Methods
    ###########################################################################
    def fit(self):
        # Fit Model -----------------------------------------------------------
        if 'plotter' in self.inp.callbacks:
            self.inp.callbacks['plotter']['test_cases'] = self.data.test

        self.model.fit(epochs=self.inp.epochs, callbacks=self.inp.callbacks)

        # Save Model ----------------------------------------------------------
        if self.inp.save_flg > 0:
            self.model.save()

        # Training history ----------------------------------------------------
        self.model.visualize_training(
            self.model.train_dir+'/history.csv', limit_y=self.inp.limit_y
        )

        # Testing phase -------------------------------------------------------
        if inp.test_flg > 0:
            plot_style = dict(plotting=self.inp.plot_flg, var_names=self.inp.var_names, labels=self.inp.labels)
            for i, test_i in enumerate(data_test):
                model.predict(test_i, test_case_num=i, scales=inp.scales, **plot_style)
                if inp.log_plot:
                    model.predict(test_i, test_case_num=i, scales=['log', 'log'], **plot_style)
                

    def update(self, ic, x):
        x = [ic, x]
        if not input_is_list:
            x = tf.concat(x, axis=1)
        Y = np.squeeze(self.net.predict(x).numpy())
        new_mix = {}
        for specie in self.mixture:
            new_mix[specie] = Y[self.sp_to_Y[specie]]
        self.mixture.update(new_mix)

    def specie_to_Y(self, mixture):
        _sp_to_Y = {}
        pos = 0
        for specie, thermo in mixture.items():
            Y_sp = thermo.specie.Y
            for i in range(len(Y_sp)):
                _sp_to_Y[specie] = [ i+pos for i in range(len(Y_sp)) ]
            pos = _sp_to_Y[specie][-1] + 1 
        return _sp_to_Y



