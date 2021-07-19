import numpy as np
import prode as pro

from hypernet.src.general import const
from hypernet.src.general import utils
from hypernet.src.thermophysicalModels.chemistryModel.surrogate.surrogate import Surrogate


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
        verbose=1
    ):
        super(DeepNet, self).__init__(mixture)
        self.inp = inputs

        # Data ================================================================
        self.data = utils.get_class(pro.data, self.inp.data['type'])(
            dim=self.inp.net['args']['dim'],
            train_data=train,
            valid_data=valid,
            test_data=test,
            training=self.inp.train_flg,
            **self.inp.data['args']
        )
        print(len(self.data.train))
        print(self.data.train[0].shape)
        print(self.data.train[0])

        # Network =============================================================
        self.net = utils.get_class(pro.networks, self.inp.net['architecture'])(
            **self.inp.net['args']
        )

        # Apply I/O Transformation --------------------------------------------
        self.net.apply_input_transform(self.inp.input_fn)
        self.net.apply_output_transform(self.inp.output_fn)

        # Model ===============================================================
        name = self.inp.name + '_' + self.inp.net['architecture']
        self.model = pro.model.Model(
            self.data, self.net, name, path=self.inp.path
        )

        # Build Model ---------------------------------------------------------
        self.model.build(
            self.net.inp_dim, verbose=verbose, load_net=self.inp.load_flg, \
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
            self.inp.callbacks['plotter']['var_names'] = self.var_names

        self.model.fit(epochs=self.inp.epochs, callbacks=self.inp.callbacks)

        # Save Model ----------------------------------------------------------
        if self.inp.save_flg > 0:
            self.model.save()

        # Training history ----------------------------------------------------
        self.model.visualize_training(
            self.model.train_dir+'/history.csv', limit_y=self.inp.limit_y
        )

        # Testing phase -------------------------------------------------------
        if self.inp.test_flg > 0:
            plot_style = dict(
                plotting=self.inp.plot_flg, var_names=self.var_names, \
                labels=self.inp.labels
            )
            for i, test_i in enumerate(self.data.test):
                # Define 'path'
                path = self.model.post_dir+'/testing/test_'+str(i+1)
                self.model.predict(
                    test_i, path=path, scales=self.inp.scales, **plot_style
                )
                if self.inp.log_plot:
                    self.model.predict(
                        test_i, path=path, scales=['log', 'log'], **plot_style
                    )

    def update(self, ic, x):
        x = [ic, x]
        if not isinstance(self.net.inp_dim, (list,tuple)):
            x = tf.concat(x, axis=1)
        Y = np.squeeze(self.net.predict(x).numpy())
        new_mix = {}
        for specie in self.mixture:
            new_mix[specie] = Y[self.sp_to_Y[specie]]
        self.mixture.update(new_mix)
