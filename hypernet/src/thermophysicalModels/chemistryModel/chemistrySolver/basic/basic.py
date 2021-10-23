import numpy as np
import pandas as pd

from hypernet.src.general import const
from hypernet.src.general import utils

from hypernet.src.thermophysicalModels.specie import specie as specieModule
from hypernet.src.thermophysicalModels.specie import partitionFun as PFModule
from hypernet.src.thermophysicalModels.specie import thermo as thermoModule
from hypernet.src.thermophysicalModels.specie import equationOfState as EOSModule



class ChemistryModel(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        mixture,
        specieThermos,
        reactionsList,
        processFlags,
        constPV='V',
        *args,
        **kwargs
    ):
        # Mixture
        self.mixture = mixture

        # Constant pressure/volume process
        if constPV == 'P':
            self.cvp = self.mixture.cp
            self.dehdY = self.mixture.dhdY
        else:
            self.cvp = self.mixture.cv
            self.dehdY = self.mixture.dedY

        # Thermodynamic specie properties
        self.spTh = specieThermos

        # Reactive processes
        self.processFlags = processFlags
        self.processIndeces = {
            'diss': [2]
            'excit': [5, 6]
        }

        # Reactions
        self.reactions = Reactions(
            specieThermos,
            reactionsList,
            self.processIndeces,
            *args,
            **kwargs
        )

        # Species
        self.nSpecies, self.specieIndeces = self.species()

    # Species details ---------------------------------------------------------
    def species(self):
        '''Get all the rates matrices for the ODE.'''
        nSpecies, specieIndeces = 0, {}
        for name, spTh_ in self.spTh.values():
            if spTh_.specie.n_at > 1:
                n = spTh_.specie.n_bins
            else:
                n = 1
            specieIndeces[name] = list(range(nSpecies, nSpecies+n))
            nSpecies += n
        return nSpecies, specieIndeces

    # Master Equation matrices ------------------------------------------------
    def matrices(self, T):
        '''Get all the rates matrices for the ODE.'''
        reac = self.reactions.update(T)
        m = self.spTh['O'].specie.m

        K_e = self.excit(self.processFlags['excit'], reac) / m
        K_d = self.diss(self.processFlags['diss'], reac) / m
        K_r = self.recom(self.processFlags['diss'], reac) / m**2 * 2
        return K_e, K_d, K_r

    def excit(self, mask, reac):
        '''Excit. & Relax. rates matrix'''

        reac = reac.loc[reac['reacIndex'].isin(self.reacIndeces['excit'])]

        # Obtain relaxation rates
        rates = rates_coeff['exch'] + rates_coeff['inel']
        idx   = list(map(tuple, np.transpose(np.where(rates == 0.0))))
        for i,j in idx:
            if i != j:
                rates[i,j] = rates[j,i] * Q_bins[j] / Q_bins[i]

        # Construct Excit. & Relax. matrix
        d = np.shape(rates)[0] + 1
        K = np.zeros((d,d), dtype=np.float64)
        if mask:
            K[:-1,:-1] = -np.diag(np.sum(rates, axis=1)) + np.transpose(rates)
        return K

    def diss(self, mask, reactions):
        '''Dissociation rates matrix'''

        # Obtain Dissociation rates
        reac = reac.loc[reac['reacIndex'].isin(self.reacIndeces['diss'])]

        # Construct Dissociation matrix
        K = np.zeros((self.n_species,self.n_species), dtype=np.float64)
        if mask:
            K[:-1,:-1] = np.diag(-rates)
            K[ -1,:-1] = rates
        return K

    def recom(self, mask, reactions):
        '''Recombination rates matrix'''

        # Obtain Recombination rates
        rates = np.multiply(rates_coeff['diss'], self.Keq_diss(Q_bins, T))

        # Construct Recombination matrix
        d = np.shape(rates)[0] + 1
        K = np.zeros((d,1), dtype=np.float64)
        if mask:
            K[:-1,0] = rates
            K[ -1,0] = -np.sum(rates)
        return K

    # Derivatives -------------------------------------------------------------
    def derivatives(self, t, y, args):
        Y, T = y[:-1], y[-1]

        # Update mixture
        self.mixture.update(
            {name: np.take(Y, idx) for idx in self.specieIndeces}
        )

        # Evaluate contributions from reactions
        dYdt = self.omegaY(Y, T)

        # Evaluate the effect on the thermodynamic system
        dTdt = self.omegaT(Y, dYdt, T)

        return np.concatenate(tuple([dYdt, dTdt]))

    def omegaY(self, Y, T):
        # Get Master Equation matrices
        K_e, K_d, K_r = self.matrices(T)

        # Evaluate contributions from each process
        excit = np.matmul(K_e, Y) * Y[-1]
        diss = np.matmul(K_d, Y) * Y[-1]
        recom = np.squeeze(K_r) * np.power(Y[-1], 3)

        return excit + diss + recom

    def omegaT(self, Y, dYdt, T):
        dYdt = {name: np.take(dYdt, idx) for idx in self.specieIndeces}
        return - self.dehdY(T, dYdt) / self.cvp(T)

    # Jacobian ----------------------------------------------------------------
    def jacobian(self, t, y, args):
        Y, T = y[:-1], y[-1]

        # Evaluate contributions from reactions
        dYdt = self.omega(Y, T)

        # Evaluate the effect on the thermodynamic system
        # >> Update mixture
        self.mixture.update(
            {name: np.take(Y, idx) for idx in self.specieIndeces}
        )
        # >> Evaluate mixture Cv
        cv_ = self.cvp(T)
        # >> Evaluate dTdt
        dYdt_ = {name: np.take(omega_, idx) for idx in self.specieIndeces}
        dedY_ = self.mixture.dedY(T, dYdt_)
        dTdt = - dedY_ / self.cvp(T)
        return excit + diss + recom


    def derivatives(self, t, y, args):
        for 
        self.mixture.update(self, mixture, var='Y'):





template<class ThermoType>
void Foam::standardChemistryModel<ThermoType>::jacobian
(
    const scalar t,
    const scalarField& c,
    const label li,
    scalarField& dcdt,
    scalarSquareMatrix& J
) const
{
    const scalar T = c[nSpecie_];
    const scalar p = c[nSpecie_ + 1];

    forAll(c_, i)
    {
        c_[i] = max(c[i], 0);
    }

    dcdt = Zero;
    J = Zero;

    // Evaluate contributions from reactions
    forAll(reactions_, ri)
    {
        const Reaction<ThermoType>& R = reactions_[ri];
        scalar omegaI, kfwd, kbwd;
        const labelList null;
        R.dwdc(p, T, c_, li, J, dcdt, omegaI, kfwd, kbwd, false, null);
        R.dwdT(p, T, c_, li, omegaI, kfwd, kbwd, J, false, null, nSpecie_);
    }

    // Evaluate the effect on the thermodynamic system ...

    // c*Cp
    scalar ccp = 0, dccpdT = 0;
    for (label i = 0; i < nSpecie_; i++)
    {
        ccp += c_[i]*specieThermos_[i].cp(p, T);
        dccpdT += c_[i]*specieThermos_[i].dcpdT(p, T);
    }

    // dT/dt
    scalar& dTdt = dcdt[nSpecie_];
    for (label i = 0; i < nSpecie_; i++)
    {
        dTdt -= dcdt[i]*specieThermos_[i].ha(p, T);
    }
    dTdt /= ccp;

    // dp/dt = 0 (pressure is assumed constant)

    // d(dTdt)/dc
    for (label i = 0; i < nSpecie_; i++)
    {
        scalar& d2Tdtdci = J(nSpecie_, i);
        for (label j = 0; j < nSpecie_; j++)
        {
            const scalar d2cjdtdci = J(j, i);
            d2Tdtdci -= d2cjdtdci*specieThermos_[j].ha(p, T);
        }
        d2Tdtdci -= specieThermos_[i].cp(p, T)*dTdt;
        d2Tdtdci /= ccp;
    }

    // d(dTdt)/dT
    scalar& d2TdtdT = J(nSpecie_, nSpecie_);
    for (label i = 0; i < nSpecie_; i++)
    {
        const scalar d2cidtdT = J(i, nSpecie_);
        d2TdtdT -=
            dcdt[i]*specieThermos_[i].cp(p, T)
          + d2cidtdT*specieThermos_[i].ha(p, T);
    }
    d2TdtdT -= dTdt*dccpdT;
    d2TdtdT /= ccp;

    // d(dpdt)/dc = 0 (pressure is assumed constant)

    // d(dpdt)/dT = 0 (pressure is assumed constant)
}