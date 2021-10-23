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

        self.chem = chemModel(
            specieThermos,
            reactionsList,
            processFlags,
            *args,
            **kwargs
        )

    # Update method -----------------------------------------------------------
    def update(self, r, T, mass):
        # Update chemistry model
        self.chem.update(T)
        # Update mixture
        self.mixture.update(
            {name: np.take(r, idx)/mass for idx in self.chem.specieIndeces}
        )

    # Derivatives -------------------------------------------------------------
    def derivatives(self, r, T, mass):

        # Evaluate contributions from reactions
        drdt = self.wr(r)

        # Evaluate the effect on the thermodynamic system
        dTdt = self.wT(drdt, T, mass)

        return np.concatenate(tuple([drdt, dTdt]))

    def wr(self, r):

        # Get Master Equation matrices
        Ke, Kd, Kr = self.chem.K

        # Evaluate contributions from each process
        wr_ = np.matmul(Ke, r) * r[-1]
        wr_ += np.matmul(Kd, r) * r[-1]
        wr_ += np.squeeze(Kr) * r[-1]**3

        return wr_

    def wT(self, drdt, T, mass):

        # Get mass fractions derivative
        dYdt = {name: np.take(drdt, idx)/mass for idx in self.specieIndeces}

        # Evaluate source term
        wT_ = - self.dehdY(T, dYdt) / self.cvp(T)

        return wT_

    # Jacobian ----------------------------------------------------------------
    def jacobian(self, r, T, mass):
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


    def dwrdr(self, r):
        # Get Master Equation matrices
        Ke, Kd, Kr = self.chem.K

        # Get Master Equation matrices
        I = np.array([[0]*self.chem.spTh['O2'].n_bins+[1]])
        wr_ = np.expand_dims(self.wr(r), 1)

        dwrdr_ = r[-1] * (Ke + Kd)
        dwrdr_ += np.matmul(wr_ + 2*r[-1]**3*Kr, I) / r[-1]

        return dwrdr_

    def dwrdT(self, r):

        # Get Master Equation matrices
        dKedT, dKddT, dKrdT = self.chem.dKdT

        # Evaluate contributions from each process
        dwrdT_ = np.matmul(dKedT, r) * r[-1]
        dwrdT_ += np.matmul(dKddT, r) * r[-1]
        dwrdT_ += np.squeeze(dKrdT) * r[-1]**3

        return dwrdT_



    def domegaYdY(self, Y, T):
        # Get Master Equation matrices
        K_e, K_d, K_r = self.chem.K

    def jac(self, t, y, arg):
        '''Jacobian calculation: jac[i,j] = df[i] / dy[j].'''
        d = np.shape(y)[0]
        J = np.zeros((d,d), dtype=np.float64)
        J[:,:-1] = (arg[0][:,:-1] + arg[1][:,:-1]) * y[-1]
        J[:, -1] = np.matmul(arg[0][:,:-1], y[:-1]) \
            + np.matmul(arg[1][:,:-1], y[:-1]) \
            + 3*np.squeeze(arg[2])*np.power(y[-1], 2)
        return J



        # Evaluate contributions from each process
        excit = np.matmul(K_e, Y) * Y[-1]
        diss = np.matmul(K_d, Y) * Y[-1]
        recom = np.squeeze(K_r) * np.power(Y[-1], 3)

        return excit + diss + recom

    def domegaYdT(self, Y, dYdt, T):
        dYdt = {name: np.take(dYdt, idx) for idx in self.specieIndeces}
        return - self.dehdY(T, dYdt) / self.cvp(T)

    def domegaTdY(self, Y, T):
        # Get Master Equation matrices
        K_e, K_d, K_r = self.matrices(T)

        # Evaluate contributions from each process
        excit = np.matmul(K_e, Y) * Y[-1]
        diss = np.matmul(K_d, Y) * Y[-1]
        recom = np.squeeze(K_r) * np.power(Y[-1], 3)

        return excit + diss + recom

    def domegaTdT(self, Y, dYdt, T):
        dYdt = {name: np.take(dYdt, idx) for idx in self.specieIndeces}
        return - self.dehdY(T, dYdt) / self.cvp(T)




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