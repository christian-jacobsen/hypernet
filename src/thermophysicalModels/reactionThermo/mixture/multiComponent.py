import numpy as np
import hypernet.src.general.const as const
from hypernet.src.thermophysicalModels.specie import Specie
from hypernet.src.thermophysicalModels.specie import Specie

class Basic(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        input_file,
        mixture
    ):
        self.mix = {}
        for s, Y in mixture.items():
            sp = Specie(input_file, s)
            sp.Y_(Y)
            self.mix[s] = Thermo(sp)
        
    # Methods
    ###########################################################################
    # Enthalpy ================================================================
    def cp(self,T):
        # [J/(kg K)]
        return self.mix_quantity('cp',T)

    def h(self,T):
        # [J/kg]
        return self.mix_quantity('h',T)

    # Internal Energy =========================================================
    def cv(self,T):
        # [J/(kg K)]
        return self.mix_quantity('cv',T)

    def e(self,T):
        # [J/kg]
        return self.mix_quantity('e',T)

    # Energy of formation -----------------------------------------------------
    def e_f(self):
        # [J/kg]
        return self.mix_quantity('e_f',T)

    # Translational Internal Energy -------------------------------------------
    def cv_tr(self):
        # [J/(kg K)]
        return self.mix_quantity('cv_tr',T)

    def e_tr(self,T):
        # [J/kg]
        return self.mix_quantity('e_tr',T)

    # Ro-Vibrational Internal Energy ------------------------------------------
    def cv_rv(self,T):
        # [J/(kg K)]
        return self.mix_quantity('cv_rv',T)

    def e_rv(self,T):
        # [J/kg]
        return self.mix_quantity('e_rv',T)


    def m(self):
        # [kg/mol]
        return 1./self.mix_quantity(1.)

    def R(self):
        # [J/(kg K)]
        return self.mix_quantity(const.URG)

    def X_i(self):
        for specie, thermo in self.mix.items():
            specie.X_(np.clip(thermo.sp.Y * self.m() / thermo.sp.m), a_min=0.)

    def n_i(self, rho):
        for specie, thermo in self.mix.items():
            specie.n_(np.clip(thermo.sp.Y * const.UNA * rho / thermo.sp.m), a_min=0.)

    def rho_i(self, rho):
        for specie, thermo in self.mix.items():
            specie.rho_(np.clip(thermo.sp.Y * rho), a_min=0.)

    def p_i(self, rho):
        for specie, thermo in self.mix.items():
            specie.p_(np.clip(thermo.sp.Y * rho), a_min=0.)


    return max(Yi*NA.value()*rho/speciesData_[speciei].W(), 0);


max(Yi*molWeightMixture(celli)/speciesData_[speciei].W(), 0);

        
        //- Update values of molar-fractions from mass-fractions for cell-set
        virtual scalar molarFraction(const label speciei, const scalar Yi, const label celli);
        
        //- Update values of molar-fractions from mass-fractions for patch
        virtual scalar molarFraction(const label speciei, const scalar Yi, const label patchi, const label facei);
        
        //- Update values of mass-fractions from molar-fractions during the initialisation
        virtual scalar massFractionFromMolarFraction(const label speciei, /*const scalar celli,*/ const scalar Xi);
        
        //- Update values of mass-fractions from partial densities
        virtual scalar massFractionFromPartialDensity(const scalar rhoi, const scalar p, const scalar Tt);
        
        //- Update values of number densities from mass fractions
        virtual scalar numberDensity(const label speciei, const scalar Yi, const scalar rho);
        
        //- Update values of partial pressures from molar fractions
        virtual scalar partialPressure(const scalar Xi, const scalar p);
        
        //- Update values of partial pressures from the equation of state
        virtual scalar partialPressureEoS(const label speciei, const scalar rhoi, const scalar Ti);
        
        //- Update values of partial densities from mass fractions
        virtual scalar partialDensity(const scalar Yi, const scalar rho);
        



    def convert_to_array(x):        
        if not isinstance(x, np.ndarray):
            x = np.array(x)
        return x

    def mix_quantity(self, quantity, args):
        var = []
        for specie, thermo in self.mix.items():
            if isinstance(quantity, str):
                var_ = getattr(thermo, quantity)(*args)
            else:
                var_ = quantity
            var_ = utils.convert_to_array(var_)
            var.append(thermo.sp.Y / thermo.sp.m * var_)
        return np.sum(np.concatenate(tuple(var)))