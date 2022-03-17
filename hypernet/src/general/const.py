import numpy as np

# Vacuum constants ------------------------------------------------------------
UEPS0       = 8.8541878128e-12                      # Permittivity [F/m]
UMU0        = 1.25663706212e-6                      # Magnetic permeability
UC0         = 1.e0/np.sqrt(UEPS0*UMU0)              # Speed of light

# Global constants ------------------------------------------------------------
UKB         = 1.38064852e-23                     	# Boltzmann's constant [J/K]
UNA         = 6.02214076e23                         # Avogadro's number [1/mol]
URG         = UKB*UNA                               # Universal gas constant [J/(mol K)]
UH          = 6.62607015e-34                        # Planck's constant [J s]
UHT         = UH/(2.e0*np.pi)                       # Reduced Plancks's constant [J s]
UAMU        = 1.66053906660e-27                     # Atomic mass unit [kg]

# Electron constants ----------------------------------------------------------
UE          = 1.602176634e-19                       # Charge [C]
GEM         = 2.e0                                  # Spin
UME         = 9.1093837015e-31                      # Mass [kg]
UMME        = UME*UNA                               # Molar mass [kg/mol]
R_EM        = URG/UMME                              # Specific gas constant [J/(kg K)]
CV_EM       = 1.5e0*R_EM                            # Constant volume specific heat [J/(kg K)]
CP_EM       = (CV_EM + R_EM)                        # Constant pressure specific heat  [J/(kg K)]
UALPHA      = UE**2/(2.e0*UH*UC0*UEPS0)             # Fine structure constant

# Conversion factors ----------------------------------------------------------
CM_to_K     = 1.e2*UH*UC0/UKB                       # [cm^-1] -> [K]
AS_to_MS    = 1.e-20                                # [A^2]   -> [m^2]
ATM_to_PA   = 101325.e0                             # [atm]   -> [Pa]
EH_to_EV    = 27.211386245988e0                     # [Eh]    -> [eV] Hartree to Electronvolt
EV_to_J     = 1.602176634e-19                       # [eV]    -> [J]  Electronvolt to Joule
EH_to_J     = EH_to_EV*EV_to_J                      # [Eh]    -> [J]  Hartree to Joule

# Hydrogen constants ----------------------------------------------------------
UIH         = 13.598434599702e0                     # Hydrogen ionization potential [eV]
UA0         = 4.e0*np.pi*UEPS0*UHT**2/(UME*UE**2)   # Bohr radius [m]
