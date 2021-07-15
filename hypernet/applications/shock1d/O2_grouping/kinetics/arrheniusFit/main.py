import os
import sys
import numpy as np
import pandas as pd
import inputs as inp
# Remove the following line if you install `hypernet` as a Python package -----
sys.path.append(inp.hypernet)
import hypernet as hy

from scipy.optimize import curve_fit
from dataGenerator import DataGenerator
from hypernet.src.general import const
from hypernet.src.general import utils
from hypernet.src.thermophysicalModels.specie import specie as specie_module

import fortranformat as ff
fmt = ff.FortranRecordWriter('(E16.8)')


# Arrhenius functions
###############################################################################
def modified_arrhenius(T, a, b, c):
    return a * np.power(T,b) * np.exp(-c/(const.URG*T))

def modified_arrhenius_lin(T, a_log, b, c):
    # a_log = np.log(a)
    return a_log + b*np.log(T) - c/(const.URG*T)

def modified_arrhenius_inv(T_inv, a, b, c):
    # T_inv = 1/T
    return a * np.power(1/T_inv,b) * np.exp(-c/const.URG*T_inv)

# Write function
###############################################################################
def write(array, name, idx, path):
    file = open(path, 'w')
    file.write('Units=cm^3/s\n')
    for i, row in enumerate(array):
        s = name[i]+':'
        row = fmt.write(row).split("\n")
        for j in row:
            s += j.strip()+','
        s += str(idx[i])+'\n'
        file.write(s)
    file.close()

# Main function
###############################################################################
@utils.main_decorator
def main():

    # Initialize species ======================================================
    utils.print_main("Initializing species")
    species = {}
    for sp_name in inp.mixture:
        species[sp_name] = specie_module.Specie(sp_name, **inp.specie)

    # Generate Data ===========================================================
    utils.print_main('Generating data')
    dataGen = DataGenerator(inp.T, species, **inp.chemistry)
    T, K = dataGen.training_data()

    # Linear Fit ==============================================================
    utils.print_main('Fitting rates')
    for j in range(K.shape[1]):
        K_log_j = np.log( K[:,j][ K[:,j] != 0. ] )
        T_j     = T[:,0][ K[:,j] != 0. ]

        param_j, _ = curve_fit(modified_arrhenius_lin, T_j, K_log_j, \
            p0=[1,1,1.e4], method='trf')
        param_j[0] = np.exp(param_j[0])
        if j == 0:
            param = np.array(param_j)
        else:
            param = np.vstack((param, np.array(param_j)))

    # Write Arrhenius coefficient =============================================
    if not os.path.exists(inp.path):
        os.makedirs(inp.path)
    path = inp.path + inp.system+'_'+inp.grouping_name
    write(param, dataGen.rates_names, dataGen.reac_idx, path)


if __name__ == "__main__":
    main()
