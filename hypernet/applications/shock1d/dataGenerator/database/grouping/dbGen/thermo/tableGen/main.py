import os
import sys
import shutil
import numpy as np
import inputs as inp
# Remove the following line if you install `hypernet` as a Python package -----
sys.path.append(inp.hypernet)
import hypernet as hy

from matplotlib import pyplot as plt
from hypernet.src.general import utils
from hypernet.src.thermophysicalModels.specie import specie as specie_module
from hypernet.src.thermophysicalModels.specie import thermo as thermo_module

import fortranformat as ff
fmt = ff.FortranRecordWriter('(E16.10)')


# Writing function
###############################################################################
def write(data, path):
    file = open(path, "w")
    for d in range(data.shape[1]):
        file.write(str(d+1)+'\n')
        file.write(fmt.write(data[:,d])+'\n')
    file.close()

# Main function
###############################################################################
@utils.main_decorator
def main():

    # Define temperature vector ===============================================
    utils.print_main("Defining T vector")
    T = np.arange(inp.T[0], inp.T[-1]+1., 25., dtype=np.float64)

    # Evaluate thermo =========================================================
    thermo, cvint, eint, qint = {}, {}, {}, {}
    for sp in inp.species:
        utils.print_main("Evaluating thermo for "+sp)
        thermo[sp] = utils.get_class(thermo_module, inp.thermo)(
                specie_module.Specie(sp, **inp.specie)
            )
        _cvint, _eint, _qint = [], [], []

        utils.print_submain("Copying thermo properties file")
        if thermo[sp].specie.n_at == 1:
            shutil.copyfile(inp.specie['thermo_path']+sp, inp.path+sp)
            file1 = open(inp.specie['thermo_path']+sp, 'r')
            file2 = open(inp.path+sp, 'w')
            l_brk = -1
            for l, line in enumerate(file1.readlines()):
                if l == l_brk:
                    break
                elif line.startswith('MODEL'):
                    file2.write('MODEL = NONE\n')
                elif line.startswith('NB_ELEC_LEVELS'):
                    file2.write('NB_ELEC_LEVELS = 1\n')
                    l_brk = l+2
                else:
                    file2.write(line)
            file1.close()
            file2.close()
        else:
            for i in range(1,thermo[sp].specie.n_bins+1):
                file1 = open(inp.specie['thermo_path']+sp, 'r')
                file2 = open(inp.path+sp+'_'+str(i), 'w')
                for line in file1.readlines():
                    if line.startswith('#Electronic'):
                        break
                    elif line.startswith('NAME'):
                        file2.write('NAME = '+sp+'_'+str(i)+'\n')
                    elif line.startswith('MODEL'):
                        file2.write('MODEL = GROUP\n')
                    elif line.startswith('THETA_ROT'):
                        file2.write('THETA_ROT = 0.0d0\n')
                    elif line.startswith('THETA_VIB'):
                        file2.write('THETA_VIB = 0.0d0\n')
                    else:
                        file2.write(line)
                file1.close()
                file2.close()

        if thermo[sp].specie.n_at > 1:
            utils.print_submain("Evaluating `e_rv`, `cv_rv` and `Q`")
            for T_i in T:
                _cvint.append(thermo[sp].cv_rv(T_i) / thermo[sp].specie.m)
                _eint.append(thermo[sp].e_rv(T_i) / thermo[sp].specie.m)
                _qint.append(thermo[sp].Q(T_i))
            cvint[sp] = np.vstack(tuple(_cvint))
            eint[sp] = np.vstack(tuple(_eint))
            qint[sp] = np.vstack(tuple(_qint))

            # Write tables ====================================================
            utils.print_submain("Writing tables")
            path = inp.path+str(thermo[sp].specie.n_bins)+'g/'
            if not os.path.exists(path):
                os.makedirs(path)

            file = open(path+'README', 'w')
            file.write('Do not erase this file')
            file.close()

            for tab in inp.tab_name['temp']:
                file = open(path+tab+'_table.dat', 'w')
                file.write(fmt.write(T))
                file.close()

            for tab in inp.tab_name['groups']:
                data = eval(tab)[sp]
                path_i = path+sp+'_'+tab+'.dat'
                write(data, path_i)


if __name__ == "__main__":
    main()
