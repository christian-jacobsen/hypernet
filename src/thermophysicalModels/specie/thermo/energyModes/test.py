from specie import specie


path = '/home/zanardi/WORKSPACE/CFD/databases/air/thermo/'
mix = ['O2', 'O']

s = specie(path, mix)
print(s.database[mix[1]])