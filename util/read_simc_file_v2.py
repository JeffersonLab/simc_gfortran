import numpy as np
from scipy.io import FortranFile
ff=FortranFile('eep_hydrogen_q8.bin','r')
nvar = int(ff.read_ints('int32'))
#names=np.empty(nvar,dtype='str')

for i in range(1,nvar+1):
    name = str(ff.read_record('a16'))
    print(i,name)
            
for i in range(1,nvar+1):
    values=float(ff.read_reals('f8'))
    print(i,values)


