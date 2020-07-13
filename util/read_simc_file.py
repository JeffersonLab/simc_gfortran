import numpy as np
#from scipy.io import FortranFile
#ff=FortranFile('eep_hydrogen_q8.bin','r','>u4')
#vals = ff.read_ints('int32')
with open('eep_hydrogen_q8.bin','rb') as f:
    header = np.fromfile(f, dtype=np.int32,count=1)
    nvar = np.fromfile(f, dtype=np.int32, count=1)
#    varname =np.empty(nvar)
#    values =np.empty(nvar)
    i=1
    while i <= nvar:
        varname=np.fromfile(f, dtype='a16', count=1)
        i=i+1
        values=np.fromfile(f, dtype='f32',count=1)
        print(i,varname,values)
print(nvar)


