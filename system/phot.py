import numpy as np
import ctypes
from ctypes import byref
clib = ctypes.CDLL("../fortran/wrapper.so")

def flux(c1, c2, rp, rm, bp, bpm, cth, sth):
    
    j = len(bp)
    
    bp, bpm, cth, sth = map(lambda a: byref((ctypes.c_double * j).from_buffer(a)), 
                            (bp, bpm, cth, sth, ))
    
    rp, rm, c1, c2 = map(lambda a: byref(ctypes.c_double(a)), (rp, rm, c1, c2))
    
    lc = ((ctypes.c_double * 8) * j).from_buffer(np.zeros((8, j)))
    j = byref(ctypes.c_int(j))
    clib.flux.restype = None
    clib.flux(c1, c2, rp, rm, bp, bpm, cth, sth, lc, j)
    return np.array(lc)