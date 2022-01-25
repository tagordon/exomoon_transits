import numpy as np
import ctypes
from ctypes import byref
clib = ctypes.CDLL("../src/wrapper.so")

def flux(c1, c2, rp, rm, bp, bpm, cth, sth):
    
    j = len(bp)
    
    bp = byref((ctypes.c_double * j).from_buffer(bp))
    bpm = byref((ctypes.c_double * j).from_buffer(bpm))
    cth = byref((ctypes.c_double * j).from_buffer(cth))
    sth = byref((ctypes.c_double * j).from_buffer(sth))
    lc = ((ctypes.c_double * 8) * j).from_buffer(np.zeros((8, j)))
    rp = byref(ctypes.c_double(rp))
    rm = byref(ctypes.c_double(rm))
    c1 = byref(ctypes.c_double(c1))
    c2 = byref(ctypes.c_double(c2))
    j = byref(ctypes.c_int(j))
    clib.flux.restype = None
    clib.flux(c1, c2, rp, rm, bp, bpm, cth, sth, lc, j)
    return np.array(lc)