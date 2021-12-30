import numpy as np
import ctypes
clib = ctypes.CDLL("../src/wrapper.so")

def flux(c1, c2, rp, rm, bp, bpm, cth, sth):
    
    bp = (ctypes.c_double * len(bp)).from_buffer(bp)
    bpm = (ctypes.c_double * len(bpm)).from_buffer(bpm)
    cth = (ctypes.c_double * len(cth)).from_buffer(cth)
    sth = (ctypes.c_double * len(sth)).from_buffer(sth)
    lc = ((ctypes.c_double * 8) * len(bp)).from_buffer(np.zeros((8, len(bp))))
    rp = ctypes.c_double(rp)
    rm = ctypes.c_double(rm)
    c1 = ctypes.c_double(c1)
    c2 = ctypes.c_double(c2)
    j = ctypes.c_int(len(bp))
    clib.LC.restype = None
    clib.LC(c1, c2, rp, rm, bp, bpm, cth, sth, lc, j)
    return np.array(lc)