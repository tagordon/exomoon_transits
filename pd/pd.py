import numpy as np
from timeit import default_timer as timer
import ctypes
import sys
sys.path.append('../gefera')
import gefera
clib = ctypes.CDLL("./photodynam.so")
clib.flux.restype = None

def flux(xp, yp, xm, ym, rp, rm, u1, u2):

    j = len(xp)
    xp = (ctypes.c_double * j).from_buffer(xp)
    yp = (ctypes.c_double * j).from_buffer(yp)
    xm = (ctypes.c_double * j).from_buffer(xm)
    ym = (ctypes.c_double * j).from_buffer(ym)
    f = (ctypes.c_double * j).from_buffer(np.zeros(j))
    j = ctypes.c_int(j)

    rp = ctypes.c_double(rp)
    rm = ctypes.c_double(rm)
    u1 = ctypes.c_double(u1)
    u2 = ctypes.c_double(u2)

    clib.flux(f, xp, yp, xm, ym, rp, rm, u1, u2, j)
    return np.array(f)

def time(t, u1, u2, rp, rm, ntimes, gefera_model):
    
    p, m = gefera_model.coords(t)
    xp, yp, _ = p
    xm, ym, _ = m
    start = timer()
    for _ in range(ntimes):
        lc = self.phot(t, u1, u2, rp, rm, bp, bpm, theta, grad=True)
    end = timer()
    return (end - start) / ntimes
    
    