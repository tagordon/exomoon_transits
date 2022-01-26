import numpy as np
import ctypes
from ctypes import byref
clib = ctypes.CDLL("../fortran/wrapper.so")

__args__ = ['a_p', 't0_p', 'ecc_p', 'per_p', 'lop_p', 'inc_p', 
            'a_m', 't0_m', 'ecc_m', 'per_m', 'lan_m', 'lop_m', 'inc_m', 'mass_m']

def validate_elements(argdict):
    
    if set(argdict.keys()) == set(__args__):
        
        argdict = {k: byref(ctypes.c_double(v)) for k, v in argdict.items()}
        idm = {v: i for i, v in enumerate(__args__)}
        args = tuple(dict(sorted(argdict.items(), 
                                 key=lambda p: idm[p[0]])).values())
        
        return args
        
    else:
        
        raise ValueError("required parameters are: ", __args__)

def impacts(t, argdict):
    
    j = len(t)
    bp, bpm, theta = tuple((ctypes.c_double * j).from_buffer(np.zeros(j)) for a in range(3))
    t = byref((ctypes.c_double * j).from_buffer(t))
    
    args = validate_elements(argdict)
    clib.impacts.restype = None
    clib.impacts(t, *args, byref(ctypes.c_int(j)), bp, bpm, theta)
    return map(np.array, (bp, bpm, theta))

    
def grad_impacts(t, argdict):
    
    j = len(t)
    bp, bpm, theta = tuple((ctypes.c_double * j).from_buffer(np.zeros(j)) for a in range(3))
    dbp, dbpm, dtheta = tuple(((ctypes.c_double * j) * 14).from_buffer(np.zeros((j, 14))) for a in range(3))
    t = byref((ctypes.c_double * j).from_buffer(t))
    
    args = validate_elements(argdict)
    clib.grad_impacts.restype = None
    clib.grad_impacts(t, *args, byref(ctypes.c_int(j)), 
                        bp, bpm, theta, dbp, dbpm, dtheta)
        
    return map(np.array, (bp, bpm, theta, dbp, dbpm, dtheta))

def coords(t, argdict):
    
    clib.coords.restype = None
    j = len(t)
    
    xp, yp, zp, xm, ym, zm = tuple((ctypes.c_double * j).from_buffer(np.zeros(j)) for a in range(6))
    t = byref((ctypes.c_double * len(t)).from_buffer(t))
    
    args = validate_elements(argdict)
    
    clib.coords(t, *args, byref(ctypes.c_int(j)), 
                xp, yp, zp, xm, ym, zm)
    
    return map(np.array, (xp, yp, zp, xm, ym, zm))
