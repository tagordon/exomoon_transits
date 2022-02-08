import numpy as np
from copy import copy
import sys
sys.path.append('../gefera')
from kep import impacts, grad_impacts
import ctypes
from ctypes import byref
clib = ctypes.CDLL("../fortran/wrapper.so")

clib.kepler_solve.restype = None
clib.kepler_solve_RPP.restype = None

uniform_draw = lambda l, h: np.random.rand() * (h - l) + l

def random_args():
    
    keys = [
        'ab', 
        'tb', 
        'eb', 
        'pb', 
        'wb', 
        'ib',      
        'am', 
        'tm', 
        'em', 
        'pm', 
        'om', 
        'wm', 
        'im', 
        'mm'
    ]
    
    ab = uniform_draw(1, 10)
    tb = uniform_draw(-10, 10)
    eb = uniform_draw(0, 0.999)
    pb = uniform_draw(10, 100)
    wb = uniform_draw(0, np.pi*2)
    ib = uniform_draw(0, np.pi/2)
    
    am = uniform_draw(0.01, 1)
    tm = uniform_draw(-10, 10)
    em = uniform_draw(0, 0.999)
    pm = uniform_draw(1, 10)
    om = uniform_draw(0, np.pi*2)
    wm = uniform_draw(0, np.pi*2)
    im = uniform_draw(0, np.pi/2)
    mm = uniform_draw(0, 1)
    
    argdict = {
        k: v for k, v in zip(
            keys, 
            (ab, tb, eb, pb, wb, ib, am, tm, em, pm, om, wm, im, mm)
        )
    }
    
    return argdict

def gradient(p):
    d = 0.000001
    t = np.linspace(0, 100, 10)
    argdict = random_args()
    ind = np.where(np.array(list(argdict.keys())) == p)[0][0]
    pargdict = copy(argdict)
    margdict = copy(argdict)
    pargdict[p] = argdict[p] + d
    margdict[p] = argdict[p] - d
    pbp, pbpm, ptheta = impacts(t, pargdict)
    mbp, mbpm, mtheta = impacts(t, margdict)
    _, _, _, dbp, dbpm, dtheta = grad_impacts(t, argdict)
    fdiff_bp = (pbp - mbp) / (2 * d)
    fdiff_bpm = (pbpm - mbpm) / (2 * d)
    fdiff_theta = (ptheta - mtheta) / (2 * d)
    
    tol = 1e-3
    
    if p == 'pb':
        tol = 0.1
    if p == 'pm':
        tol = 0.5
    
    assert np.all(np.isclose(fdiff_bp, dbp[ind], atol=tol))
    assert np.all(np.isclose(fdiff_bpm, dbpm[ind], atol=tol))
    assert np.all(np.isclose(fdiff_theta, dtheta[ind], atol=tol))

def RPP_vs_newton(M, ecc):
    
    j = len(M)
    M = byref((ctypes.c_double * j).from_buffer(M))
    
    for e in ecc:
        e = byref(ctypes.c_double(e))
        cosf_newton = (ctypes.c_double * j).from_buffer(np.zeros(j))
        sinf_newton = (ctypes.c_double * j).from_buffer(np.zeros(j))
        cosf_RPP = (ctypes.c_double * j).from_buffer(np.zeros(j))
        sinf_RPP = (ctypes.c_double * j).from_buffer(np.zeros(j))
        newton = clib.kepler_solve(M, e, cosf_newton, sinf_newton, byref(ctypes.c_int(j)))
        RPP = clib.kepler_solve_RPP(M, e, cosf_RPP, sinf_RPP, byref(ctypes.c_int(j)))
        assert np.all(np.isclose(np.array(sinf_RPP), np.array(sinf_newton), atol=1e-12))
        assert np.all(np.isclose(np.array(cosf_RPP), np.array(cosf_newton), atol=1e-12))

