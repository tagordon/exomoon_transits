import numpy as np
import ctypes
from ctypes import byref
clib = ctypes.CDLL("../src/wrapper.so")

def grad_kepler_solve(M, e):
    
    clib.grad_kepler_solve.restype = None
    j = len(M)
    
    M = byref((ctypes.c_double * j).from_buffer(M))
    e = byref((ctypes.c_double(e)))
    cosf = (ctypes.c_double * j).from_buffer(np.zeros(j))
    sinf = (ctypes.c_double * j).from_buffer(np.zeros(j))
    f_e = (ctypes.c_double * j).from_buffer(np.zeros(j))
    f_M = (ctypes.c_double * j).from_buffer(np.zeros(j))
    
    j = byref(ctypes.c_int(j))
    
    clib.grad_kepler_solve(M, e, cosf, sinf, f_e, f_M, j)
    
    return np.array(cosf), np.array(sinf), np.array(f_e), np.array(f_M)

def kepler_solve(M, e):
    
    clib.kepler_solve.restype = None
    j = len(M)
    
    M = byref((ctypes.c_double * j).from_buffer(M))
    e = byref(ctypes.c_double(e))
    cosf = (ctypes.c_double * j).from_buffer(np.zeros(j))
    sinf = (ctypes.c_double * j).from_buffer(np.zeros(j))
    
    j = byref(ctypes.c_int(j))
    
    clib.kepler_solve(M, e, cosf, sinf, j)
    
    return np.array(cosf), np.array(sinf)

def kepler_solve_RPP(M, e):
    
    clib.kepler_solve_RPP.restype = None
    j = len(M)
    
    M = byref((ctypes.c_double * j).from_buffer(M))
    e = byref(ctypes.c_double(e))
    cosf = (ctypes.c_double * j).from_buffer(np.zeros(j))
    sinf = (ctypes.c_double * j).from_buffer(np.zeros(j))
    
    j = byref(ctypes.c_int(j))
    
    clib.kepler_solve_RPP(M, e, cosf, sinf, j)
    
    return np.array(cosf), np.array(sinf)

def impacts(t, ap, t0p, ep, Pp, wp, ip, am, t0m, em, Pm, Om, wm, im, mm):
        
    clib.impacts.restype = None
    j = len(t)
    
    bp = (ctypes.c_double * j).from_buffer(np.zeros(j))
    bpm = (ctypes.c_double * j).from_buffer(np.zeros(j))
    theta = (ctypes.c_double * j).from_buffer(np.zeros(j))
    
    t = byref((ctypes.c_double * j).from_buffer(t))
    j = byref(ctypes.c_int(j))
    
    ap = byref(ctypes.c_double(ap))
    t0p = byref(ctypes.c_double(t0p))
    ep = byref(ctypes.c_double(ep))
    Pp = byref(ctypes.c_double(Pp))
    wp = byref(ctypes.c_double(wp))
    ip = byref(ctypes.c_double(ip))
    am = byref(ctypes.c_double(am))
    t0m = byref(ctypes.c_double(t0m))
    em = byref(ctypes.c_double(em))
    Pm = byref(ctypes.c_double(Pm))
    Om = byref(ctypes.c_double(Om))
    wm = byref(ctypes.c_double(wm))
    im = byref(ctypes.c_double(im))
    mm = byref(ctypes.c_double(mm))
    
    clib.impacts(t, ap, t0p, ep, Pp, wp, ip, am, t0m, em, Pm, Om, wm, im, mm, j, 
                       bp, bpm, theta)
    
    return np.array(bp), np.array(bpm), np.array(theta)
    
def grad_impacts(t, ap, t0p, ep, Pp, wp, ip, am, t0m, em, Pm, Om, wm, im, mm):
    
    clib.grad_impacts.restype = None
    j = len(t)
    
    bp = (ctypes.c_double * j).from_buffer(np.zeros(j))
    bpm = (ctypes.c_double * j).from_buffer(np.zeros(j))
    theta = (ctypes.c_double * j).from_buffer(np.zeros(j))
    
    dbp = ((ctypes.c_double * j) * 14).from_buffer(np.zeros((j, 14)))
    dbpm = ((ctypes.c_double * j) * 14).from_buffer(np.zeros((j, 14)))
    dtheta = ((ctypes.c_double * j) * 14).from_buffer(np.zeros((j, 14)))
    
    t = byref((ctypes.c_double * j).from_buffer(t))
    j = byref(ctypes.c_int(j))
    
    ap = byref(ctypes.c_double(ap))
    t0p = byref(ctypes.c_double(t0p))
    ep = byref(ctypes.c_double(ep))
    Pp = byref(ctypes.c_double(Pp))
    wp = byref(ctypes.c_double(wp))
    ip = byref(ctypes.c_double(ip))
    am = byref(ctypes.c_double(am))
    t0m = byref(ctypes.c_double(t0m))
    em = byref(ctypes.c_double(em))
    Pm = byref(ctypes.c_double(Pm))
    Om = byref(ctypes.c_double(Om))
    wm = byref(ctypes.c_double(wm))
    im = byref(ctypes.c_double(im))
    mm = byref(ctypes.c_double(mm))
    
    clib.grad_impacts(t, ap, t0p, ep, Pp, wp, ip, am, t0m, em, Pm, Om, wm, im, mm, j, 
                       bp, bpm, theta, dbp, dbpm, dtheta)
    
    return np.array(bp), np.array(bpm), np.array(theta), np.array(dbp), np.array(dbpm), np.array(dtheta)

def coords(t, ap, t0p, ep, Pp, wp, ip, am, t0m, em, Pm, Om, wm, im, mm):
    
    clib.coords.restype = None
    j = len(t)
    
    xp = (ctypes.c_double * len(t)).from_buffer(np.zeros(j))
    yp = (ctypes.c_double * len(t)).from_buffer(np.zeros(j))
    zp = (ctypes.c_double * len(t)).from_buffer(np.zeros(j))
    xm = (ctypes.c_double * len(t)).from_buffer(np.zeros(j))
    ym = (ctypes.c_double * len(t)).from_buffer(np.zeros(j))
    zm = (ctypes.c_double * len(t)).from_buffer(np.zeros(j))
    t = byref((ctypes.c_double * len(t)).from_buffer(t))
    j = byref(ctypes.c_int(j))
    
    ap = byref(ctypes.c_double(ap))
    t0p = byref(ctypes.c_double(t0p))
    ep = byref(ctypes.c_double(ep))
    Pp = byref(ctypes.c_double(Pp))
    wp = byref(ctypes.c_double(wp))
    ip = byref(ctypes.c_double(ip))
    am = byref(ctypes.c_double(am))
    t0m = byref(ctypes.c_double(t0m))
    em = byref(ctypes.c_double(em))
    Pm = byref(ctypes.c_double(Pm))
    Om = byref(ctypes.c_double(Om))
    wm = byref(ctypes.c_double(wm))
    im = byref(ctypes.c_double(im))
    mm = byref(ctypes.c_double(mm))
    
    clib.coords(t, ap, t0p, ep, Pp, wp, ip, am, t0m, em, Pm, Om, wm, im, mm, j, xp, yp, zp, xm, ym, zm)
    return np.array(xp), np.array(yp), np.array(zp), np.array(xm), np.array(ym), np.array(zm)