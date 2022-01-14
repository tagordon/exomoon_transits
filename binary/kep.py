import numpy as np
import ctypes
clib = ctypes.CDLL("../src/wrapper.so")

def grad_kepler_solve(M, e):
    
    clib.grad_kepler_solve.restype = None
    j = len(M)
    
    M = (ctypes.c_double * j).from_buffer(M)
    e = ctypes.c_double(e)
    cosf = (ctypes.c_double * j).from_buffer(np.zeros(j))
    sinf = (ctypes.c_double * j).from_buffer(np.zeros(j))
    f_e = (ctypes.c_double * j).from_buffer(np.zeros(j))
    f_M = (ctypes.c_double * j).from_buffer(np.zeros(j))
    
    j = ctypes.c_int(j)
    
    clib.grad_kepler_solve_(M, e, cosf, sinf, f_e, f_M, j)
    
    return np.array(cosf), np.array(sinf), np.array(f_e), np.array(f_M)

def kepler_solve(M, e):
    
    clib.kepler_solve.restype = None
    j = len(M)
    
    M = (ctypes.c_double * j).from_buffer(M)
    e = ctypes.c_double(e)
    cosf = (ctypes.c_double * j).from_buffer(np.zeros(j))
    sinf = (ctypes.c_double * j).from_buffer(np.zeros(j))
    
    j = ctypes.c_int(j)
    
    clib.kepler_solve_(M, e, cosf, sinf, j)
    
    return np.array(cosf), np.array(sinf)

def impacts(t, ms, t0p, ep, Pp, Op, wp, ip, mp, t0m, em, Pm, Om, wm, im, mm):
        
    clib.impacts.restype = None
    j = len(t)
    
    bp = (ctypes.c_double * j).from_buffer(np.zeros(j))
    bpm = (ctypes.c_double * j).from_buffer(np.zeros(j))
    theta = (ctypes.c_double * j).from_buffer(np.zeros(j))
    
    t = (ctypes.c_double * j).from_buffer(t)
    j = ctypes.c_int(j)
    
    ms = ctypes.c_double(ms)
    t0p = ctypes.c_double(t0p)
    ep = ctypes.c_double(ep)
    Pp = ctypes.c_double(Pp)
    Op = ctypes.c_double(Op)
    wp = ctypes.c_double(wp)
    ip = ctypes.c_double(ip)
    mp = ctypes.c_double(mp)
    t0m = ctypes.c_double(t0m)
    em = ctypes.c_double(em)
    Pm = ctypes.c_double(Pm)
    Om = ctypes.c_double(Om)
    wm = ctypes.c_double(wm)
    im = ctypes.c_double(im)
    mm = ctypes.c_double(mm)
    
    clib.impacts_(t, ms, t0p, ep, Pp, Op, wp, ip, mp, t0m, em, Pm, Om, wm, im, mm, j, 
                       bp, bpm, theta)
    
    return np.array(bp), np.array(bpm), np.array(theta)
    
def grad_impacts(t, ms, t0p, ep, Pp, Op, wp, ip, mp, t0m, em, Pm, Om, wm, im, mm):
    
    clib.grad_impacts_.restype = None
    j = len(t)
    
    bp = (ctypes.c_double * j).from_buffer(np.zeros(j))
    bpm = (ctypes.c_double * j).from_buffer(np.zeros(j))
    theta = (ctypes.c_double * j).from_buffer(np.zeros(j))
    
    dbp = ((ctypes.c_double * 15) * j).from_buffer(np.zeros((15, j)))
    dbpm = ((ctypes.c_double * 15) * j).from_buffer(np.zeros((15, j)))
    dtheta = ((ctypes.c_double * 15) * j).from_buffer(np.zeros((15, j)))
    
    t = (ctypes.c_double * j).from_buffer(t)
    j = ctypes.c_int(j)
    
    ms = ctypes.c_double(ms)
    t0p = ctypes.c_double(t0p)
    ep = ctypes.c_double(ep)
    Pp = ctypes.c_double(Pp)
    Op = ctypes.c_double(Op)
    wp = ctypes.c_double(wp)
    ip = ctypes.c_double(ip)
    mp = ctypes.c_double(mp)
    t0m = ctypes.c_double(t0m)
    em = ctypes.c_double(em)
    Pm = ctypes.c_double(Pm)
    Om = ctypes.c_double(Om)
    wm = ctypes.c_double(wm)
    im = ctypes.c_double(im)
    mm = ctypes.c_double(mm)
    
    clib.grad_impacts_(t, ms, t0p, ep, Pp, Op, wp, ip, mp, t0m, em, Pm, Om, wm, im, mm, j, 
                       bp, bpm, theta, dbp, dbpm, dtheta)
    
    return np.array(bp), np.array(bpm), np.array(theta), np.array(dbp), np.array(dbpm), np.array(dtheta)

def coords(t, ms, t0p, ep, Pp, Op, wp, ip, mp, t0m, em, Pm, Om, wm, im, mm):
    
    clib.coords_.restype = None
    j = len(t)
    
    xp = (ctypes.c_double * len(t)).from_buffer(np.zeros(j))
    yp = (ctypes.c_double * len(t)).from_buffer(np.zeros(j))
    zp = (ctypes.c_double * len(t)).from_buffer(np.zeros(j))
    xm = (ctypes.c_double * len(t)).from_buffer(np.zeros(j))
    ym = (ctypes.c_double * len(t)).from_buffer(np.zeros(j))
    zm = (ctypes.c_double * len(t)).from_buffer(np.zeros(j))
    t = (ctypes.c_double * len(t)).from_buffer(t)
    j = ctypes.c_int(j)
    
    ms = ctypes.c_double(ms)
    t0p = ctypes.c_double(t0p)
    ep = ctypes.c_double(ep)
    Pp = ctypes.c_double(Pp)
    Op = ctypes.c_double(Op)
    wp = ctypes.c_double(wp)
    ip = ctypes.c_double(ip)
    mp = ctypes.c_double(mp)
    t0m = ctypes.c_double(t0m)
    em = ctypes.c_double(em)
    Pm = ctypes.c_double(Pm)
    Om = ctypes.c_double(Om)
    wm = ctypes.c_double(wm)
    im = ctypes.c_double(im)
    mm = ctypes.c_double(mm)
    
    clib.coords_(t, ms, t0p, ep, Pp, Op, wp, ip, mp, t0m, em, Pm, Om, wm, im, mm, j, xp, yp, zp, xm, ym, zm)
    return np.array(xp), np.array(yp), np.array(zp), np.array(xm), np.array(ym), np.array(zm)