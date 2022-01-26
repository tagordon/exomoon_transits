import numpy as np
from astropy import constants as ac

from kep import impacts, coords, grad_impacts
from phot import flux

__all__ = ['star', 'planet', 'mon', 'system']

class star:
    
    ms = ac.M_sun.value / ac.M_earth.value
    
    def __init__(self, c1, c2):
        
        self.c1 = c1
        self.c2 = c2
    
class planet:
    
    def __init__(self, r, a, t0, ecc, per, lop, inc):
        
        self.r = r
        self.a = a
        self.t = t0
        self.e = ecc
        self.p = per
        self.w = lop
        self.i = inc
        
        self.args = (a, t0, ecc, per, lop, inc)
        argnames = ['a', 't0', 'ecc', 'per', 'lop', 'inc']
        self.argnames = tuple([n + '_p' for n in argnames])
        self.argdict = dict(zip(self.argnames, self.args))
    
class moon:
    
    def __init__(self, r, a, t0, ecc, per, lan, lop, inc, mass):
        
        self.r = r
        self.a = a
        self.t = t0
        self.e = ecc
        self.p = per
        self.O = lan
        self.w = lop
        self.i = inc
        self.m = mass
        
        self.args = (a, t0, ecc, per, lan, lop, inc, mass)
        argnames = ['a', 't0', 'ecc', 'per', 'lan', 'lop', 'inc', 'mass']
        self.argnames = tuple([n + '_m' for n in argnames])
        self.argdict = dict(zip(self.argnames, self.args))

class system:
    
    def __init__(self, star, planet, moon):
        
        self.star = star
        self.planet = planet
        self.moon = moon
        self.arglist = ()
        self._lc = None
        self._grad = {}
        self._primary_coords = None
        self._secondary_coords = None
        self._computed = False
        self._computed_coords = False
        
    def compute(self, t, grad=False, coords=False):
        
        p = self.planet
        m = self.moon
        
        if coords:
            xp, yp, zp = coords(t, *(p.args + m.args))
            self._primary_coords = (xp, yp, zp)
            self._secondary_coords = (xm, ym, zm)
            self._computed_coords = True
        else:
            if grad:
                bp, bpm, theta, dbp, dbpm, dtheta = grad_impacts(t, {**p.argdict, **m.argdict})
                lc = flux(self.star.c1, self.star.c2, p.r, m.r, bp, bpm, np.cos(theta), np.sin(theta)).T
                self._lc = lc[ 0]
                f_bp = lc[3]
                f_bpm = lc[4]
                f_theta = lc[5]
                
                df = f_bp * dbp + f_bpm * dbpm + f_theta * dtheta
                self._grad = {(p.argnames + m.argnames)[i]: df[i] for i in range(np.shape(df)[0])}
                self._grad['rm'] = lc[1]
                self._grad['rp'] = lc[2]
                self._grad['c1'] = lc[6]
                self._grad['c2'] = lc[7]
            else:
                #bp, bpm, theta = impacts(t, *(p.args + m.args))
                bp, bpm, theta = impacts(t, {**p.argdict, **m.argdict})
                self._lc = flux(self.star.c1, self.star.c2, p.r, m.r, bp, bpm, np.cos(theta), np.sin(theta))[:, 0]
        
        self._computed = True
            
    def lightcurve(self):
        
        if self._computed:
            return self._lc
        else:
            raise Exception("System lightcurve has not been computed. ", 
                            "Call system.compute(t) to compute the lightcurve.")
            
    def grad(self):
        if self._computed:
            return self._grad
        else:
            raise Exception("System lightcurve has not been computed. ", 
                            "Call system.compute(t) to compute the lightcurve.")
            
    def primary_coords(self):
        
        if self._computed_coords:
            return self._primary_coords
        else:
            raise Exception("System coordinates have not been computed. ",
                            "Call system.compute(t, coords=True) to compute the coordinates")
    
    def secondary_coords(self):
        
        if self._computed_coords:
            return self._secondary_coords
        else:
            raise Exception("System coordinates have not been computed. ",
                            "Call system.compute(t, coords=True) to compute the coordinates")
            
    def loglike(self, y, sigma):
        
        mu = self.lightcurve()
        s2 = sigma * sigma
        return -0.5 * np.sum((y - mu) ** 2 / s2 + np.log(s2))