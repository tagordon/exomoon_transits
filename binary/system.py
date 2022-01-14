import numpy as np
from astropy import constants as ac

from kep import impacts, coords, grad_impacts
from phot import flux

__all__ = ['star', 'rock', 'system']

class star:
    
    ms = ac.M_sun.value / ac.M_earth.value
    
    def __init__(self, mass, c1, c2):
        
        self.m = mass * self.ms
        self.c1 = c1
        self.c2 = c2
    
class rock:
    
    def __init__(self, ident, r, mass, t0, ecc, per, lan, lop, inc):
        
        self.ident = ident
        self.r = r
        self.m = mass
        self.t = t0
        self.e = ecc
        self.p = per
        self.O = lan
        self.w = lop
        self.i = inc
        
        self.args = (t0, ecc, per, lop, lan, inc, mass)
        argnames = ['t0', 'ecc', 'per', 'lop', 'lan', 'inc', 'mass']
        self.argnames = tuple([n + str(ident) for n in argnames])

class system:
    
    def __init__(self, star, primary, secondary):
        
        self.star = star
        self.primary = primary
        self.secondary = secondary
        self.arglist = ()
        self._lc = None
        self._grad = {}
        self._primary_coords = None
        self._secondary_coords = None
        self._computed = False
        self._computed_coords = False
        
    def compute(self, t, grad=False, coords=False):
        
        p = self.primary
        s = self.secondary
        
        if coords:
            xp, yp, zp = coords(t, self.star.m, *(p.args + s.args))
            self._primary_coords = (xp, yp, zp)
            self._secondary_coords = (xm, ym, zm)
            self._computed_coords = True
        else:
            if grad:
                bp, bpm, theta, dbp, dbpm, dtheta = grad_impacts(t, self.star.m, *(p.args + s.args))
                lc = flux(self.star.c1, self.star.c2, p.r, s.r, bp, bpm, np.cos(theta), np.sin(theta))
                self._lc = lc[:, 0]
                f_bp = lc[:, 3]
                f_bpm = lc[:, 4]
                f_theta = lc[:, 5]
                
                df = f_bp[:, None] * dbp + f_bpm[:, None] * dbpm + f_theta[:, None] * dtheta
                self._grad = {(('ms',) + p.argnames + s.argnames)[i]: df[:, i] for i in range(np.shape(df)[1])}
                self._grad['rm'] = lc[:, 1]
                self._grad['rp'] = lc[:, 2]
                self._grad['c1'] = lc[:, 6]
                self._grad['c2'] = lc[:, 7]
            else:
                bp, bpm, theta = impacts(t, self.star.m, *(p.args + s.args))
                self._lc = flux(self.star.c1, self.star.c2, p.r, s.r, bp, bpm, np.cos(theta), np.sin(theta))[:, 0]
        
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