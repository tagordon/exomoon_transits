import numpy as np
from astropy import constants as ac

from kep import input_coords, coords
from phot import flux

__all__ = ['star', 'rock', 'system']

class star:
    
    ms = ac.M_sun.value / ac.M_earth.value
    
    def __init__(self, mass, c1, c2):
        
        self.m = mass * self.ms
        self.c1 = c1
        self.c2 = c2
    
class rock:
    
    def __init__(self, r, mass, t0, ecc, per, lan, lop, inc):
        
        self.r = r
        self.m = mass
        self.t = t0
        self.e = ecc
        self.p = per
        self.O = lan
        self.w = lop
        self.i = inc

class system:
    
    def __init__(self, star, primary, secondary):
        
        self.star = star
        self.primary = primary
        self.secondary = secondary
        self._lc = None
        self._primary_coords = None
        self._secondary_coords = None
        self._computed = False
        self._computed_lc = False
        
    def compute(self, t, lc=True):
        
        p = self.primary
        s = self.secondary
        
        xp, yp, zp, xm, ym, zm, bp, bpm, theta = input_coords(t, self.star.m, 
                                      p.t, p.e, p.p, p.O, p.w, p.i, p.m, 
                                      s.t, s.e, s.p, s.O, s.w, s.i, s.m, 
                                      coords=True)
            
        self._primary_coords = (xp, yp, zp)
        self._secondary_coords = (xm, ym, zm)
        
        if lc:
            self._lc = flux(self.star.c1, self.star.c2, p.r, s.r, 
                            bp, bpm, np.cos(theta), np.sin(theta)) 
            self._computed_lc = True
        
        self._computed = True
            
    def lightcurve(self):
        
        if self._computed_lc:
            return self._lc[:, 0]
        else:
            raise Exception("System lightcurve has not been computed. ", 
                            "Call system.compute(t, lc=True) to compute the lightcurve.")
            
    def derivatives(self):
        if self._computed_lc:
            return self._lc[:, 1:]
        else:
            raise Exception("System lightcurve has not been computed. ", 
                            "Call system.compute(t, lc=True) to compute the lightcurve.")
            
    def primary_coords(self):
        
        if self._computed:
            return self._primary_coords
        else:
            raise Exception("System coordinates have not been computed. ",
                            "Call system.compute(t) to compute the coordinates")
    
    def secondary_coords(self):
        
        if self._computed:
            return self._secondary_coords
        else:
            raise Exception("System coordinates have not been computed. ",
                            "Call system.compute(t) to compute the coordinates")