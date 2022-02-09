import numpy as np
import ctypes
from ctypes import byref
clib = ctypes.CDLL("wrapper.so")

import sys
sys.path.append('../gefera')
from kep import impacts, coords
from phot import flux

__all__ = ['BarycenterOrbit', 'MoonOrbit', 'System']
    
class BarycenterOrbit:
    
    def __init__(self, ab, tb, eb, pb, wb, ib):
        
        self.pnames = [
            'ab', 
            'tb', 
            'eb', 
            'pb', 
            'wb', 
            'ib'
        ]
        self.pdict = dict(
            zip(
                self.pnames,
                (ab, tb, eb, pb, wb, ib)
            )
        )
    
class MoonOrbit:
    
    def __init__(self, am, tm, em, pm, om, wm, im, mm):
        
        self.pnames = [
            'am', 
            'tm', 
            'em', 
            'pm', 
            'om', 
            'wm', 
            'im', 
            'mm'
        ]
        self.pdict = dict(
            zip(
                self.pnames, 
                (am, tm, em, pm, om, wm, im, mm)
            )
        )

class System:
    
    def __init__(self, bo, mo):
        
        self.bo = bo
        self.mo = mo
        self._lc = None
        self._grad = {}
        self._primary_coords = None
        self._secondary_coords = None
        self._computed = False
        self._computed_coords = False
    
    def lightcurve(self, t, u1, u2, rp, rm):
        
        au_r = 215.03215567054764
        
        bo = self.bo
        mo = self.mo
        
        bp, bpm, theta = impacts(
            t, 
            {**bo.pdict, **mo.pdict}
        )
        _, _, _, xm, ym, _ = coords(
            t, 
            {**bo.pdict, **mo.pdict}
        )
        bm = np.sqrt(xm**2 + ym**2) * au_r
        
        lc_nomoon = flux(
            u1, 
            u2, 
            rp, 
            0.0, 
            bp, 
            bpm, 
            np.cos(theta), 
            np.sin(theta)
        )[:, 0]
        
        j = len(t)
        bp = (ctypes.c_double * j).from_buffer(bp)
        bm = (ctypes.c_double * j).from_buffer(bm)
        bpm = (ctypes.c_double * j).from_buffer(bpm)
        f = (ctypes.c_double * j).from_buffer(np.zeros(j))

        clib.flux.restype = None

        clib.flux(
            f, 
            bp, 
            bm, 
            bpm, 
            ctypes.c_double(rp), 
            ctypes.c_double(rm), 
            ctypes.c_double(0.0), 
            ctypes.c_double(u1 + 2 * u2), 
            ctypes.c_double(0.0), 
            ctypes.c_double(-u2),
            ctypes.c_int(j)
        )
        lc_noplanet = f
        
        self._lc = lc_noplanet + lc_nomoon  
        
        return self._lc