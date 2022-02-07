import numpy as np

from kep import impacts, coords, grad_impacts
from phot import flux

__all__ = ['BarycenterOrbit', 'MoonOrbit', 'System']
    
class BarycenterOrbit:
    
    """
    The orbit of the planet-moon barycenter around the star. 
    
    Args:
        ab: Semimajor axis
        tb: Time of periastron passage
        eb: Eccentricity
        pb: Period
        wb: Argument of periastron (in radians)
        ib: Inclination (in radians)
    """
    
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
    
    """
    The orbit of the moon around the planet.
    
    Args:
        am: Semimajor axis
        tm: Time of periastron passage
        em: Eccentricity
        pm: Period
        om: Longitude of ascending node (in radians)
        wm: Argument of periastron (in radians)
        im: Inclination (in radians)
        mm: Moon/planet mass ratio
    """
    
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
    
    """
    Class representing the three-body star-planet-moon system. 
    
    
    .. note:: The dynamics of the system are approximated by assuming 
    no gravitational interaction between the moon and the star. 
    In other words, the planet/moon system is treated as a single 
    body for purposes of computing the orbit of their mutual 
    barycenter around the star, and then the orbit of the moon 
    around the planet is computed without considering star-moon 
    interactions. This apprximation should be valid in the limit 
    that the moon's orbit is well within the planet's Hill sphere. 
    
    Args:
        bo (BarycenterOrbit): The orbit of the planet-moon 
                barycenter with respect to the star.
        mo (MoonOrbit): The orbit of the moon with respect to the 
                planet. 
    """
    
    def __init__(self, bo, mo):
        
        self.bo = bo
        self.mo = mo
        self._lc = None
        self._grad = {}
        self._primary_coords = None
        self._secondary_coords = None
        self._computed = False
        self._computed_coords = False
        
    def coords(self, t):
        
        """
        Get the coordinates of the planet and moon.
        
        Args:
            t: Times at which the coordinates should be computed.
            
        Returns:
            pc: Coorinates of the planet as a tuple of arrays (x, y, z)
            mc: Coordinates of the moon as a tuple of arrays (x, y, z)
        """
        
        bo = self.bo
        mo = self.mo
        
        xp, yp, zp, xm, ym, zm = coords(t, {**bo.pdict, **mo.pdict})
        return (xp, yp, zp), (xm, ym, zm)
    
    def lightcurve(self, t, u1, u2, rp, rm, grad=False):
        
        """
        Get the lightcurve resulting from a transit of the moon/planet system.
        
        Args: 
            t: Times at which the flux should be computed
            rp: Radius of the planet
            rm: Radius of the moon
            u1: The first limb-darkening parameter
            u2: The second limb-darkening parameter
            grad (bool): If True, compute the gradient of the lightcurve.
                Default is False. 
                
        Returns:
            lc: The lightcurve
            grad: A dict containing the derivatives of the 
                lightcurve with respect to each of the input parameters.
            
        """
        
        bo = self.bo
        mo = self.mo
        
        if grad:
            bp, bpm, theta, dbp, dbpm, dtheta = grad_impacts(
                t, 
                {**bo.pdict, **mo.pdict}
            )
            lc = flux(
                u1, 
                u2, 
                rp, 
                rm, 
                bp, 
                bpm, 
                np.cos(theta), 
                np.sin(theta)
            ).T
                
            self._lc = lc[0]
            f_bp = lc[3]
            f_bpm = lc[4]
            f_theta = lc[5]
                
            df = (
                f_bp * dbp 
                + f_bpm * dbpm 
                + f_theta * dtheta
            )
                
            self._grad = {
                (bo.pnames + mo.pnames)[i]: df[i] 
                for i in range(np.shape(df)[0])
            }
                
            self._grad['rm'] = lc[1]
            self._grad['rp'] = lc[2]
            self._grad['u1'] = lc[6]
            self._grad['u2'] = lc[7]
            
            return self._lc, self._grad
                
        else:
            bp, bpm, theta = impacts(
                t, 
                {**bo.pdict, **mo.pdict}
            )
            self._lc = flux(
                u1, 
                u2, 
                rp, 
                rm, 
                bp, 
                bpm, 
                np.cos(theta), 
                np.sin(theta)
            )[:, 0]
            
            return self._lc
            
    def loglike(self, y, t, u1, u2, rp, rm, sigma):
        
        """
        Get the log-likelihood of the lightcurve.
        
        Args:
            y: A vector of observations to compute the likelihood with 
                respect to. 
            t: Times at which the flux should be computed
            rp: Radius of the planet
            rm: Radius of the moon
            u1: The first limb-darkening parameter
            u2: The second limb-darkening parameter
            sigma: The standard deviation of the model
        """
        
        mu = self.lightcurve(t, u1, u2, rp, rm)
        s2 = sigma * sigma
        return -0.5 * np.sum((y - mu) ** 2 / s2 + np.log(s2))