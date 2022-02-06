import aesara
from .ops import ImpactsOp, FluxOp

__args__ = [
    'ap', 
    'tp', 
    'ep', 
    'pp', 
    'wp', 
    'ip',     
    'am', 
    'tm', 
    'em', 
    'pm', 
    'om', 
    'wm', 
    'im', 
    'mm'
]

def validate_elements(argdict):
    
    if set(argdict.keys()) == set(__args__):
        
        idm = {
            v: i 
            for i, v in enumerate(__args__)
        }
        
        args = tuple(
            dict(
                sorted(
                    argdict.items(), 
                    key=lambda p: idm[p[0]]
                )
            ).values()
        )
        return args
        
    else:
        raise ValueError(
            "required parameters are: ", 
            __args__
        )

class star:
        
    def __init__(self, u1, u2):
        
        self.u1 = u1
        self.u2 = u2
    
class planet:
    
    def __init__(self, rp, ap, tp, ep, pp, wp, ip):
        
        self.r = rp
        self.pnames = [
            'ap', 
            'tp', 
            'ep', 
            'pp', 
            'wp', 
            'ip'
        ]
        self.pdict = dict(
            zip(
                self.pnames,
                (ap, tp, ep, pp, wp, ip)
            )
        )
    
class moon:
    
    def __init__(self, rm, am, tm, em, pm, om, wm, im, mm):
        
        self.r = rm
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
        
class system:
    
    def __init__(self, star, planet, moon):
        
        self.star = star
        self.planet = planet
        self.moon = moon
        
    def light_curve(self, t):
        
        impacts = ImpactsOp(t)
        flux = FluxOp()
        
        impact_params = validate_elements(
            {
                **self.planet.pdict, 
                **self.moon.pdict
            }
        )
        
        out = impacts(*impact_params)
        bp, bpm, theta = out[0], out[1], out[2]
        return flux(
            self.star.u1, 
            self.star.u2, 
            self.planet.r, 
            self.moon.r, 
            bp, bpm, theta
        )[0]