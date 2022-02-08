import pytest
import numpy as np

import sys
sys.path.append('../gefera')
import testing
from kep import coords

np.random.seed(42)

params = [
    'ab', 'tb', 'eb', 'pb', 'wb', 'ib',      
    'am', 'tm', 'em', 'pm', 'om', 'wm', 'im', 'mm'
]

def test_kepler():

    n = 100
    M = np.random.rand(n) * np.pi * 2
    ecc = np.random.rand(n)
    ecc[-1] = 0.99999999
    ecc[-2] = 0.0

    testing.RPP_vs_newton(M, ecc)
    
def test_gradient():
    
    for _ in range(10):
        map(
            lambda p: testing.gradient(p), 
            params
        )
        
def test_coords():
    
    argdict = testing.random_args()
    argdict['ab'] = 0.0
    am = argdict['am']
    argdict['eb'] = 0.0
    argdict['em'] = 0.0
    
    t = np.linspace(0, 100, 10)
    xp, yp, zp, xm, ym, zm = coords(t, argdict)
    assert np.all(
        np.isclose(
            (xm - xp)**2 + (ym - yp)**2 + (zm - zp)**2, 
            am ** 2, 
            atol=1e-10
        )
    )
    
    argdict = testing.random_args()
    argdict['am'] = 0.0
    ab = argdict['ab']
    argdict['eb'] = 0.0
    argdict['em'] = 0.0
    
    t = np.linspace(0, 100, 10)
    xp, yp, zp, xm, ym, zm = coords(t, argdict)
    assert np.all(
        np.isclose(
            xm**2 + ym**2 + zm**2, 
            ab ** 2, 
            atol=1e-10
        )
    )
    assert np.all(
        np.isclose(
            xp**2 + yp**2 + zp**2, 
            ab ** 2, 
            atol=1e-10
        )
    )
    