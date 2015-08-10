#!/usr/bin/env python2.7

from core import Dipole, DipolePack, World, Event
import numpy as np
from mayavi import mlab
from tvtk.tools import visual
Qel = 3e-7

d1 = Dipole(np.array([1e-2,1e-3,0.]),np.array([1e-2,0.,0.]),Qel,1e-2,1e-2)
d2 = Dipole(np.array([0.,0.,0.]),np.array([-1e-3,0.,0.]),Qel,1e-2,1e-2)
d3 = Dipole(np.array([-1e-2,2e-3,0.]),np.array([-1e-2,1e-3,0.]),Qel,1e-2,1e-2)
d4 = Dipole(np.array([0.,2e-3,-2e-2]),np.array([0.,3e-3,-2e-2]),Qel,1e-2,1e-2)
d5 = Dipole(np.array([0.,-2e-3,2e-2]),np.array([0.,-3e-3,2e-2]),Qel,1e-2,1e-2)
d6 = Dipole(np.array([-1e-3,-4e-3,3e-2]),np.array([-1e-3,-3e-3,3e-2]),Qel,1e-2,1e-2)

p1 = DipolePack(d1,d2,d3,d4)

w1 = World(200,200,200)

t = 10 # czas symulacji (zycie swiata)
dt = 0.0001 # skok sym.

e1 = Event( t, dt, p1, w1 )
e1.setoutput("output.dat")  

e1.run(True)
