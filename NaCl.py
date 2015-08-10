#!/usr/bin/env python2.7
##########################################
#
#   Proba symulacji krysztalu NaCL:
#
#   totalna klapa - model sie nie nadaje
#   jest wolny i niepoprawny fizycznie
#   (nie uwzglednia wiazan elektronowych)
#   
##########################################
from core import Dipole, DipolePack, World, Event
import numpy as np
from mayavi import mlab
from tvtk.tools import visual

# The cell-edge length is 564.1 pm. So, Na<-- 282.05pm -->Cl

Q = 1.6e-17 # [C]
mNa = 1 # [kg]
mCl = 1 # [kg]
d = 282.05e-12 # bond length [m]
 
NaCl = []

# ilosc dipoli w wymiarze
x =  4
y =  4
z =  4  #ilosc warstw

t = d #translacja 
# 1:Na , 2:Cl
for i in range(0,z+1):
	for j in range(0,y+1):
		for k in range(0,x+1):
			if i%2==0:
				if j%2 == 0: NaCl.append(Dipole(   np.array([k*d*2+d+t,j*d*2+t,i*d*2+t])  ,  np.array([k*d*2+t,j*d*2+t,i*d*2+t]) , Q , mNa, mCl  ))
				else:        NaCl.append(Dipole(   np.array([k*d*2+t,j*d*2+t,i*d*2+t])    ,  np.array([k*d*2+d+t,j*d*2+t,i*d*2+t]) , Q , mNa, mCl  ))
			else:
				if j%2 == 0: NaCl.append(Dipole(   np.array([k*d*2+t,j*d*2+t,i*d*2+t])    ,  np.array([k*d*2+d+t,j*d*2+t,i*d*2+t]) , Q , mNa, mCl  ))
				else:        NaCl.append(Dipole(   np.array([k*d*2+d+t,j*d*2+t,i*d*2+t])  ,  np.array([k*d*2+t,j*d*2+t,i*d*2+t]) , Q , mNa, mCl  ))    
p1 = DipolePack(*NaCl)

w1 = World(200,200,200)

time = 10 # czas symulacji (zycie swiata)
dt = 0.0001 # skok sym.

e1 = Event( time, dt, p1, w1 )
e1.setoutput("output.dat")  

e1.run() # True - wizulaizuje, False - bez wizyalizacji
