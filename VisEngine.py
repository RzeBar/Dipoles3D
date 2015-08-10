#!/usr/bin/env python2.7
import numpy as np
from mayavi.mlab import *
from mayavi import mlab
from tvtk.tools import visual


def Arrow_From_A_to_B(x1, y1, z1, x2, y2, z2):
    ar1=visual.arrow(x=x1, y=y1, z=z1)
    ar1.length_cone=0.4

    arrow_length=np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
    ar1.actor.scale=[arrow_length, arrow_length, arrow_length]
    ar1.pos = ar1.pos/arrow_length
    ar1.axis = [x2-x1, y2-y1, z2-z1]
    return ar1

class VisData:
	''' Wizualizuje symulacje '''
	def __init__(self,Dlist): # czerwony + , niebieski -
		self.VisPoints = [ [visual.sphere(x=D.r1[0],y=D.r1[1],z=D.r1[2],color=visual.color.red,radius=D.l/2**0.5),visual.sphere(x=D.r2[0],y=D.r2[1],z=D.r2[2],color=visual.color.blue,radius=D.l/2**0.5)] for D in Dlist]

	
	def NextFrame(self,Dlist):
		for i in range(len(Dlist)):
			self.VisPoints[i][0].pos = Dlist[i].r1
			self.VisPoints[i][1].pos = Dlist[i].r2


	
