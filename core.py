#!/usr/bin/env python2.7

import numpy as np
import math

from VisEngine import *
from tvtk.tools import visual

class Dipole:
	''' Zmienne do inicjalizacji:
		r1 - wektor wodzacy dla ladunku + [m] typu numpy.array()
		r2 - wektor wodzacy dla ladunku - [m] typu numpy.array()
		q - ladunek [C] typu float()
		m1 - masa punktowa posiadajaca ladunek + [kg] typu float()
		m2 - masa punktowa posiadajaca ladunek - [kg] typu float()
			  
	    Zmienne do obliczen:
		l - odleglosc miedzy ladunkami
		rCM - wektor wodzacy srodka masy
		mi - wektor momentu dipolowego
		v - vektor predkosci
		o - vektor predkosci katowej
		E1 - vektor natezenia pola w punkcie dla ladunku +
		E2 - vektor natezenia pola w punkcie dla ladunku -
	'''
	def __init__(self,r1,r2,q,m1,m2):
		self.r1 = r1
		self.r2 = r2
		self.q = q
		self.m1 = m1
		self.m2 = m2
		self.l = np.sqrt(np.subtract(self.r1, self.r2).dot(np.subtract(self.r1, self.r2))) # dlugosc to np.sqrt(x.dot(x))
		self.rCM = np.divide(np.add(np.multiply(self.r1,self.m1),np.multiply(self.r2,self.m2)),np.add(self.m1,self.m2))
		self.mi = np.multiply(self.q,np.subtract(self.r2,self.rCM))
		self.v = np.array([0.,0.,0.]) # predkosc postepowa
		self.o = np.array([0.,0.,0.]) # predkosc katowa
		self.I = GetInertia(r1,r2,m1,m2,self.rCM)	#moment bezwladnosci dipola 
		self.E1 = np.array([0.,0.,0.])
		self.E2 = np.array([0.,0.,0.])
		self.q1 = q
		self.q2 = -q
	def __str__(self):
		return  "____ ...:::: DIPOL ::::... ____"+"\n"+"r1 = "+str(self.r1)+"\n"+"r2 = "+str(self.r2)+"\n"+"q = "+str(self.q)+"\n"+"m1 = "+str(self.m1)+"\n"+"m2 = "+str(self.m2)+"\n"+"-------------------------------"

#################################################################################################################

class DipolePack:
	''' Paczka dipoli - argumenty typu "Dipole" '''
	def __init__(self,*items):
		self.items = [ i for i in items ]
	def unpack(self):
		''' zwraca liste z dipolami '''
		return self.items
	def append(self,*items):
		''' dodaje obiekty typu Dipole do listy '''
		for i in range(len(items)): self.items.append(items[i])
	def clean(self):
		''' czysci paczke dipoli'''
		self.items=[]
	
#################################################################################################################

class World:
	'''
	Jesli dipol wyjdzie poza swiat, konczy sie jego zycie
	Swiat jest prostopadloscianem x,y,z o srodku wyznaczonym przez wektor r
		x - dlugosc polowkowa w x [m] typu float()
		y - dlugosc polowkowa w y [m] typu float()
		z - dlugosc polowkowa w z [m] typu float()
		r - wektor wodzacy centrum ([m])  typu numpy.array()
	'''
	def __init__(self,x,y,z, r = np.array([0.,0.,0.])):
		self.x = x
		self.y = y
		self.z = z
		self.r = r

#################################################################################################################
def GetInertia(r1,r2,m1,m2,rCM):
	''' Zwraca moment bezwladnosci dipola '''
	R1 = r1-rCM	#wektor w ukl. srodka ciezkosci
	R2 = r2-rCM
	i1 = R1.dot(R1)*m1
	i2 = R2.dot(R2)*m2
	return i1+i2


#################################################################################################################

def GetElectricField(x,D):
	''' Zwraca wektor natezenia pola el. pochodzacy od dipola D w pkt wyznaczonym przez wektor x
		x - obiekt numpy.array()
		D - obiekt Dipole'''

	k = 9e10 # const.

	r1 = np.subtract(x,D.r1) # wektor od ladunku + dipola D do pkt x
	dl_r1 = np.sqrt(r1.dot(r1)) # dl wektora r1

	r2 = np.subtract(x,D.r2) # wektor od ladunku - dipola D do pkt x
	dl_r2 = np.sqrt(r2.dot(r2)) # dl wektora r2

	e1 = np.multiply(k*D.q1/dl_r1**3,r1) # wektor natezenia pochodzacy od ladunku + dipola D
	e2 = np.multiply(k*D.q2/dl_r2**3,r2) # wektor natezenia pochodzacy od ladunku - dipola D


	return np.add(e1,e2)

#################################################################################################################


def GetRotateProduct(vector,theta,axis=np.array([0,0,1])):
	''' Obraca vector (np.array()) o kat theta (float()) wokol osi 'axis' (np.array()) zgodnie z formula Eulera-Rodrigues'a i zasada prawoskretnosci '''
	axis = np.divide(axis,math.sqrt(np.dot(axis,axis)))
	a = math.cos(theta/2)
	b,c,d = -axis*math.sin(theta/2)
	rot_matrix = np.array([	[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
					[2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
					[2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]	])

	return np.dot(vector,rot_matrix)

#################################################################################################################
def GetParaVector(v1,v2):
	''' Zwraca rzut wektora v1 na wektor v2 '''
	return np.dot(np.dot(v1,v2)/np.dot(v2,v2),v2 )

#################################################################################################################
def GetOrthoVector(v1,v2):
	''' Zwraca skladowa prostopadla wektora v1 do wektora v2 '''
	return v1-np.dot(np.dot(v1,v2)/np.dot(v2,v2),v2 )

#################################################################################################################
def LetElectricForceWork(D,t):
	''' Symuluje dzialanie sily na dipol D w czasie t
		D - obiekt typu Dipole '''

	#_____ruch postepowy_____#
	F1para = GetParaVector(np.dot(D.q1,D.E1),D.mi)
	F2para = GetParaVector(np.dot(D.q2,D.E2),D.mi)
	Fwpara = F1para+F2para
	a = Fwpara/(D.m1+D.m2)

	v = a*t
	D.v += v
#	print(D.v)
	r = D.v*t

	D.r1 += r
	D.r2 += r

	D.rCM = np.divide(np.add(np.multiply(D.r1,D.m1),np.multiply(D.r2,D.m2)),np.add(D.m1,D.m2))


	#_____ruch obrotowy______#
	F1ortho = GetOrthoVector(np.dot(D.q1,D.E1),D.mi)
	F2ortho = GetOrthoVector(np.dot(D.q2,D.E2),D.mi)

	R1 = D.r1-D.rCM	#wektor w ukl. srodka ciezkosci
	R2 = D.r2-D.rCM

	#dzialanie sil ortho w jednym kierunkju
	if np.sqrt((F1ortho+F2ortho).dot((F1ortho+F2ortho))) > np.sqrt(F1ortho.dot(F1ortho)) and np.sqrt((F1ortho+F2ortho).dot((F1ortho+F2ortho))) > np.sqrt(F2ortho.dot(F2ortho)) :
		if np.sqrt(F1ortho.dot(F1ortho)) > np.sqrt(F2ortho.dot(F2ortho)):
			a = 2.*F2ortho/(D.m1+D.m2)	
			F2ortho=0*F2ortho
			F1ortho=F1ortho-F2ortho
			
		else:
			a = 2.*F1ortho/(D.m1+D.m2)	
			F1ortho=0*F1ortho
			F2ortho=F2ortho-F1ortho
	
		v = a*t
		D.v += v

		r = D.v*t

		D.r1 += r
		D.r2 += r

		D.rCM = np.divide(np.add(np.multiply(D.r1,D.m1),np.multiply(D.r2,D.m2)),np.add(D.m1,D.m2))	





	M1 = np.cross(R1,F1ortho)	# moment sily
	M2 = np.cross(R2,F2ortho)
	Mw = M1+M2 	# wypadkowy moment sily

	epsilon = Mw/D.I

	omega = epsilon*t 
	D.o += omega

	alpha = t*np.sqrt(D.o.dot(D.o)) # kat to dlugosc wektora omega razy maluskie dt
	if alpha != 0:
		R1 = GetRotateProduct(R1,alpha,omega) # obrot wektorow w ukl srodka masy dipola
		R2 = GetRotateProduct(R2,alpha,omega)

		D.r1 = D.rCM + R1 # powrot do wspolrzednych globalnych
		D.r2 = D.rCM + R2

		#D.l = np.sqrt(np.subtract(D.r1, D.r2).dot(np.subtract(D.r1, D.r2))); print(D.l) # sprawdzenie czy dipol sie nie rozjezdza
		D.mi = np.multiply(D.q,np.subtract(D.r2,D.rCM)) # aktualizacja momentu dipolowego

#################################################################################################################
def IsInWorld(D,W):
	''' Zwraca zmienna typu bool okreslajaca czy dipol D jest w swiecie W '''
	q1_is_in = True if D.r1[0]<(W.x+W.r[0]) and D.r1[1]<(W.y+W.r[1])  and D.r1[2]<(W.z+W.r[2]) else False
	q2_is_in = True if D.r2[0]<(W.x+W.r[0]) and D.r2[1]<(W.y+W.r[1])  and D.r2[2]<(W.z+W.r[2]) else False
	#print(q1_is_in,q2_is_in)
	return q1_is_in and q2_is_in

#################################################################################################################
class Event:
	''' Uruchamia symulacje na czas 't' [s], ze skokiem 'dt' [s] dla dipoli 'dpack' w swiecie 'world' '''
	def __init__(self,t,dt,dpack,world):
		self.t = t			# czas symulacji [s]
		self.dt = dt			# skok [s]
		self.dlist = dpack.unpack() 	# lista z dipolami
		self.world = world		# obiekt typu World
		self.time = 0			# zegar symulacji
		self.fname = None 		# nazwa pliku do zrzutu danych
		self.fobj = None		# zmienna pliku do zapisu

	def __step(self,vis=False):
		## liczenie natezenia w kazdym punkcie
		for to_mod_dipole in self.dlist:
			to_mod_dipole.E1=np.array([0.,0.,0.])
			to_mod_dipole.E2=np.array([0.,0.,0.])
			for get_dipole in self.dlist:
				if not to_mod_dipole == get_dipole:
					to_mod_dipole.E1=np.add(to_mod_dipole.E1,GetElectricField(to_mod_dipole.r1,get_dipole))
					to_mod_dipole.E2=np.add(to_mod_dipole.E2,GetElectricField(to_mod_dipole.r2,get_dipole))
					
		## poddajemy dipole dzialaniom sil 
		for i in range(len(self.dlist)):
			if IsInWorld(self.dlist[i],self.world):		# sprawdzenie czy dipol jest w swiecie
				LetElectricForceWork(self.dlist[i],self.dt)
			else:
				raise 0	# jesli wyszedl poza swiat - koniec sym.

		## wizualizacja stepu
		if vis: 
			print "Time: "+str(self.time)+"  Refreshing..."
			self.visdata.NextFrame(self.dlist)



		## zapis danych do pliku 
		self.fobj.write(str(self.time)+"\n")

		## END ##


	def run(self,vis=False):
		''' Domyslnie vizualizacje sa wylaczone - vis = True uruchamia vizualizacje '''
		#zabezpieczenie przed overloadingiem , zero division itp
		np.seterr(all="raise") # w razie przekroczenia zakresu wyrzuci yjatek
		if vis: 
			f = mlab.figure(size=(600,600))
			visual.set_viewer(f)
			self.visdata = VisData(self.dlist)# lista punktow do wyswietlenia - jesli vis = True bedzie aktualizowana i wyswietlana

		if not self.fname == None:
			self.fobj = open(self.fname,'w')

		while self.time<self.t:
			try:
				self.__step(vis) # jesli koniec rzuci wyjatek
				self.time+=self.dt
			except: 
				print("Jeden z dipoli opuscil swiat - koniec symulacji")
				return 1

		if self.fobj:
			self.fobj.close()

	def setoutput(self,fname):		# nazwa pliku w postaci str
		self.fname = fname


#################################################################################################################



