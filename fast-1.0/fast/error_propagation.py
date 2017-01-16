# -*- coding: utf-8 -*-

#************************************************************************
#       Copyright (C) 2014 - 2017 Oscar Gerardo Lazo Arjona             *
#              <oscar.lazo@correo.nucleares.unam.mx>                    *
#                                                                       *
#  This file is part of FAST.                                           *
#                                                                       *
#  FAST is free software: you can redistribute it and/or modify         *
#  it under the terms of the GNU General Public License as published by *
#  the Free Software Foundation, either version 3 of the License, or    *
#  (at your option) any later version.                                  *
#                                                                       *
#  FAST is distributed in the hope that it will be useful,              *
#  but WITHOUT ANY WARRANTY; without even the implied warranty of       *
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
#  GNU General Public License for more details.                         *
#                                                                       *
#  You should have received a copy of the GNU General Public License    *
#  along with FAST.  If not, see <http://www.gnu.org/licenses/>.        *
#                                                                       *
#************************************************************************

#The class Measurement will be defined as an extention of floats
#that have associated standard deviations or errors (x,sigmax).
#Basic operations on these numbers will be defined such that one
#can comfortably make operations on these numbers to obtain the
#corresponding errors.

from math import sqrt, log

class Measurement(object):
	def __init__(self,value,sigma):
		self.value=float(value)
		self.sigma=sigma
	
	def __str__(self):
		return '('+str(self.value)+', '+str(self.sigma)+')'
	
	#Multiplication.
	def __mul__(self,other,cov=0.0):
		#Scalar multiplication
		if isinstance(other,float) or isinstance(other,int):
			return Measurement(other*self.value, abs(other)*self.sigma)
		#Measurement multiplication
		elif isinstance(other,Measurement):
			sigmaf = self.value**2 * other.sigma**2
			sigmaf+= other.value**2 * self.sigma**2
			sigmaf+= 2*self.value*other.value*cov
			sigmaf = sqrt(sigmaf)
			return Measurement(self.value*other.value,sigmaf)			
	def __rmul__(self,other):
		return self.__mul__(other)
	
	#Addition.
	def __add__(self,other,cov=0.0):
		#Scalar addition
		if isinstance(other,float) or isinstance(other,int):
			return Measurement(other+self.value, self.sigma)
		#Measurement addition
		elif isinstance(other,Measurement):
			sigmaf = self.sigma**2 + other.sigma**2 + 2*cov
			sigmaf = sqrt(sigmaf)
			return Measurement(self.value + other.value,sigmaf)			
	def __radd__(self,other):
		return self.__add__(other)

	#Substraction.
	def __sub__(self,other,cov=0.0):
		#Scalar substraction
		if isinstance(other,float) or isinstance(other,int):
			return Measurement(-other+self.value, self.sigma)
		#Measurement substraction
		elif isinstance(other,Measurement):
			sigmaf = self.sigma**2 + other.sigma**2 - 2*cov
			sigmaf = sqrt(sigmaf)
			return Measurement(self.value - other.value,sigmaf)			
	def __rsub__(self,other):
		if isinstance(other,float) or isinstance(other,int):
			other=Measurement(other,0.0)

		return other.__sub__(self)

	#Division.
	def __div__(self,other,cov=0.0):
		#Scalar division.
		if isinstance(other,float) or isinstance(other,int):
			other=Measurement(other,0.0)
		#Measurement division.
		sigmaf = (self.sigma/self.value)**2 
		sigmaf+= (other.sigma/other.value)**2 - 2*cov/(self.value*other.value)
		sigmaf = sqrt(sigmaf)

		return Measurement(self.value / other.value,sigmaf)			
	def __rdiv__(self,other):
		if isinstance(other,float) or isinstance(other,int):
			other=Measurement(other,0.0)
		
		return other.__div__(self)
	
	#Negative.
	def __neg__(self):
		return Measurement(-self.value,self.sigma)
	
	#Power.
	def __pow__(self,other,cov=0.0):
		#Scalar power.
		if isinstance(other,float) or isinstance(other,int):
			other=Measurement(other,0.0)
		#Measurement power.
		sigmaf = (other.value*self.sigma/self.value)**2 
		sigmaf+= (log(self.value)*other.sigma)**2
		sigmaf+= 2*other.value*log(self.value)*cov/self.value
		sigmaf = sqrt(sigmaf)
		
		return Measurement(self.value ** other.value , sigmaf)
	def __rpow__(self,other):
		if isinstance(other,float) or isinstance(other,int):
			other=Measurement(other,0.0)
		
		return other.__pow__(self)

def rel(m):
	return m.sigma/m.value
P=100 #uW
err=0.2

P=Measurement(P,P*err)
a=Measurement(4.7,0.0)
k=5

E02=k*P/a**2
E0=(k*P/a**2)**0.5

#print a,rel(a)
#print a**2,rel(a**2)
#print 1/a**2,rel(1/a**2)

#~ print P,P.sigma/P.value
#~ print E02,E02.sigma/E02.value
print E02,rel(E02)
print E0, rel(E0)
