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

sage_included = 'sage' in globals().keys()

#Physical constants (SI units):
from scipy.constants import physical_constants
from scipy.constants import pi as Pi # Every one's favourite constant
#                                                         units
c   =physical_constants["speed of light in vacuum"][0]  # m/s
a0  =physical_constants["Bohr radius"][0]               # m
hbar=physical_constants["Planck constant over 2 pi"][0] # J s
e   =physical_constants["elementary charge"][0]         # C
mu0 =physical_constants["mag. constant"][0]             # N / A^2

if not sage_included:
    from math import pi,sqrt,cos,sin
    Pi=pi
else:
    Pi=pi.n()

class PlaneWave(object):
    r"""This class implements plane waves propagating in an arbitrary direction and with an arbitrary (and well-defined) polarization. It takes as input:
    
    phi .- The azimutal angle in spherical coordinates of the wave vector.
    theta .- The polar angle in spherical coordinates of the wave vector.
    alpha .- The angle between the fast axis of a half-wave plate and 
             an incident linearly polarized beam.
    beta .- The angle between the fast axis of a quarter-wave plate and 
             an incident linearly polarized beam.
    
    The easiest way to understand what this means is to use the function
    draw_lasers to get a 3d view of the direction and polarization of the
    wave.
    
    An explicit way to understand what alpha and beta are is to see how the
    overall form of the electric field is contructed. We begin with a
    polarization vector for a plane wave propagating towards z:
    
    >>> from sympy import symbols, Matrix, I, exp, sin, cos, pprint
    >>> E0=symbols("E0",real=True)
    >>> t,x,y,z,omega,kx,ky,kz=symbols("t,x,y,z,omega,kx,ky,kz",real=True)
    >>> phi, theta, alpha, beta = symbols("phi, theta, alpha, beta",real=True)
    >>> ep=Matrix([cos(2*beta),I*sin(2*beta),0])

    >>> pprint(ep,use_unicode=False)
    [ cos(2*beta) ]
    [             ]
    [I*sin(2*beta)]
    [             ]
    [      0      ]
    
    Where beta specifies the circularity of the polarization. We define 
    the following rotation matrices
    
    >>> R1=Matrix([[cos(2*alpha),-sin(2*alpha),0],[sin(2*alpha),cos(2*alpha),0],[0,0,1]])
    >>> pprint(R1,use_unicode=False)
    [cos(2*alpha)  -sin(2*alpha)  0]
    [                              ]
    [sin(2*alpha)  cos(2*alpha)   0]
    [                              ]
    [     0              0        1]
    
    >>> R2=Matrix([[cos(theta),0, sin(theta)],[0,1,0],[-sin(theta),0,cos(theta)]])
    >>> pprint(R2,use_unicode=False)
    [cos(theta)   0  sin(theta)]
    [                          ]
    [     0       1      0     ]
    [                          ]
    [-sin(theta)  0  cos(theta)]
    
    >>> R3=Matrix([[cos(phi),-sin(theta),0],[ sin(theta),cos(theta),0],[0,0,1]])
    >>> pprint(R3,use_unicode=False)
    [ cos(phi)   -sin(theta)  0]
    [                          ]
    [sin(theta)  cos(theta)   0]
    [                          ]
    [    0            0       1]
    
    Then we create a general polarization vector applying this rotation to ep
    
    >>> epsilonp=R3*R2*R1*ep
    
    And we build the plane wave in the following way
    
    >>> plane_wave=E0/2 * exp(I*(kx*x+ky*y+kz*z-omega*t)) * epsilonp
    >>> plane_wave=plane_wave + plane_wave.conjugate()
    
    """
    def __init__(self,phi,theta,alpha,beta,omega=1,E0=1,color='blue',symbolical=True):
        self.phi=phi
        self.theta=theta
        self.alpha=alpha
        self.beta=beta
        self.E0=E0
        self.omega=omega
        if sage_included:
            self.color=color
        else:
            if color=='blue':
                self.color='b'
            elif color=='red':
                self.color='r'
            elif color=='green':
                self.color='g'

        
        if sage_included:
			t=var('t',domain=RR)
			#We calculate the electric field
			E=[(cos(alpha)*cos(phi)*cos(theta) + sin(alpha)*sin(phi))*E0*cos(omega*t + beta) -
			 (cos(phi)*cos(theta)*sin(alpha) - cos(alpha)*sin(phi))*E0*cos(omega*t)]
			E.append((cos(alpha)*cos(theta)*sin(phi) - cos(phi)*sin(alpha))*E0*cos(omega*t + beta) -
			 (cos(theta)*sin(alpha)*sin(phi) + cos(alpha)*cos(phi))*E0*cos(omega*t))
			E.append(-E0*cos(omega*t + beta)*cos(alpha)*sin(theta) + E0*cos(omega*t)*sin(alpha)*sin(theta))
			self.E=E
                
        #We calculate the laser's cofficients Y^l(+)_p and Y^l(-)_p
        if sage_included:
            Ypm1=1/Integer(4)*sqrt(Integer(2))*(cos(alpha)*cos(beta)*cos(phi)*cos(theta) - I*cos(alpha)*cos(phi)*cos(theta)*sin(beta) - I*cos(alpha)*cos(beta)*cos(theta)*sin(phi) - cos(alpha)*cos(theta)*sin(beta)*sin(phi) + I*cos(beta)*cos(phi)*sin(alpha) - cos(phi)*cos(theta)*sin(alpha) + cos(phi)*sin(alpha)*sin(beta) + cos(beta)*sin(alpha)*sin(phi) + I*cos(theta)*sin(alpha)*sin(phi) - I*sin(alpha)*sin(beta)*sin(phi) + I*cos(alpha)*cos(phi) + cos(alpha)*sin(phi))

            Yp0 =-1/Integer(2)*(cos(alpha)*cos(beta) - I*cos(alpha)*sin(beta) - sin(alpha))*sin(theta)

            Yp1 =-1/Integer(4)*sqrt(Integer(2))*(cos(alpha)*cos(beta)*cos(phi)*cos(theta) - I*cos(alpha)*cos(phi)*cos(theta)*sin(beta) + I*cos(alpha)*cos(beta)*cos(theta)*sin(phi) + cos(alpha)*cos(theta)*sin(beta)*sin(phi) - I*cos(beta)*cos(phi)*sin(alpha) - cos(phi)*cos(theta)*sin(alpha) - cos(phi)*sin(alpha)*sin(beta) + cos(beta)*sin(alpha)*sin(phi) - I*cos(theta)*sin(alpha)*sin(phi) - I*sin(alpha)*sin(beta)*sin(phi) - I*cos(alpha)*cos(phi) + cos(alpha)*sin(phi))

            Ymm1=1/Integer(4)*sqrt(Integer(2))*(cos(alpha)*cos(beta)*cos(phi)*cos(theta) + I*cos(alpha)*cos(phi)*cos(theta)*sin(beta) - I*cos(alpha)*cos(beta)*cos(theta)*sin(phi) + cos(alpha)*cos(theta)*sin(beta)*sin(phi) + I*cos(beta)*cos(phi)*sin(alpha) - cos(phi)*cos(theta)*sin(alpha) - cos(phi)*sin(alpha)*sin(beta) + cos(beta)*sin(alpha)*sin(phi) + I*cos(theta)*sin(alpha)*sin(phi) + I*sin(alpha)*sin(beta)*sin(phi) + I*cos(alpha)*cos(phi) + cos(alpha)*sin(phi))

            Ym0 =-1/Integer(2)*(cos(alpha)*cos(beta) + I*cos(alpha)*sin(beta) - sin(alpha))*sin(theta)

            Ym1 =-1/Integer(4)*sqrt(Integer(2))*(cos(alpha)*cos(beta)*cos(phi)*cos(theta) + I*cos(alpha)*cos(phi)*cos(theta)*sin(beta) + I*cos(alpha)*cos(beta)*cos(theta)*sin(phi) - cos(alpha)*cos(theta)*sin(beta)*sin(phi) - I*cos(beta)*cos(phi)*sin(alpha) - cos(phi)*cos(theta)*sin(alpha) + cos(phi)*sin(alpha)*sin(beta) + cos(beta)*sin(alpha)*sin(phi) - I*cos(theta)*sin(alpha)*sin(phi) + I*sin(alpha)*sin(beta)*sin(phi) - I*cos(alpha)*cos(phi) + cos(alpha)*sin(phi))
            
            
            
            Xpm1= 1/Integer(2)*cos(alpha)*cos(beta)*cos(phi)*cos(theta) - 1/Integer(2)*I*cos(alpha)*cos(phi)*cos(theta)*sin(beta) - 1/Integer(2)*cos(phi)*cos(theta)*sin(alpha) + 1/Integer(2)*cos(beta)*sin(alpha)*sin(phi) - 1/Integer(2)*I*sin(alpha)*sin(beta)*sin(phi) + 1/Integer(2)*cos(alpha)*sin(phi)

            Xp0 = 1/Integer(2)*cos(alpha)*cos(beta)*cos(theta)*sin(phi) - 1/Integer(2)*I*cos(alpha)*cos(theta)*sin(beta)*sin(phi) - 1/Integer(2)*cos(beta)*cos(phi)*sin(alpha) + 1/Integer(2)*I*cos(phi)*sin(alpha)*sin(beta) - 1/Integer(2)*cos(theta)*sin(alpha)*sin(phi) - 1/Integer(2)*cos(alpha)*cos(phi)

            Xp1 =-1/Integer(2)*cos(alpha)*cos(beta)*sin(theta) + 1/Integer(2)*I*cos(alpha)*sin(beta)*sin(theta) + 1/Integer(2)*sin(alpha)*sin(theta)
            
            Xmm1= 1/Integer(2)*cos(alpha)*cos(beta)*cos(phi)*cos(theta) + 1/Integer(2)*I*cos(alpha)*cos(phi)*cos(theta)*sin(beta) - 1/Integer(2)*cos(phi)*cos(theta)*sin(alpha) + 1/Integer(2)*cos(beta)*sin(alpha)*sin(phi) + 1/Integer(2)*I*sin(alpha)*sin(beta)*sin(phi) + 1/Integer(2)*cos(alpha)*sin(phi)

            Xm0 = 1/Integer(2)*cos(alpha)*cos(beta)*cos(theta)*sin(phi) + 1/Integer(2)*I*cos(alpha)*cos(theta)*sin(beta)*sin(phi) - 1/Integer(2)*cos(beta)*cos(phi)*sin(alpha) - 1/Integer(2)*I*cos(phi)*sin(alpha)*sin(beta) - 1/Integer(2)*cos(theta)*sin(alpha)*sin(phi) - 1/Integer(2)*cos(alpha)*cos(phi)

            Xm1 =-1/Integer(2)*cos(alpha)*cos(beta)*sin(theta) - 1/Integer(2)*I*cos(alpha)*sin(beta)*sin(theta) + 1/Integer(2)*sin(alpha)*sin(theta)
            
            
        else:
            #~ Ypm1=0.25*sqrt(2)*(cos(alpha)*cos(beta)*cos(phi)*cos(theta) - 1j*cos(alpha)*cos(phi)*cos(theta)*sin(beta) - 1j*cos(alpha)*cos(beta)*cos(theta)*sin(phi) - cos(alpha)*cos(theta)*sin(beta)*sin(phi) + 1j*cos(beta)*cos(phi)*sin(alpha) - cos(phi)*cos(theta)*sin(alpha) + cos(phi)*sin(alpha)*sin(beta) + cos(beta)*sin(alpha)*sin(phi) + 1j*cos(theta)*sin(alpha)*sin(phi) - 1j*sin(alpha)*sin(beta)*sin(phi) + 1j*cos(alpha)*cos(phi) + cos(alpha)*sin(phi))
#~ 
            #~ Yp0 =-0.5*(cos(alpha)*cos(beta) - 1j*cos(alpha)*sin(beta) - sin(alpha))*sin(theta)
#~ 
            #~ Yp1 =-0.25*sqrt(2)*(cos(alpha)*cos(beta)*cos(phi)*cos(theta) - 1j*cos(alpha)*cos(phi)*cos(theta)*sin(beta) + 1j*cos(alpha)*cos(beta)*cos(theta)*sin(phi) + cos(alpha)*cos(theta)*sin(beta)*sin(phi) - 1j*cos(beta)*cos(phi)*sin(alpha) - cos(phi)*cos(theta)*sin(alpha) - cos(phi)*sin(alpha)*sin(beta) + cos(beta)*sin(alpha)*sin(phi) - 1j*cos(theta)*sin(alpha)*sin(phi) - 1j*sin(alpha)*sin(beta)*sin(phi) - 1j*cos(alpha)*cos(phi) + cos(alpha)*sin(phi))
#~ 
            #~ Ymm1=0.25*sqrt(2)*(cos(alpha)*cos(beta)*cos(phi)*cos(theta) + 1j*cos(alpha)*cos(phi)*cos(theta)*sin(beta) - 1j*cos(alpha)*cos(beta)*cos(theta)*sin(phi) + cos(alpha)*cos(theta)*sin(beta)*sin(phi) + 1j*cos(beta)*cos(phi)*sin(alpha) - cos(phi)*cos(theta)*sin(alpha) - cos(phi)*sin(alpha)*sin(beta) + cos(beta)*sin(alpha)*sin(phi) + 1j*cos(theta)*sin(alpha)*sin(phi) + 1j*sin(alpha)*sin(beta)*sin(phi) + 1j*cos(alpha)*cos(phi) + cos(alpha)*sin(phi))
#~ 
            #~ Ym0 =-0.5*(cos(alpha)*cos(beta) + 1j*cos(alpha)*sin(beta) - sin(alpha))*sin(theta)
#~ 
            #~ Ym1 =-0.25*sqrt(2)*(cos(alpha)*cos(beta)*cos(phi)*cos(theta) + 1j*cos(alpha)*cos(phi)*cos(theta)*sin(beta) + 1j*cos(alpha)*cos(beta)*cos(theta)*sin(phi) - cos(alpha)*cos(theta)*sin(beta)*sin(phi) - 1j*cos(beta)*cos(phi)*sin(alpha) - cos(phi)*cos(theta)*sin(alpha) + cos(phi)*sin(alpha)*sin(beta) + cos(beta)*sin(alpha)*sin(phi) - 1j*cos(theta)*sin(alpha)*sin(phi) + 1j*sin(alpha)*sin(beta)*sin(phi) - 1j*cos(alpha)*cos(phi) + cos(alpha)*sin(phi))

            ################################

            Ypm1 = 0.5*sqrt(2)*((cos(2*alpha)*cos(phi)*cos(theta) - sin(2*alpha)*sin(phi))*cos(2*beta) - 1j*(cos(phi)*cos(theta)*sin(2*alpha) + cos(2*alpha)*sin(phi))*sin(2*beta)) - 0.5*1j*sqrt(2)*((cos(2*alpha)*cos(theta)*sin(phi) + cos(phi)*sin(2*alpha))*cos(2*beta) - 1j*(cos(theta)*sin(2*alpha)*sin(phi) - cos(2*alpha)*cos(phi))*sin(2*beta))

            Yp0 = -cos(2*alpha)*cos(2*beta)*sin(theta) + 1j*sin(2*alpha)*sin(2*beta)*sin(theta)

            Yp1 = -0.5*sqrt(2)*((cos(2*alpha)*cos(phi)*cos(theta) - sin(2*alpha)*sin(phi))*cos(2*beta) - 1j*(cos(phi)*cos(theta)*sin(2*alpha) + cos(2*alpha)*sin(phi))*sin(2*beta)) - 0.5*1j*sqrt(2)*((cos(2*alpha)*cos(theta)*sin(phi) + cos(phi)*sin(2*alpha))*cos(2*beta) - 1j*(cos(theta)*sin(2*alpha)*sin(phi) - cos(2*alpha)*cos(phi))*sin(2*beta))


            Ymm1 = 0.5*sqrt(2)*((cos(2*alpha)*cos(phi)*cos(theta) - sin(2*alpha)*sin(phi))*cos(2*beta) + 1j*(cos(phi)*cos(theta)*sin(2*alpha) + cos(2*alpha)*sin(phi))*sin(2*beta)) - 0.5*1j*sqrt(2)*((cos(2*alpha)*cos(theta)*sin(phi) + cos(phi)*sin(2*alpha))*cos(2*beta) + 1j*(cos(theta)*sin(2*alpha)*sin(phi) - cos(2*alpha)*cos(phi))*sin(2*beta))

            Ym0 = -cos(2*alpha)*cos(2*beta)*sin(theta) - 1j*sin(2*alpha)*sin(2*beta)*sin(theta)

            Ym1 = -0.5*sqrt(2)*((cos(2*alpha)*cos(phi)*cos(theta) - sin(2*alpha)*sin(phi))*cos(2*beta) + 1j*(cos(phi)*cos(theta)*sin(2*alpha) + cos(2*alpha)*sin(phi))*sin(2*beta)) - 0.5*1j*sqrt(2)*((cos(2*alpha)*cos(theta)*sin(phi) + cos(phi)*sin(2*alpha))*cos(2*beta) + 1j*(cos(theta)*sin(2*alpha)*sin(phi) - cos(2*alpha)*cos(phi))*sin(2*beta))



            #I=1j

            #~ Xpm1= 0.5*cos(alpha)*cos(beta)*cos(phi)*cos(theta) - 0.5*1j*cos(alpha)*cos(phi)*cos(theta)*sin(beta) - 0.5*cos(phi)*cos(theta)*sin(alpha) + 0.5*cos(beta)*sin(alpha)*sin(phi) - 0.5*1j*sin(alpha)*sin(beta)*sin(phi) + 0.5*cos(alpha)*sin(phi)
#~ 
            #~ Xp0 = 0.5*cos(alpha)*cos(beta)*cos(theta)*sin(phi) - 0.5*1j*cos(alpha)*cos(theta)*sin(beta)*sin(phi) - 0.5*cos(beta)*cos(phi)*sin(alpha) + 0.5*1j*cos(phi)*sin(alpha)*sin(beta) - 0.5*cos(theta)*sin(alpha)*sin(phi) - 0.5*cos(alpha)*cos(phi)
#~ 
            #~ Xp1 =-0.5*cos(alpha)*cos(beta)*sin(theta) + 0.5*1j*cos(alpha)*sin(beta)*sin(theta) + 0.5*sin(alpha)*sin(theta)
            #~ 
            #~ Xmm1= 0.5*cos(alpha)*cos(beta)*cos(phi)*cos(theta) + 0.5*1j*cos(alpha)*cos(phi)*cos(theta)*sin(beta) - 0.5*cos(phi)*cos(theta)*sin(alpha) + 0.5*cos(beta)*sin(alpha)*sin(phi) + 0.5*1j*sin(alpha)*sin(beta)*sin(phi) + 0.5*cos(alpha)*sin(phi)
#~ 
            #~ Xm0 = 0.5*cos(alpha)*cos(beta)*cos(theta)*sin(phi) + 0.5*1j*cos(alpha)*cos(theta)*sin(beta)*sin(phi) - 0.5*cos(beta)*cos(phi)*sin(alpha) - 0.5*1j*cos(phi)*sin(alpha)*sin(beta) - 0.5*cos(theta)*sin(alpha)*sin(phi) - 0.5*cos(alpha)*cos(phi)
#~ 
            #~ Xm1 =-0.5*cos(alpha)*cos(beta)*sin(theta) - 0.5*1j*cos(alpha)*sin(beta)*sin(theta) + 0.5*sin(alpha)*sin(theta)
			

        if not symbolical and sage_included:

			#~ Xpm1=Xpm1.N()
			#~ Xp0 = Xp0.N()
			#~ Xp1 = Xp1.N()
#~ 
			#~ Xmm1=Xmm1.N()
			#~ Xm0 = Xm0.N()
			#~ Xm1 = Xm1.N()


			Ypm1=Ypm1.N()
			Yp0 = Yp0.N()
			Yp1 = Yp1.N()

			Ymm1=Ymm1.N()
			Ym0 = Ym0.N()
			Ym1 = Ym1.N()

        #We set the components of the polarization vectors epsilon^(+-) in the helicity basis.
        self.Yp=[Ypm1,Yp0,Yp1]
        self.Ym=[Ymm1,Ym0,Ym1]
        
        #~ self.Xp=[Xpm1,Xp0,Xp1]
        #~ self.Xm=[Xmm1,Xm0,Xm1]
        
        if not sage_included:
			for i in range(3):
				
				a=self.Yp[i]; rea=a.real; ima=a.imag
				b=self.Ym[i]; reb=b.real; imb=b.imag
				
				if abs(rea)<1e-15: rea=0
				if abs(ima)<1e-15: ima=0
				if abs(reb)<1e-15: reb=0
				if abs(imb)<1e-15: imb=0
				
				self.Yp[i]=rea+1j*ima
				self.Ym[i]=reb+1j*imb


    def __str__(self):
        s ='Laser with phi='+str(self.phi)+', theta='+str(self.theta)
        s+=', alpha='+str(self.alpha)+', beta='+str(self.beta)
        s+=', E0='+str(self.E0)
        return s
    #~ def plot(self,dist_to_center=1,**kwds):
		#~ 
        #~ Ex=[]; Ey=[]; Ez=[]
        #~ print 333
        #~ if not sage_included:
			#~ k=[cos(self.phi)*sin(self.theta),sin(self.phi)*sin(self.theta),cos(self.theta)]
			#~ 
			#~ Nt=1000
			#~ tstep=7*pi/4/(Nt-1)
			#~ alpha=self.alpha
			#~ beta =self.beta
			#~ phi  =self.phi
			#~ theta=self.theta
			#~ omega=1
#~ #			print 222
			#~ 
			#~ for i in range(Nt):
				#~ t=i*tstep
#~ #				Ex+=[(cos(self.alpha)*cos(self.phi)*cos(self.theta) + sin(self.alpha)*sin(self.phi))*cos(self.beta + t) -
#~ #        (cos(self.phi)*cos(self.theta)*sin(self.alpha) - cos(self.alpha)*sin(self.phi))*cos(t)-dist_to_center*k[0]]
#~ #				Ey+=[(cos(self.alpha)*cos(self.theta)*sin(self.phi) - cos(self.phi)*sin(self.alpha))*cos(self.beta + t) -
#~ #        (cos(self.theta)*sin(self.alpha)*sin(self.phi) + cos(self.alpha)*cos(self.phi))*cos(t)-dist_to_center*k[1]]
#~ #				Ez+=[-cos(self.alpha)*cos(self.beta + t)*sin(self.theta) + cos(t)*sin(self.alpha)*sin(self.theta)-dist_to_center*k[2]]
			#~ 
			#~ #############################
			#~ Ex+=[(cos(2*alpha)*cos(phi)*cos(theta) - sin(2*alpha)*sin(phi))*cos(omega*t)*cos(2*beta) - (cos(phi)*cos(theta)*sin(2*alpha) + cos(2*alpha)*sin(phi))*sin(omega*t)*sin(2*beta)]
#~ 
			#~ Ey+=[(cos(2*alpha)*cos(theta)*sin(phi) + cos(phi)*sin(2*alpha))*cos(omega*t)*cos(2*beta) - (cos(theta)*sin(2*alpha)*sin(phi) - cos(2*alpha)*cos(phi))*sin(omega*t)*sin(2*beta)]
#~ 
			#~ Ez+=[-cos(omega*t)*cos(2*alpha)*cos(2*beta)*sin(theta) + sin(omega*t)*sin(2*alpha)*sin(2*beta)*sin(theta)	]
#~ 
			#~ k=[cos(phi)*sin(theta),sin(phi)*sin(theta),cos(theta)]
			#~ return [Ex,Ey,Ez]+[-dist_to_center*k[0],-dist_to_center*k[1],-dist_to_center*k[2]]
            #~ #raise NotImplementedError,"Ploting the laser's electric field is not implemented outside Sage."
        #~ else:
            #~ t=var('t')
            #~ E=((cos(self.alpha)*cos(self.phi)*cos(self.theta) + sin(self.alpha)*sin(self.phi))*cos(self.beta + t) -
            #~ (cos(self.phi)*cos(self.theta)*sin(self.alpha) - cos(self.alpha)*sin(self.phi))*cos(t),
            #~ (cos(self.alpha)*cos(self.theta)*sin(self.phi) - cos(self.phi)*sin(self.alpha))*cos(self.beta + t) -
            #~ (cos(self.theta)*sin(self.alpha)*sin(self.phi) + cos(self.alpha)*cos(self.phi))*cos(t), 
            #~ -cos(self.alpha)*cos(self.beta + t)*sin(self.theta) + cos(t)*sin(self.alpha)*sin(self.theta))
#~ 
            #~ E=vector(E)
#~ 
            #~ k=vector([cos(self.phi)*sin(self.theta),sin(self.phi)*sin(self.theta),cos(self.theta)])
            #~ E=E-dist_to_center*k
            #~ plo =point3d(E.subs(t=7*pi/4),size=15,color=self.color,**kwds)
            #~ plo+=parametric_plot3d(E,(t,0,7*pi/4),thickness=4,color=self.color,**kwds)
            #~ plo+=arrow(-dist_to_center*k,-(dist_to_center-1)*k,color=self.color,**kwds)
            #~ plo+=sphere(opacity=0)
            #~ return plo

class MotField(object):
	def __init__(self,relative_intensities,parity=1,color='blue'):
		lx   = PlaneWave(   Pi, Pi/2, 0,  parity*Pi/8,symbolical=False,color=color)
		lx_r = PlaneWave(    0, Pi/2, 0,  parity*Pi/8,symbolical=False,color=color)
		
		ly   = PlaneWave(-Pi/2, Pi/2, 0,  parity*Pi/8,symbolical=False,color=color)
		ly_r = PlaneWave( Pi/2, Pi/2, 0,  parity*Pi/8,symbolical=False,color=color)
		
		lz   = PlaneWave( 0, Pi, 0,  -parity*Pi/8,symbolical=False,color=color)
		lz_r = PlaneWave( 0,  0, 0,  -parity*Pi/8,symbolical=False,color=color)
		
		self.lx=lx; self.ly=ly; self.lz=lz
		self.lx_r=lx_r;	self.ly_r=ly_r;	self.lz_r=lz_r
		self.beams=[lx,ly,lz,lx_r,ly_r,lz_r]
		
		Yp=[0 for i in range(3)]; Ym=[0 for i in range(3)]
		for i in range(3):
			Yp[i]+=relative_intensities[0]*lx.Yp[i]
			Yp[i]+=relative_intensities[1]*ly.Yp[i]
			Yp[i]+=relative_intensities[2]*lz.Yp[i]
			Yp[i]+=relative_intensities[3]*lx_r.Yp[i]
			Yp[i]+=relative_intensities[4]*ly_r.Yp[i]
			Yp[i]+=relative_intensities[5]*lz_r.Yp[i]

			Ym[i]+=relative_intensities[0]*lx.Ym[i]
			Ym[i]+=relative_intensities[1]*ly.Ym[i]
			Ym[i]+=relative_intensities[2]*lz.Ym[i]
			Ym[i]+=relative_intensities[3]*lx_r.Ym[i]
			Ym[i]+=relative_intensities[4]*ly_r.Ym[i]
			Ym[i]+=relative_intensities[5]*lz_r.Ym[i]
		
		self.Yp=Yp
		self.Ym=Ym
		
	def plot(self,**kwds):
		plo  = self.lx.plot(  dist_to_center=3,**kwds)
		plo += self.ly.plot(  dist_to_center=3,**kwds)
		plo += self.lz.plot(  dist_to_center=3,**kwds)
		plo += self.lx_r.plot(dist_to_center=3,**kwds)
		plo += self.ly_r.plot(dist_to_center=3,**kwds)
		plo += self.lz_r.plot(dist_to_center=3,**kwds)
		return plo
	
def electric_field_amplitude_gaussian(P,sigmax,sigmay=None,Omega=1.0e6):
	'''This function returns the value of E0 (the amplitude of the electric field)
	at the center of a laser beam of power P (in Watts) and a gaussian intensity
	distribution of standard deviations sigmax, sigmay (in meters). The value of E0 is given in rescaled units
	according to the frequency scale  Omega (in Hertz) understood as absolute frequency
	(as opposed to angular frequency).'''
	e0=hbar*Omega/(e*a0) #This is the electric field scale.
	
	if sigmay==None: sigmay=sigmax
	return sqrt((c*mu0*P)/(2*Pi))/sqrt(sigmax*sigmay)/e0

def electric_field_amplitude_top(P,a,Omega=1.0e6):
	'''This function returns the value of E0 (the amplitude of the electric field)
	at the center of a laser beam of power P (in Watts) and a top-hat intensity
	distribution of radius a (in meters). The value of E0 is given in rescaled units
	according to the frequency scale  Omega (in Hertz) understood as absolute frequency
	(as opposed to angular frequency).'''
	e0=hbar*Omega/(e*a0) #This is the electric field scale.
	
	return sqrt((c*mu0*P)/(Pi*a**2))/e0

def electric_field_amplitude_intensity(s0,Omega=1.0e6):
	'''This function returns the value of E0 (the amplitude of the electric field)
	at a given saturation parameter s0=I/I0, where I0=2.50399 mW/cm^2 is the
	saturation intensity of the D2 line of Rubidium for linearly polarized light.'''
	e0=hbar*Omega/(e*a0) #This is the electric field scale.
	
	I0=2.50399 #mW/cm^2
	I0=1.66889451102868 #mW/cm^2
	
	I0=I0/1000*(100**2) #W/m^2
	r_ciclic=4.226983616875483 #a0
	gamma_D2=2*Pi*6.065e6/Omega # The decay frequency of the D2 line.
	E0_sat=gamma_D2/r_ciclic/sqrt(2.0)

	E0_sat=E0_sat*e0
	I0=E0_sat**2/2/c/mu0

    #return sqrt(c*mu0*s0*I0/2)/e0
	#return sqrt(c*mu0*s0*I0)/e0
	return sqrt(2*c*mu0*s0*I0)/e0
	
