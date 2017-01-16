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

if not sage_included:
	from sympy.core.numbers import Rational as Integer
	#from fractions import Fraction as Integer
	from math import sqrt,pi
	from sympy.physics.wigner import wigner_3j,wigner_6j

class State(object):
    def __init__(self,element,isotope,n,l,j,f=None,m=None):
        '''This class is for atomic states, it allows to specify any atomic state of the form 
|n,l,j>, |n,l,j,f>, |n,l,j,f,m> which are refered to respectively as states with fine 
hyperfine and magnetic detail.'''

        self.element=element			
    
        #We go from chemical notation to quantum numbers.
        if str(l)=='S': l=0
        if str(l)=='P': l=1
        if str(l)=='D': l=2
        if str(l)=='F': l=3

        #We check the value of l.
        lperm=range(0,n)
        if l not in lperm: raise ValueError, 'l = '+str(l)+' is not allowed.'
    
        #We check the value of j.
        jmin=abs(l-Integer(1)/Integer(2)); nj=int(2*min(l,Integer(1)/Integer(2))+1)
        jperm=[jmin]+[jmin + ii for ii in range(1,nj)]
        if j not in jperm: raise ValueError, 'j = '+str(j)+' is not allowed.'

        self.isotope=isotope
        self.n=n
        self.l=l
        self.j=j
        self.f=f
        self.m=m
        self.quantum_numbers=[isotope,n,l,j]			

        #~ if isotope==85:
			#~ i=Integer(5)/Integer(2)
        #~ elif isotope==87:
			#~ i=Integer(3)/Integer(2)
        #~ else:
			#~ raise ValueError,str(isotope)+' is not a valid isotope.'

        
        #We calculate the energy of the state.
        ################################################################
        #All tables are given in (cm^-1).
        c=299792458
        if element=="Rb":
            if isotope==85:
                i=5/Integer(2)
                #        N, L,     K       , E (cm^-1),      A (cm^-1)      B (cm^-1)
                nivfin=[[5, 0, 1/Integer(2), 0.0000000, 0.033753721   , 0.00000000],
                        [5, 1, 1/Integer(2), 12578.950, 0.004026000   , 0.00000000],
                        [5, 1, 3/Integer(2), 12816.545, 0.000835200   , 0.00086760],
                        
                        [4, 2, 5/Integer(2), 19355.203, 0.000168000   , 0.00000000],
                        [4, 2, 3/Integer(2), 19355.649, 0.000244000   , 0.00000000],
                        
                        [6, 0, 1/Integer(2), 20132.460, 239.18e6/(c*100), 0.00000000],
                        
                        [6, 1, 1/Integer(2), 23715.081, 0.001305000   , 0.00000000],
                        [6, 1, 3/Integer(2), 23792.591, 0.000272300   , 0.00027320],
                        
                        [5, 2, 3/Integer(2), 25700.536, 0.000140834   , 0.00006373],
                        [5, 2, 5/Integer(2), 25703.498,-0.000073090   , 0.00008940],
                        
                        [7, 0, 1/Integer(2), 26311.437, 0.003159000   , 0.00000000],
                        [7, 1, 1/Integer(2), 27835.020, 0.000590000   , 0.00000000],
                        [7, 1, 3/Integer(2), 27870.110, 0.000123800   , 0.00012300]]

            elif isotope==87:
                i=3/Integer(2)
                #        N, L,     K       , E (cm^-1),      A (cm^-1)      B (cm^-1)
                nivfin=[[5, 0, 1/Integer(2), 0.0000000, 0.113990236053642, 0.00000000],
                        [5, 1, 1/Integer(2), 12578.950, 0.013650000000000, 0.00000000],
                        [5, 1, 3/Integer(2), 12816.545, 0.002825900000000, 0.00041684],
                        
                        [4, 2, 5/Integer(2), 19355.203,-16.801e6/(c*100), 3.645e6/(c*100)],
                        [4, 2, 3/Integer(2), 19355.649, 24.750e6/(c*100), 2.190e6/(c*100)],

                        [6, 0, 1/Integer(2), 20132.460, 807.66e6/(c*100), 0.00000000],
                        
                        [6, 1, 1/Integer(2), 23715.081, 0.004421700000000, 0.00000000],
                        [6, 1, 3/Integer(2), 23792.591, 0.000924000000000, 0.00013200],
                        
                        [5, 2, 3/Integer(2), 25700.536, 0.000481340000000, 0.00003109],
                        [5, 2, 5/Integer(2), 25703.498, -0.00024886000000, 0.00004241],
                        
                        [7, 0, 1/Integer(2), 26311.437, 0.010664000000000, 0.00000000],
                        [7, 1, 1/Integer(2), 27835.020, 0.001999000000000, 0.00000000],
                        [7, 1, 3/Integer(2), 27870.110, 0.000419300000000, 0.00005700]]

            else:
                s="The isotope "+str(isotope)+str(element)+" is not in the database."
                raise ValueError,s

        elif element=="Cs":
            if isotope==133:
                i=7/Integer(2)
                #        N, L,     K       , E (cm^-1),                     A (cm^-1)      B (cm^-1)
                nivfin=[[6, 0, 1/Integer(2), 0.0000000               , 2.2981579425e9/(c*100), 0.0],# A is exact.
                        [6, 1, 1/Integer(2), 335.116048807e12/(c*100),     291.9201e6/(c*100), 0.0],
                        [6, 1, 3/Integer(2), 351.72571850e12 /(c*100),     50.28827e6/(c*100),-493.4e3/(c*100)]]# C: 560.0e0/(c*100)
            else:
                s="The isotope "+str(isotope)+str(element)+" is not in the database."
                raise ValueError,s
        else:
            s="The element "+str(element)+" is not in the database."
            raise ValueError,s
            
		
		#We rewrite the table in Hz
        nivfin=[ nivfin[ii][:3] + [nivfin[ii][3]*c*100 ] +[nivfin[ii][4]*c*100 ] +[nivfin[ii][5]*c*100 ] for ii in range(len(nivfin)) ]

        #We check whether the quantum numbers given are in the database.
        #We find the energy of the state up to fine structure.
        in_database=False
        for ii in range(len(nivfin)):
            #print [n,l,j],nivfin[ii][:3]
            if [n,l,j]==nivfin[ii][:3]:
                nufin=nivfin[ii][3]; A=nivfin[ii][4]; B=nivfin[ii][5]
                in_database=True
                break
        if not in_database:
            s ="The values of n,l,k: "+str(n)+", "+str(l)+", "+str(j)
            s+=" are not in the database for "+str(isotope)+str(element)+"."
            raise ValueError,s


        #We check the value of f.
        fmin=int(abs(j-i)); nf=int(2*min(j,i)+1)
        fperm=[fmin]+[fmin + ii for ii in range(1,nf)]        
        if f != None:
            if f in fperm:
				self.quantum_numbers.append(f)
				mperm=[ii for ii in range(-f,f+1)]
				if m != None:
					if m in mperm:
						self.quantum_numbers.append(m)
					else:
						raise ValueError, 'm = '+str(m)+' is not allowed.'
				else:
					self.mperm=mperm
            else:
            	raise ValueError, 'f = '+str(f)+' is not allowed.'
        else:
            self.fperm=fperm



        #print fperm
        #print nufin,A,B
        #Establecemos la corrección hiperfina de la energía o la estructura hiperfina.
        if f != None:
            K=f*(f+1)-i*(i+1)-j*(j+1)
            deltanu_hfin=A*K/Integer(2)
            if j != Integer(1)/Integer(2):
				deltanu_hfin= deltanu_hfin + B*( 3*K*(K+1)/Integer(2) -2*i*(i+1)*j*(j+1) )/(4*i*(2*i-1)*j*(2*j-1))
            deltanu_hfin=float(deltanu_hfin)
            self.nu= nufin + deltanu_hfin
            
        else:
            estructura_hiperfina=[]
            for f in fperm:
                K=f*(f+1)-i*(i+1)-j*(j+1)
                deltanu_hfin=A*K/Integer(2)
                if j != Integer(1)/Integer(2):
					deltanu_hfin= deltanu_hfin + B*( 3*K*(K+1)/Integer(2) -2*i*(i+1)*j*(j+1) )/(4*i*(2*i-1)*j*(2*j-1))
                estructura_hiperfina.append( (deltanu_hfin,f) )
            estructura_hiperfina.sort()
            
            self.estructura_hiperfina=estructura_hiperfina
            self.nu=nufin
        

    def __repr__(self):
        if self.l==0:
            l='S'
        elif self.l==1:
            l='P'
        elif self.l==2:
            l='D'
        elif self.l==3:
            l='F'
        else:
            l=str(self.l)
            
        if self.f==None:
            s= str(self.isotope)+self.element+' '+str(self.n)+l+'_'+str(self.j)
        else:
            s= str(self.isotope)+self.element+' '+str(self.n)+l+'_'+str(self.j)+'^'+str(self.f)
        
        if self.m !=None:
			s += ','+str(self.m)
        return s

    def _latex_(self):
        if self.l==0:
            l='S'
        elif self.l==1:
            l='P'
        elif self.l==2:
            l='D'
        elif self.l==3:
            l='F'
        else:
            l=str(self.l)
        
        if self.f==None:
            s= '^{'+str(self.isotope)+'}\\mathrm{'+self.element+'}\\ '+str(self.n)+l+'_{'+str(self.j)+'}'
        else:
            s= '^{'+str(self.isotope)+'}\\mathrm{'+self.element+'}\\ '+str(self.n)+l+'_{'+str(self.j)+'}^{'+str(self.f)+'}'
        
        if self.m != None:
			s= s[:-1] + ','+str(self.m)+'}'
        
        return s
    
    def __eq__(self,other):
		return self.quantum_numbers==other.quantum_numbers

class Transition(object):
	def __init__(self,e1,e2):
		'''This class describes a transition between different states of the same level of detail.'''
		if e1.isotope != e2.isotope:
			raise ValueError,'Transitions between diferent isotopes are forbidden.'
		
		self.e1=e1; self.e2=e2
		if e1.f==None or e2.f==None:
			l1=e1.l; l2=e2.l
			j1=e1.j; j2=e2.j
			allowed=True
			if not (abs(j2-j1) in [-1,0,1]): allowed=False
			if j1==0 and j2==0: allowed=False
			if not (abs(l2-l1) in [-1 ,1]): allowed=False
			if l1==0 and l2==0: allowed=False
			
			#In dipolar electric transitions states must have oposite parity.
			#In particular, transitions between identical fine states are forbidden.
			if e1==e2: allowed=False
			self.allowed=allowed
			#raise ValueError, 'Las transiciones deben ser entre estados hiperfinos (con una f específica).'
		else:			
			if self.e1.f - self.e2.f in [-1,0,1]:
				if self.e1 != self.e2:
					if self.e1.nu > self.e2.nu:
						if self.e1.m==None and self.e2.m==None:
							self.allowed=True
						elif self.e1.m -self.e2.m in [-1,0,1]:
							self.allowed=True
						else:
							self.allowed=False
					else:
						self.allowed=False
				else:
					self.allowed=False
			else:
				self.allowed=False
		
		self.nu=self.e1.nu-self.e2.nu
		if self.nu ==0:
			self.wavelength=float('inf')
		else:
			self.wavelength= 299792458.0 / self.nu
		self.omega=self.nu*2*pi
		
		#We find the Einstein A and B coefficients up to the fine structure
		#according to literature (in Hz).
		n1=e1.n; l1=e1.l; j1=e1.j
		n2=e2.n; l2=e2.l; j2=e2.j
		ord1=[n1,l1,j1,n2,l2,j2]
		ord2=[n2,l2,j2,n1,l1,j1]
		
		####################################################################
		# This is the database of Einstein A coefficients.
		# The quantum numbers of the lower states followed by those of excited states
		# followed by Einstein A coefficients in Hz (angular frequency).
		pairs= [[5, 1, Integer(1)/Integer(2),   5, 0, Integer(1)/Integer(2), 3.61031827750539e7],
				[5, 1, Integer(3)/Integer(2),   5, 0, Integer(1)/Integer(2), 3.81075188880442e7],
				[4, 2, Integer(5)/Integer(2),   5, 1, Integer(3)/Integer(2), 1.06814150222053e7],
				[4, 2, Integer(3)/Integer(2),   5, 1, Integer(1)/Integer(2), 9.42477796076938e6],
				[4, 2, Integer(3)/Integer(2),   5, 1, Integer(3)/Integer(2), 1.88495559215388e6],
				[6, 0, Integer(1)/Integer(2),   5, 1, Integer(1)/Integer(2), 1.09704415463356e7],#The branching ratios here are not well referenced.
				[6, 0, Integer(1)/Integer(2),   5, 1, Integer(3)/Integer(2), 1.09704415463356e7],#The branching ratios here are not well referenced.
				[6, 1, Integer(3)/Integer(2),   5, 0, Integer(1)/Integer(2), 1.87867240684670e6],
				[6, 1, Integer(3)/Integer(2),   4, 2, Integer(5)/Integer(2), 179699.099785336],
				[6, 1, Integer(3)/Integer(2),   4, 2, Integer(3)/Integer(2), 1.61729189806803e6],
				[6, 1, Integer(3)/Integer(2),   6, 0, Integer(1)/Integer(2), 4.49247749463340e6],
				[5, 2, Integer(5)/Integer(2),   5, 1, Integer(3)/Integer(2), 3.10264947105589e6],
				[5, 2, Integer(5)/Integer(2),   6, 1, Integer(3)/Integer(2), 1.09012008442504e6]]
		
		self.einsteinA=0.0
		#if ord1 == ord2: self.einsteinA=0.0
		if self.allowed: self.einstein=None
		for pair in pairs:
			if pair[:-1]==ord1 or pair[:-1]==ord2:
				self.einsteinA=pair[-1]
		
		
		
		#~ #The numbers between 5S1/2 -> 5P3/2
		#~ 
		#~ 
		#~ pair1=[5,0,Integer(1)/Integer(2),   5,1,Integer(3)/Integer(2)]
		#~ #print ord1,ord2,pair1
		#~ pair1=ord1==pair1 or ord2==pair1
		#~ #print pair1
		#~ 
		#~ A1=6.0666e6
		#~ ####################################################################
		#~ 
		#~ if pair1:
			#~ self.einsteinA=A1
		#~ elif ord1==ord2:
			#~ self.einsteinA=0.0
		#~ else:
			#~ self.einsteinA=None
		
	
	def __repr__(self):
		if self.allowed==True:
			return self.e1.__repr__()+' -----> '+self.e2.__repr__()
		elif self.allowed==False:
			return self.e1.__repr__()+' --/--> '+self.e2.__repr__()
		else:
			return self.e1.__repr__()+' --?--> '+self.e2.__repr__()

	def _latex_(self):
		if self.allowed==True:
			return self.e1._latex_()+'\\ \\rightarrow \\ '+self.e2._latex_()
		elif self.allowed==False:
			return self.e1._latex_()+'\\ \\nrightarrow \\ '+self.e2._latex_()
		else:
			return self.e1._latex_()+'\\ \\rightarrow^? \\ '+self.e2._latex_()

        
		return self.e1._latex_()+'\\ \\nleftrightarrow \\ '+self.e2._latex_()

	def __eq__(self,other):
		return self.e1==other.e1 and self.e2==other.e2

def split_fine_to_hyperfine(state):

	if type(state)==list:
		mag=[]
		for s in state:
			mag+=split_fine_to_hyperfine(s)
		return mag
	
	if len(state.quantum_numbers) != 4:
		s=str(state)+' no es un estado fino.'
		raise ValueError,s
	
	return [State(state.element,state.isotope,state.n,state.l,state.j,f) for f in state.fperm]
	
def split_hyperfine_to_magnetic(state):
	if type(state)==list:
		mag=[]
		for s in state:
			mag+=split_hyperfine_to_magnetic(s)
		return mag


	if len(state.quantum_numbers) != 5:
		raise ValueError,str(state)+' no es un estado hiperfino.'
	
	#print state.j,type(state.j)
	
	return [State(state.element,state.isotope,state.n,state.l,state.j,state.f,m) for m in state.mperm]

def split_fine_to_magnetic(state):
	if type(state)==list:
		mag=[]
		for s in state:
			mag+=split_fine_to_magnetic(s)
		return mag
	
	if len(state.quantum_numbers) != 4:
		raise ValueError,str(state)+' no es un estado fino.'
	
	hip=split_fine_to_hyperfine(state)
	
	mag=[]
	for ll in [split_hyperfine_to_magnetic(h) for h in hip]:
		mag+=ll
	return mag

def order_by_energy(states):
    iso=states[0].isotope
    for ss in states:
        if ss.isotope != iso:
            raise ValueError,'We have a wrong isotope in this list: '+str(ss)
    aux=[(ee.nu,ee) for ee in states]
    aux.sort()
    
    return [i[1] for i in aux ]

def make_list_of_states(states,structure=None,verbose=1):
	if structure==None:
		return order_by_energy(states)
	elif structure=='hyperfine':
		l1=order_by_energy(states)
		l2=split_fine_to_hyperfine(l1)
		l3=order_by_energy(l2)
		if l2 != l3:
			if verbose>0: print 'Warning: the ordering of the hyperfine states has been changed to ensure they are ordered by ascending energy.'
		return l3
	elif structure=='magnetic':
		l1=order_by_energy(states)
		l2=split_fine_to_hyperfine(l1)
		l3=order_by_energy(l2)
		if l2 != l3:
			if verbose>0: print 'Warning: the ordering of the hyperfine states has been changed to ensure they are ordered by ascending energy.'
		l4=split_hyperfine_to_magnetic(l3)
		#print l4[0].j,type(l4[0].j)
		return l4

def calculate_omega_matrix(states,Omega=1):
	"""This function recieves a list of states and returns the corresponding
	omega_ij matrix, rescaled to units of Omega. These are understood to be
	absolute frequencies (as opposed angular frequencies)."""
	
	Pi = 3.14159265358979
	N=len(states)
	omega=[[ 2*Pi*(states[i].nu-states[j].nu)/Omega for j in range(N)] for i in range(N)]

	return omega

def get_einstein_A_matrix(fine_states,Omega=1):
	
	einsteinA=[[0.0 for jj in range(len(fine_states))] for ii in range(len(fine_states))]
	for ii in range(len(fine_states)):
		i=fine_states[ii]
		for jj in range(ii):
			j=fine_states[jj]
			t=Transition(i,j)
			Aij=t.einsteinA

			if Aij==None:
				s ='The Einstein coeficcient A_ij between '+str(i)+' and '+str(j)+' is unknown, '
				s+='please add it to the database in the code of the class "Transition" to proceed.'
				raise NotImplementedError,s
			einsteinA[ii][jj]=Aij/Omega
			einsteinA[jj][ii]=-Aij/Omega
	
	return einsteinA

def calculate_gamma_matrix(magnetic_states, Omega=1):
	#index_list_fine, index_list_hyperfine = calculate_boundaries(fine_states, full_magnetic_states)
	Ne=len(magnetic_states)

	if magnetic_states[0].isotope==85:
		II=Integer(5)/Integer(2)
	elif magnetic_states[0].isotope==87:
		II=Integer(3)/Integer(2)
	#II=float(II)
	
	gamma=[ [0.0 for j in range(Ne)] for i in range(Ne)]
	for i in range(Ne):
		for j in range(i):
			ei=magnetic_states[i]
			ej=magnetic_states[j]
			einsteinAij=Transition(ei,ej).einsteinA
			
			if einsteinAij != 0:
				ji=ei.j; jj=ej.j
				fi=ei.f; fj=ej.f
				mi=ei.m; mj=ej.m
				
				gammaij =(2.0*ji+1)
				gammaij*=(2.0*fi+1)
				gammaij*=(2.0*fj+1)
				
				gammaij*=float(wigner_6j(ji,fi,II,fj,jj,1)**2)
				gammaij*=sum([float(wigner_3j(fj,1,fi,-mj,q,mi)**2) for q in [-1,0,1]])
				
				gammaij*=einsteinAij/Omega
				gammaij=float(gammaij)
				
				gamma[i][j]= gammaij
				gamma[j][i]=-gammaij
	
	return gamma

def calculate_reduced_matrix_elements(fine_states):
	'''This function calculates the reduced matrix elments <N,L,J||T^1(r)||N',L',J'> given a list of fine states.'''
	#We calculate the reduced matrix elements starting from the list of fine_states
	#In SI units
	epsilon0 = 8.854187817e-12  #The permitivity of vacuum
	c        = 299792458.0      #The speed of light
	h        = 1.054571726e-34  #Plank's reduced constant
	e        = 1.602176565e-19  #The charge of the electron
	me       = 9.10938291e-31   #The electron rest mass
	Pi       = 3.14159265358979 #Everyone's favourite constant

	#The factor composed of physical quantities.
	factor=sqrt(3*c**3*me**2*e**2/(16*Pi*epsilon0*h**3))	
	#We read the database to obtain the Einstein A coefficients in Hz.
	einsteinA=get_einstein_A_matrix(fine_states)
	#We read the database to obtain the transition frequencies in Hz.
	omega_fine=calculate_omega_matrix(fine_states)
	
	reduced_matrix_elements=[[0.0 for jj in range(len(fine_states))] for ii in range(len(fine_states))]
	
	for ii in range(len(fine_states)):
		i=fine_states[ii]
		for jj in range(ii):
			j=fine_states[jj]
			t=Transition(i,j)
			einsteinAij=einsteinA[ii][jj]
			omega0=omega_fine[ii][jj]
			
			#The formula is valid only for i =/= j so that omega0 =/= 0.
			#Because fine states are asumed to be ordered by their energies we can asume that
			# i decays in j.
			Ji=i.j; Jj=j.j

			rij=(2.0*Ji+1)/sqrt(2.0*Jj+1)*sqrt(einsteinAij/omega0**3)
			rij=factor*rij
		
			reduced_matrix_elements[ii][jj]=rij
			#We add the matrix elements on the other side of the diagonal.
			reduced_matrix_elements[jj][ii]=rij*(-1)**(Ji-Jj) 
	
	return reduced_matrix_elements

def calculate_r_matrices(fine_states, reduced_matrix_elements):
	
	magnetic_states=make_list_of_states(fine_states,'magnetic',verbose=0)
	index_list_fine, index_list_hyperfine = calculate_boundaries(fine_states, magnetic_states)
	
	Ne=len(magnetic_states)
		
	r=[[[0.0 for j in range(Ne)] for i in range(Ne)] for p in range(3)]
	
	if fine_states[0].isotope==85:
		II=Integer(5)/Integer(2)
	elif fine_states[0].isotope==87:
		II=Integer(3)/Integer(2)

	for p in [-1,0,1]:
		for i in range(Ne):
			ei=magnetic_states[i]
			ii=fine_index(i, index_list_fine)
			
			for j in range(Ne):
				ej=magnetic_states[j]
				jj=fine_index(j, index_list_fine)
				
				reduced_matrix_elementij=reduced_matrix_elements[ii][jj]
				if reduced_matrix_elementij != 0:
					
					ji=ei.j; jj=ej.j
					fi=ei.f; fj=ej.f
					mi=ei.m; mj=ej.m
										
					rpij =(-1)**(fi-mi)
					#print fi,1,jj,-mi,p, mj
					rpij*=wigner_3j(fi,1,fj,-mi,p, mj)
					
					rpij*=(-1)**(fj+ji+1+II)
					rpij*=sqrt(2*fj+1)
					rpij*=sqrt(2*fi +1)
					rpij*=wigner_6j(ji,jj,1,fj,fi,II)
					
					rpij*=reduced_matrix_elementij
					
					r[p+1][i][j]=float(rpij)
	return r

def calculate_matrices(states,Omega=1):
	'''This function calculates the matrices omega_ij, gamma_ij and r_pij given a list
	of atomic states. The states can be arbitrarily in their fine, hyperfine or magnetic detail.'''
	
	#We check that all states belong to the same isotope
	iso=states[0].isotope
	for state in states[1:]:
		if state.isotope != iso:
			raise ValueError,'All states must belong to the same isotope.'
	
	#We find the fine states involved in the problem.
	fine_states=find_fine_states(states)

	#We find the full magnetic states. The matrices will be first calculated for the complete
	#problem and later reduced to include only the states of interest.
	full_magnetic_states=make_list_of_states(fine_states,'magnetic',verbose=0)

	#We calculate the indices corresponding to each sub matrix of fine and hyperfine levels.
	#index_list_fine, index_list_hyperfine = calculate_boundaries(fine_states, full_magnetic_states)
	
	#We calculate the frequency differences between states.
	omega_full=calculate_omega_matrix(full_magnetic_states,Omega)
	#We calculate the matrix gamma

	gamma_full=calculate_gamma_matrix(full_magnetic_states, Omega)
	#We calculate the reduced matrix elements
	reduced_matrix_elements=calculate_reduced_matrix_elements(fine_states)
	
	#We calculate the matrices r_-1, r_0, r_1
	r_full=calculate_r_matrices(fine_states,reduced_matrix_elements)
		
	#Reduction to be implemented
	omega=omega_full
	r=r_full
	gamma=gamma_full
	return omega,gamma,r

def calculate_boundaries(fine_states,full_magnetic_states):
	"""This function calculates the boundary indices of each fine state and each hyperfine state
	within a list of magnetic states. The output is a list of tuples (a,b) with a the starting index
	of a state and b it's ending index."""
	N_magnetic=len(full_magnetic_states)
	
	#We calculate the boundaries of the various detail levels
	#First we will make a list of indices of that will tell where each fine level begins and ends.
	fq=full_magnetic_states[0].quantum_numbers[:4]
	index_list_fine=[]; start_fine=0
	#And another list of indices that will tell where each hyperfine level begins and ends.
	hq=full_magnetic_states[0].quantum_numbers[:5]
	index_list_hyperfine=[]; start_hyperfine=0
	
	for i in range(N_magnetic):
		magnetic=full_magnetic_states[i]
		if magnetic.quantum_numbers[:4]!=fq:
			index_list_fine+=[(start_fine,i)]
			start_fine=i
			fq=magnetic.quantum_numbers[:4]

		if magnetic.quantum_numbers[:5]!=hq:
			index_list_hyperfine+=[(start_hyperfine,i)]
			start_hyperfine=i
			hq=magnetic.quantum_numbers[:5]
			

		if i==N_magnetic-1:
			index_list_fine+=[(start_fine,i+1)]
			index_list_hyperfine+=[(start_hyperfine,i+1)]
	return index_list_fine,index_list_hyperfine

def fine_index(magnetic_index,index_list_fine):
	""""""
	N_fine=len(index_list_fine)
	for i in range(N_fine):
		if index_list_fine[i][0] <= magnetic_index < index_list_fine[i][1]:
			return i

def quaver(isotope,p,  J,F,M,  Jp,Fp,Mp, numeric=False,verbose=False):
	if isotope==85:
		II=Integer(5)/Integer(2)
	elif isotope==87:
		II=Integer(3)/Integer(2)
	
	qu =(-1)**(F-M+J+II+Fp+1)
	qu*=sqrt((2*F+1)*(2*Fp+1))
	if not sage_included:
		II=float(II)
		J=float(J)
		Jp=float(Jp)
	qu*=wigner_3j(F,1,Fp,  -M,p,Mp)*wigner_6j(II,J,F,  1,Fp,Jp)
	
	if numeric: return float(qu)
	else: return qu

#~ def calculate_gamma_matrix(full_magnetic_states, reduced_matrix_elements, index_list_fine, Omega=1):
	#~ epsilon0 = 8.854187817e-12  #The permitivity of vacuum
	#~ c        = 299792458.0      #The speed of light
	#~ h        = 1.054571726e-34  #Plank's reduced constant
	#~ e        = 1.602176565e-19  #The charge of the electron
	#~ Pi       = 3.14159265358979
	#~ N_magnetic=len(full_magnetic_states)
	#~ gamma=[]
	#~ for i in range(N_magnetic):
		#~ row=[]
		#~ ISO,N,L,J,F,M=full_magnetic_states[i].quantum_numbers
		#~ ii=fine_index(i, index_list_fine)
		#~ for j in range(N_magnetic):
			#~ ISOp,Np,Lp,Jp,Fp,Mp=full_magnetic_states[j].quantum_numbers
			#~ jj=fine_index(j, index_list_fine)
			#~ 
			#~ t=Transition(full_magnetic_states[i],full_magnetic_states[j])
			#~ omegaij=2*Pi*t.nu
			#~ red=reduced_matrix_elements[ii][jj]
			#~ s=sum([quaver(ISO,p, J, F, M, Jp, Fp, Mp,numeric=True)**2 for p in [-1,0,1]])
			#~ row+=[4*e**2*omegaij**3*red**2*s/(3*epsilon0*c**3*h*Omega)]
		#~ gamma+=[row]
	#~ return gamma

def find_fine_states(magnetic_states):
	fine_states=[]
	for state in magnetic_states:
		fq=state.quantum_numbers[:4]
		fine_state=State(state.element,fq[0],fq[1],fq[2],fq[3])
		if fine_state not in fine_states: fine_states+=[fine_state]
	return fine_states

def exclude_states(omega,gamma,r,Lij,states,excluded_states):
	"""This function takes the matrices and excludes the states listed in excluded_states."""
	Ne=len(omega)
	excluded_indices=[i for i in range(Ne) if states[i] in excluded_states]

	omega_new=[]; gamma_new=[]; r_new=[[],[],[]]; Lij_new=[]
	for i in range(Ne):
		row_om=[]; row_ga=[]; row_L=[]
		for j in range(Ne):
			if j not in excluded_indices:
				row_om+=[omega[i][j]]
				row_ga+=[gamma[i][j]]
				row_L +=[Lij[i][j]]
		if i not in excluded_indices:
			omega_new+=[row_om]
			gamma_new+=[row_ga]
			Lij_new+=[row_L]
	
	for p in range(3):
		for i in range(Ne):
			row_r=[]
			for j in range(Ne):
				if j not in excluded_indices:
					row_r+=[r[p][i][j]]
			if i not in excluded_indices:
				r_new[p]+=[row_r]
	
	states_new=[states[i] for i in range(Ne) if i not in excluded_indices]
	
	return omega_new, gamma_new, r_new, Lij_new, states_new

def reduce_magnetic_to_hyperfine(omega,gamma,r,Lij,magnetic_states,hyperfine_states,isotropic_r=False):
	#We find the fine states involved in the problem.
	fine_states=find_fine_states(magnetic_states)

	#We calculate the indices corresponding to each sub matrix of fine and hyperfine levels.
	index_list_fine,index_list_hyperfine=calculate_boundaries(fine_states,magnetic_states)
	
	#We determine which indices will be reduced
	index_list_to_reduce0=[]
	for a,b in index_list_hyperfine:
		qn = magnetic_states[a].quantum_numbers[:5]
		for hs in hyperfine_states:
			if hs.quantum_numbers == qn:
				index_list_to_reduce0+=[(a,b)]
	
	Ne_magnetic=len(magnetic_states)
	index_list_to_reduce=[]
	for i in range(Ne_magnetic):
		band=True
		for a,b in index_list_to_reduce0:
			if a <= i < b:
				if (a,b) not in index_list_to_reduce:
					index_list_to_reduce+=[(a,b)]
				band=False
				break
		if band:
			index_list_to_reduce+=[(i,i+1)]
	
	#We calculate the x,y,z matrices in order to reduce them.
	x=[[   (r[0][i][j]-r[2][i][j])/sqrt(2.0) for j in range(Ne_magnetic)] for i in range(Ne_magnetic)]
	y=[[1j*(r[0][i][j]+r[2][i][j])/sqrt(2.0) for j in range(Ne_magnetic)] for i in range(Ne_magnetic)]
	z=[[r[1][i][j]                           for j in range(Ne_magnetic)] for i in range(Ne_magnetic)]

	#We make the reduction
	Ne_reduced = len(index_list_to_reduce)
	omega_red=[]; gamma_red=[]; xyz=[[],[],[]]; Lij_red=[]; states_red=[]
	for ii in range(Ne_reduced):
		aii,bii=index_list_to_reduce[ii]
		omega_row=[]; gamma_row=[]; x_row,y_row,z_row=[],[],[]; Lij_row=[]
		for jj in range(Ne_reduced):
			ajj,bjj=index_list_to_reduce[jj]
			
			giijj=0; xiijj,yiijj,ziijj=0,0,0; Liijj=[]
			for i in range(aii,bii):
				for j in range(ajj,bjj):
					#We can simply take whichever submatrix element for omega.
					oiijj=omega[i][j]
					#We sum the elements of gamma.
					giijj+=gamma[i][j]
					#We make the quadrature sum the elements of x,y,z.
					xiijj+=x[i][j].conjugate()*x[i][j]
					yiijj+=y[i][j].conjugate()*y[i][j]
					ziijj+=z[i][j].conjugate()*z[i][j]
					#For Lij we include whichever l exists in the submatrix in the reduced Lij.
					for l in Lij[i][j]:
						if l not in Liijj: Liijj+=[l]

			omega_row+=[oiijj]
			gamma_row+=[giijj]
			x_row+=[sqrt(xiijj.real)]
			y_row+=[sqrt(yiijj.real)]
			z_row+=[sqrt(ziijj.real)]
			Lij_row+=[Liijj]
		#For the states we check wether there is a reduction, and add the apropiate states.
		element_n=magnetic_states[aii].element
        qn=magnetic_states[aii].quantum_numbers[:5]
        hs=State(element_n,qn[0], qn[1], qn[2], qn[3], qn[4])
        if bii-aii==1:
			if split_hyperfine_to_magnetic([hs])==[magnetic_states[aii]]:
				if hs in hyperfine_states:
					states_red+=[hs]
				else:
					states_red+=[magnetic_states[aii]]
			else:
				states_red+=[magnetic_states[aii]]
        else:
			qn=magnetic_states[aii].quantum_numbers[:5]
			element_n=magnetic_states[aii].element
			states_red+=[State(element_n,qn[0], qn[1], qn[2], qn[3], qn[4])]
        omega_red+=[omega_row]
        gamma_red+=[gamma_row]
        xyz[0]+=[x_row]; xyz[1]+=[y_row]; xyz[2]+=[z_row]
        Lij_red+=[Lij_row]

	#We calculate the reduced r from xyz
	if isotropic_r:
		#We can impose the isotropy of r
		xyz[0]=xyz[2][:]; xyz[1]=xyz[2][:]
	rm1_red=[[ ( xyz[0][i][j]-1j*xyz[1][i][j])/sqrt(2.0) for j in range(Ne_reduced)] for i in range(Ne_reduced)]
	r0_red = xyz[2]
	rp1_red=[[ (-xyz[0][i][j]-1j*xyz[1][i][j])/sqrt(2.0) for j in range(Ne_reduced)] for i in range(Ne_reduced)]
	r_red=[rm1_red, r0_red, rp1_red]

	return omega_red,gamma_red,r_red,Lij_red,states_red

def calculate_reduced_matrix_elements_0(fine_states):
	'''This function calculates the reduced matrix elments <N,L,J||T^1(r)||N',L',J'> given a list of fine states.'''
	#We calculate the reduced matrix elements starting from the list of fine_states
	#In SI units
	epsilon0 = 8.854187817e-12  #The permitivity of vacuum
	c        = 299792458.0      #The speed of light
	h        = 1.054571726e-34  #Plank's reduced constant
	e        = 1.602176565e-19  #The charge of the electron
	me       = 9.10938291e-31   #The electron rest mass
	Pi       = 3.14159265358979 #Everyone's favourite constant

	#The factor composed of physical quantities.
	factor=sqrt(3*c**3*me**2*e**2/(16*Pi*epsilon0*h**3))	
	#We read the database to obtain the Einstein A coefficients in Hz.
	einsteinA=get_einstein_A_matrix(fine_states)
	#We read the database to obtain the transition frequencies in Hz.
	omega_fine=calculate_omega_matrix(fine_states)
	
	reduced_matrix_elements=[[0.0 for jj in range(len(fine_states))] for ii in range(len(fine_states))]
	
	for ii in range(len(fine_states)):
		i=fine_states[ii]
		for jj in range(ii):
			j=fine_states[jj]
			t=Transition(i,j)
			einsteinAij=einsteinA[ii][jj]
			omega0=omega_fine[ii][jj]
			
			#The formula is valid only for i =/= j so that omega0 =/= 0.
			#Because fine states are asumed to be ordered by their energies we can asume that
			# i decays in j.
			Ji=i.j; Jj=j.j

			rij=sqrt((2.0*Ji+1)/(2*Jj+1))*sqrt(einsteinAij/omega0**3)
			rij=factor*rij
		
			reduced_matrix_elements[ii][jj]=rij
			#We add the matrix elements on the other side of the diagonal.
			reduced_matrix_elements[jj][ii]=rij*(-1)**(Jj-Ji) 
	
	return reduced_matrix_elements

def calculate_reduced_matrix_elements_steck(fine_states):
	'''This function calculates the reduced matrix elments <N,L,J||T^1(r)||N',L',J'> given a list of fine states.'''
	#We calculate the reduced matrix elements starting from the list of fine_states
	#In SI units
	epsilon0 = 8.854187817e-12  #The permitivity of vacuum
	c        = 299792458.0      #The speed of light
	h        = 1.054571726e-34  #Plank's reduced constant
	e        = 1.602176565e-19  #The charge of the electron
	me       = 9.10938291e-31   #The electron rest mass
	Pi       = 3.14159265358979 #Everyone's favourite constant

	#The factor composed of physical quantities.
	factor=sqrt(3*c**3*me**2*e**2/(16*Pi*epsilon0*h**3))	
	#We read the database to obtain the Einstein A coefficients in Hz.
	einsteinA=get_einstein_A_matrix(fine_states)
	#We read the database to obtain the transition frequencies in Hz.
	omega_fine=calculate_omega_matrix(fine_states)
	
	reduced_matrix_elements=[[0.0 for jj in range(len(fine_states))] for ii in range(len(fine_states))]
	
	for ii in range(len(fine_states)):
		i=fine_states[ii]
		for jj in range(ii):
			j=fine_states[jj]
			t=Transition(i,j)
			einsteinAij=einsteinA[ii][jj]
			omega0=omega_fine[ii][jj]
			
			#The formula is valid only for i =/= j so that omega0 =/= 0.
			#Because fine states are asumed to be ordered by their energies we can asume that
			# i decays in j.
			Ji=i.j; Jj=j.j

			rij=sqrt((2.0*Ji+1)/(2*Jj+1))*sqrt(einsteinAij/omega0**3)
			rij=factor*rij
		
			reduced_matrix_elements[ii][jj]=rij
			#We add the matrix elements on the other side of the diagonal.
			reduced_matrix_elements[jj][ii]=rij*(-1)**(Jj-Ji) * sqrt(2.0*Ji+1)/sqrt(2.0*Jj+1)
	
	return reduced_matrix_elements

def calculate_r_matrices_0(fine_states, reduced_matrix_elements):
	
	full_magnetic_states=make_list_of_states(fine_states,'magnetic',verbose=0)
	index_list_fine, index_list_hyperfine = calculate_boundaries(fine_states, full_magnetic_states)
	
	N_magnetic=len(full_magnetic_states)
	r=[]
	for p in [-1,0,1]:
		mat=[]
		for i in range(N_magnetic):
			row=[]
			ISO,N,L,J,F,M=full_magnetic_states[i].quantum_numbers
			ii=fine_index(i, index_list_fine)
			for j in range(N_magnetic):
				ISOp,Np,Lp,Jp,Fp,Mp=full_magnetic_states[j].quantum_numbers
				jj=fine_index(j, index_list_fine)
				
				red=reduced_matrix_elements[ii][jj]
				qu=quaver(ISO,p, J, F, M, Jp, Fp, Mp,numeric=True)
				row+=[red*qu]
			mat+=[row]
		r+=[mat]
	return r

def calculate_r_matrices_steck(fine_states, reduced_matrix_elements):
	
	magnetic_states=make_list_of_states(fine_states,'magnetic',verbose=0)
	index_list_fine, index_list_hyperfine = calculate_boundaries(fine_states, magnetic_states)
	
	Ne=len(magnetic_states)
		
	r=[[[0.0 for j in range(Ne)] for i in range(Ne)] for p in range(3)]
	
	if fine_states[0].isotope==85:
		II=Integer(5)/Integer(2)
	elif fine_states[0].isotope==87:
		II=Integer(3)/Integer(2)

	for p in [-1,0,1]:
		for i in range(Ne):
			ei=magnetic_states[i]
			ii=fine_index(i, index_list_fine)
			
			for j in range(Ne):
				ej=magnetic_states[j]
				jj=fine_index(j, index_list_fine)
				
				reduced_matrix_elementij=reduced_matrix_elements[ii][jj]
				if reduced_matrix_elementij != 0:
					
					ji=ei.j; jj=ej.j
					fi=ei.f; fj=ej.f
					mi=ei.m; mj=ej.m
										
					rpij =(-1)**(fj-1-mi)
					rpij*=sqrt(2*fi+1)
					rpij*=wigner_3j(fj,1,fi,mj,p,-mi)
					
					rpij*=(-1)**(fj+ji+1+II)
					rpij*=sqrt(2*fj+1)
					rpij*=sqrt(2*ji +1)
					rpij*=wigner_6j(ji,jj,1,fj,fi,II)
					
					rpij*=reduced_matrix_elementij
					
					r[p+1][i][j]=float(rpij)
	return r


if sage_included:
	var('S P D F')
else:
	S='S'
	P='P'
	D='D'
	F='F'

k_B=1.3806488e-23 #Boltzman's constant in J/K
uma=1.660538782e-27#Atomic mass unit in kg
m_Rb85=84.911789732*uma#Rb85 mass in kg
m_Rb87=86.909180520*uma#Rb87 mass in kg
