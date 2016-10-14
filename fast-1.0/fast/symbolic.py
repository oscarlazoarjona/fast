# -*- coding: utf-8 -*-
#Oscar Gerardo Lazo Arjona

from sympy import Symbol,Matrix,symbols
from sympy import I,conjugate
from sympy import sin,cos,exp,sqrt,pi
from sympy import pprint
from sympy import simplify
from sympy import KroneckerDelta

def define_density_matrix(Ne,explicitly_hermitian=False,normalized=False):
    if Ne>9:
		comma=","
		name=r"\rho"
		open_brace="}"
		close_brace="}"
    else:
		comma=""
		name="rho"
		open_brace=""
		close_brace=""

    rho=[]; re_rho=[]; im_rho=[]
    for i in range(Ne):
        row_rho=[]; row_re_rho=[]; row_im_rho=[]
        for j in range(Ne):
            if i==j:
                row_rho   +=[ Symbol(            name+open_brace+str(i+1)+comma+str(j+1)+close_brace,positive=True )]
                row_re_rho+=[ Symbol(            name+open_brace+str(i+1)+comma+str(j+1)+close_brace,positive=True )]
                row_im_rho+=[ 0 ]
            elif i>j:
                row_rho   +=[ Symbol(            name+open_brace+str(i+1)+comma+str(j+1)+close_brace               )]
                row_re_rho+=[ Symbol(r"\mathfrak{R}\rho_{"+str(i+1)+comma+str(j+1)+"}",    real=True )]
                row_im_rho+=[ Symbol(r"\mathfrak{I}\rho_{"+str(i+1)+comma+str(j+1)+"}",    real=True )]
            else:
                if explicitly_hermitian:
                    row_rho   +=[ conjugate(Symbol(            name+open_brace+str(j+1)+comma+str(i+1)+close_brace               ))]
                    row_re_rho+=[           Symbol(r"\mathfrak{R}\rho_{"+str(j+1)+comma+str(i+1)+"}",    real=True ) ]
                    row_im_rho+=[          -Symbol(r"\mathfrak{I}\rho_{"+str(j+1)+comma+str(i+1)+"}",    real=True ) ]
                else:
                    row_rho   +=[           Symbol(            name+open_brace+str(i+1)+comma+str(j+1)+close_brace              )]
                    row_re_rho+=[           Symbol(r"\mathfrak{R}\rho_{"+str(i+1)+comma+str(j+1)+"}",    real=True )]
                    row_im_rho+=[           Symbol(r"\mathfrak{I}\rho_{"+str(i+1)+comma+str(j+1)+"}",    real=True )]
        
        rho+=[row_rho]; re_rho+=[row_re_rho]; im_rho+=[row_im_rho]
    
    if normalized:
        rho11=1-sum([ rho[i][i] for i in range(1,Ne)])
        rho[0][0]   =rho11
        re_rho[0][0]=rho11
    
    
    rho   =Matrix(   rho)
    re_rho=Matrix(re_rho)
    im_rho=Matrix(im_rho)
    return rho,re_rho,im_rho

def define_laser_variables(Nl):
    E0         =[Symbol(  r"E_0^"+str(l+1),    real=True ) for l in range(Nl)]
    omega_laser=[Symbol(r"omega^"+str(l+1),    real=True ) for l in range(Nl)]
    return E0,omega_laser

def polarization_vector(phi,theta,alpha,beta,p):
    epsilon=Matrix([cos(2*beta),p*I*sin(2*beta),0])
    
    R1=Matrix([[ cos(2*alpha), -sin(2*alpha), 0         ],
               [ sin(2*alpha),  cos(2*alpha), 0         ],
               [  0          ,   0          , 1         ]])

    R2=Matrix([[ cos(theta)  ,  0           , sin(theta) ],
               [  0          ,  1           , 0          ],
               [-sin(theta)  ,  0           , cos(theta) ]])

    R3=Matrix([[ cos(phi)    , -sin(phi)    , 0         ],
               [ sin(phi)    ,  cos(phi)    , 0         ],
               [  0          ,   0          , 1         ]])
 
    return R3*R2*R1*epsilon

def cartesian_to_helicity(vector):
    return Matrix([ (vector[0]-I*vector[1])/sqrt(2), vector[2], -(vector[0]+I*vector[1])/sqrt(2) ])
def helicity_to_cartesian(vector):
    return Matrix([(vector[0]-vector[2])/sqrt(2), I*(vector[0]+vector[2])/sqrt(2), vector[1]])
def helicity_dot_product(v1,v2):
    return -v1[2]*v2[0] +v1[1]*v2[1]-v1[0]*v2[2]
def cartesian_dot_product(v1,v2):
    return  v1[0]*v2[0] +v1[1]*v2[1]+v1[2]*v2[2]

    
def define_r_components(Ne,explicitly_hermitian=False,helicity=False,real=True):
    if Ne>9: comma=","
    else: comma=""
    
    if helicity:
        names=["r_{-1;","r_{0;","r_{+1;"]
    else:
        names=["x","y","z"]
    
    r=[]
    if helicity:
        for p in range(3):
            r_comp=[]
            for i in range(Ne):
                r_row=[]
                for j in range(Ne):
                    if i==j:
                        r_row   +=[ 0 ]
                    elif i>j:
                        r_row   +=[           Symbol( names[p  ]+str(i+1)+comma+str(j+1)+"}" ,real=real) ]
                    elif explicitly_hermitian:
                        sign=int((-1)**(p-1))
                        r_row   +=[sign*conjugate(Symbol( names[2-p]+str(j+1)+comma+str(i+1)+"}" ,real=real))]
                    else:
                        r_row   +=[           Symbol( names[p  ]+str(i+1)+comma+str(j+1)+"}" ,real=real) ]
                r_comp+=[r_row]
            r_comp=Matrix(r_comp)
            r+=[r_comp]

    else:
        for p in range(3):
            r_comp=[]
            for i in range(Ne):
                r_row=[]
                for j in range(Ne):
                    if i==j:
                        r_row   +=[ 0 ]
                    elif i>j:
                        r_row   +=[           Symbol( names[p]+r"_{"+str(i+1)+comma+str(j+1)+"}" ,real=real) ]
                    elif explicitly_hermitian:
                        r_row   +=[ conjugate(Symbol( names[p]+r"_{"+str(j+1)+comma+str(i+1)+"}" ,real=real))]
                    else:
                        r_row   +=[           Symbol( names[p]+r"_{"+str(i+1)+comma+str(j+1)+"}" ,real=real) ]
                r_comp+=[r_row]
            r_comp=Matrix(r_comp)
            r+=[r_comp]
            
    return r

def define_frequencies(Ne,explicitly_antisymmetric=False):
	
	omega_level=[Symbol('omega_'+str(i+1),real=True) for i in range(Ne)]
	
	if Ne>9:
		comma=","
		open_brace= "{"
		close_brace="}"
	else:
		comma=""
		open_brace= ""
		close_brace=""

	omega=[]; gamma=[]
	for i in range(Ne):
		row_omega=[]; row_gamma=[]
		for j in range(Ne):
			if i==j:
				om=0; ga=0
			elif i>j:
				om= Symbol(r"omega_"+open_brace+str(i+1)+comma+str(j+1)+close_brace,real=True)
				ga= Symbol(r"gamma_"+open_brace+str(i+1)+comma+str(j+1)+close_brace,real=True)
			elif explicitly_antisymmetric:
				om=-Symbol(r"omega_"+open_brace+str(j+1)+comma+str(i+1)+close_brace,real=True)
				ga=-Symbol(r"gamma_"+open_brace+str(j+1)+comma+str(i+1)+close_brace,real=True)
			else:
				om= Symbol(r"omega_"+open_brace+str(i+1)+comma+str(j+1)+close_brace,real=True)
				ga= Symbol(r"gamma_"+open_brace+str(i+1)+comma+str(j+1)+close_brace,real=True)
				

			row_omega+=[om]
			row_gamma+=[ga]
			
		omega+=[row_omega]
		gamma+=[row_gamma]
		
	omega=Matrix(omega)
	gamma=Matrix(gamma)
	
	return omega_level,omega,gamma

def delta_greater(i,j):
    if i>j: return 1
    else: return 0

def delta_lesser(i,j):
    if i<j: return 1
    else: return 0

def bra(i,Ne):
    if i not in range(1,Ne+1):
        raise ValueError,"i must be in [1 .. Ne]."
    return Matrix([KroneckerDelta(i-1,j) for j in range(Ne)]).transpose()
def ket(i,Ne):
    if i not in range(1,Ne+1):
        raise ValueError,"i must be in [1 .. Ne]."
    return Matrix([KroneckerDelta(i-1,j) for j in range(Ne)])

def lindblad_operator(A,rho):
    return A*rho*A.adjoint() - (A.adjoint()*A*rho + rho*A.adjoint()*A)/2

def lindblad_terms(gamma,rho,Ne):
    L=zeros(Ne)
    for i in range(Ne):
        for j in range(i):
            L += gamma[i,j]*lindblad_operator(ket(j+1,Ne)*bra(i+1,Ne),rho)
    return L
