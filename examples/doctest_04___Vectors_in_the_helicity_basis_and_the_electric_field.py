# -*- coding: utf-8 -*-
# Copyright (C) 2017 Oscar Gerardo Lazo Arjona
# mailto: oscar.lazo@correo.nucleares.unam.mx

__doc__ = r"""

# The helicity basis and the dot product

>>> from sympy import init_printing
>>> init_printing() 
>>> print_ascii=True
>>> #print_ascii=False

>>> from sympy import Symbol,Matrix,symbols
>>> from sympy import I,conjugate
>>> from sympy import sin,cos,exp,sqrt,pi
>>> from sympy import pprint
>>> from sympy import simplify
    
>>> from fast.symbolic import define_laser_variables, polarization_vector
>>> from fast.symbolic import cartesian_to_helicity, helicity_to_cartesian, helicity_dot_product
>>> from fast.symbolic import define_r_components
>>> from fast.misc import fprint

Two vectors defined through their components in the cartesian basis,

>>> a=Matrix(symbols("a_x a_y a_z"))
>>> b=Matrix(symbols("b_x b_y b_z"))
>>> fprint([a,b], print_ascii=print_ascii)
 [a_x]  [b_x] 
 [   ]  [   ] 
[[a_y], [b_y]]
 [   ]  [   ] 
 [a_z]  [b_z] 



have their dot product is defined as

>>> fprint( a.dot(b),  print_ascii=print_ascii)
a_x*b_x + a_y*b_y + a_z*b_z



The vectors can be expressed through their components $a_{-1},\\ a_{0},\\ a_{+1}$ and $b_{-1},\\ b_{0},\\ b_{+1}$ in the helicity basis as

>>> a_helicity=cartesian_to_helicity(a)
>>> b_helicity=cartesian_to_helicity(b)
    
>>> fprint( [a_helicity, b_helicity], print_ascii=print_ascii)
 [  ___               ]  [  ___               ] 
 [\/ 2 *(a_x - I*a_y) ]  [\/ 2 *(b_x - I*b_y) ] 
 [------------------- ]  [------------------- ] 
 [         2          ]  [         2          ] 
 [                    ]  [                    ] 
[[        a_z         ], [        b_z         ]]
 [                    ]  [                    ] 
 [  ___               ]  [  ___               ] 
 [\/ 2 *(-a_x - I*a_y)]  [\/ 2 *(-b_x - I*b_y)] 
 [--------------------]  [--------------------] 
 [         2          ]  [         2          ] 



The dot product of two vectors in the helicity basis can be prooved to be $\\vec{a}\\cdot\\vec{b}=- a_{-1}b_{+1} +a_{0}   b_{0} -a_{+1} b_{-1}$

>>> fprint( helicity_dot_product(a_helicity,b_helicity), print_ascii=print_ascii)
          (-a_x - I*a_y)*(b_x - I*b_y)   (a_x - I*a_y)*(-b_x - I*b_y)
a_z*b_z - ---------------------------- - ----------------------------
                       2                              2              



>>> fprint( simplify( helicity_dot_product(a_helicity,b_helicity) ), print_ascii=print_ascii)
a_x*b_x + a_y*b_y + a_z*b_z



# The electric field
We define a few important symbols

>>> t,c=symbols("t c",positive=True)
>>> X,Y,Z=symbols("X Y Z",real=True)
>>> R=Matrix([X,Y,Z])
>>> fprint( [t,c,R], print_ascii=print_ascii)
       [X] 
       [ ] 
[t, c, [Y]]
       [ ] 
       [Z] 



We will specify the electric field associated to a plane wave with arbitrary amplitude and frequency

>>> E0,omega_laser=define_laser_variables(1)
>>> fprint( [E0,omega_laser], print_ascii=print_ascii)
[[E_{01}], [varpi_1]]



propagating through an arbitrary wave vector $\\vec{k}$

>>> phi,theta,alpha,beta=symbols("phi theta alpha beta")
    
>>> k=omega_laser[0]/c*Matrix([cos(phi)*sin(theta),sin(phi)*sin(theta),cos(theta)])
>>> fprint(k, print_ascii=print_ascii)
[varpi_1*sin(theta)*cos(phi)]
[---------------------------]
[             c             ]
[                           ]
[varpi_1*sin(phi)*sin(theta)]
[---------------------------]
[             c             ]
[                           ]
[    varpi_1*cos(theta)     ]
[    ------------------     ]
[            c              ]



with an arbitrary polarization,

>>> ep=polarization_vector(phi,theta,alpha,beta, 1)
>>> em=polarization_vector(phi,theta,alpha,beta,-1)
>>> fprint(ep, print_ascii=print_ascii)
[(-sin(2*alpha)*sin(phi) + cos(2*alpha)*cos(phi)*cos(theta))*cos(2*beta) + I*(-sin(2*alpha)*cos(phi)*cos(theta) - sin(
[                                                                                                                     
[(sin(2*alpha)*cos(phi) + sin(phi)*cos(2*alpha)*cos(theta))*cos(2*beta) + I*(-sin(2*alpha)*sin(phi)*cos(theta) + cos(2
[                                                                                                                     
[                                    I*sin(2*alpha)*sin(2*beta)*sin(theta) - sin(theta)*cos(2*alpha)*cos(2*beta)      
<BLANKLINE>
phi)*cos(2*alpha))*sin(2*beta)]
                              ]
*alpha)*cos(phi))*sin(2*beta) ]
                              ]
                              ]



>>> fprint(em, print_ascii=print_ascii)
[(-sin(2*alpha)*sin(phi) + cos(2*alpha)*cos(phi)*cos(theta))*cos(2*beta) - I*(-sin(2*alpha)*cos(phi)*cos(theta) - sin(
[                                                                                                                     
[(sin(2*alpha)*cos(phi) + sin(phi)*cos(2*alpha)*cos(theta))*cos(2*beta) - I*(-sin(2*alpha)*sin(phi)*cos(theta) + cos(2
[                                                                                                                     
[                                   -I*sin(2*alpha)*sin(2*beta)*sin(theta) - sin(theta)*cos(2*alpha)*cos(2*beta)      
<BLANKLINE>
phi)*cos(2*alpha))*sin(2*beta)]
                              ]
*alpha)*cos(phi))*sin(2*beta) ]
                              ]
                              ]



The electric field is given by

>>> arg=k.dot(R)-omega_laser[0]*t
>>> E=E0[0]/2*(ep*exp(+arg) + em*exp( -arg))

Which can be seen to be orthogonal to the propagation direction.

>>> print simplify(E.dot(k))
0



We can simplify this expression by noting that in the cases of our interest the wavelengths of the fields (780 nm and 776 nm) are much longer than the atomic radius (248 pm for rubidium) the spacial variation of the field is relatively insignificant, so a single point to evaluate this field can be taken at $\\vec{R}=0$, and thus the field can be taken as

>>> arg=-omega_laser[0]*t
>>> E=E0[0]/2*(ep*exp(+arg) + em*exp( -arg))

# The position operator in the helicity basis.
We can write think of the position operator as a vector of operators $\\vec{\\hat{r}}$ that act on a Hilbert space of a given dimension. In this case we will use dimension 2.

>>> r_cartesian=define_r_components(2)
>>> fprint(r_cartesian, print_ascii=print_ascii)
 [  0     x_{12}]  [  0     y_{12}]  [  0     z_{12}] 
[[              ], [              ], [              ]]
 [x_{21}    0   ]  [y_{21}    0   ]  [z_{21}    0   ] 



>>> r_helicity=define_r_components(2,helicity=True)
>>> fprint(r_helicity, print_ascii=print_ascii)
 [    0      r_{-1;12}]  [   0      r_{0;12}]  [    0      r_{+1;12}] 
[[                    ], [                  ], [                    ]]
 [r_{-1;21}      0    ]  [r_{0;21}     0    ]  [r_{+1;21}      0    ] 



>>> r_helicity21=Matrix([r_helicity[0][1,0],r_helicity[1][1,0],r_helicity[2][1,0]])
>>> fprint(r_helicity21, print_ascii=print_ascii)
[r_{-1;21}]
[         ]
[r_{0;21} ]
[         ]
[r_{+1;21}]



We take $\\vec{r}_{ij}$ in the helicity basis, and transforming it to the cartesian and taking the complex conjugate basis we get $\\vec{r}_{ji}$, which in turn can be taken back to the helicity basis to obtain $\\vec{r}_{ji}$ in terms of $\\vec{r}_{ij}$ in the helicity basis.

>>> fprint( simplify( cartesian_to_helicity(conjugate(helicity_to_cartesian(r_helicity21))) ), print_ascii=print_ascii)
[-r_{+1;21}]
[          ]
[ r_{0;21} ]
[          ]
[-r_{-1;21}]



Check whether this actually makes $\\vec{\\hat{r}}$ hermitian.

>>> r_helicity=define_r_components(2,helicity=True,explicitly_hermitian=True)
>>> fprint(r_helicity, print_ascii=print_ascii)
 [    0      -r_{+1;21}]  [   0      r_{0;21}]  [    0      -r_{-1;21}] 
[[                     ], [                  ], [                     ]]
 [r_{-1;21}      0     ]  [r_{0;21}     0    ]  [r_{+1;21}      0     ] 



>>> r_helicity21=Matrix([r_helicity[0][1,0],r_helicity[1][1,0],r_helicity[2][1,0]])
>>> r_helicity12=Matrix([r_helicity[0][0,1],r_helicity[1][0,1],r_helicity[2][0,1]])
>>> fprint( [r_helicity21,r_helicity12], print_ascii=print_ascii)
 [r_{-1;21}]  [-r_{+1;21}] 
 [         ]  [          ] 
[[r_{0;21} ], [ r_{0;21} ]]
 [         ]  [          ] 
 [r_{+1;21}]  [-r_{-1;21}] 



>>> fprint( simplify( helicity_to_cartesian(r_helicity21)-helicity_to_cartesian(r_helicity12).conjugate() ), print_ascii=print_ascii)
[0]
[ ]
[0]
[ ]
[0]



So yes, this makes the $\\vec{\\hat{r}}$ operator hermitian.

[]

"""
__doc__=__doc__.replace("+IGNORE_PLOT_STEP1", "+ELLIPSIS\n[<...>]")
__doc__=__doc__.replace("+IGNORE_PLOT_STEP2", "+ELLIPSIS\n<...>")
__doc__=__doc__.replace("+IGNORE_PLOT_STEP3", "+ELLIPSIS\n(...)")
__doc__=__doc__.replace("+IGNORE_PLOT_STEP5", "+ELLIPSIS\nText(...)")
__doc__=__doc__.replace("+IGNORE_PLOT_STEP4", "\n")
