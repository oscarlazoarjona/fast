# -*- coding: utf-8 -*-
# Copyright (C) 2017 Oscar Gerardo Lazo Arjona
# mailto: oscar.lazoarjona@physics.ox.ac.uk

__doc__ = r"""

>>> from fast import *
>>> init_printing()
>>> use_unicode=True; use_unicode=False

>>> from sympy import sin,cos,exp,sqrt,pi,zeros,I
>>> from numpy import array

We define the number of states and of radiation fields.

>>> Ne=3
>>> Nl=2

We define the variables related to the laser field.

>>> E0,omega_laser=define_laser_variables(Nl)
>>> pprint(E0,use_unicode=use_unicode)
[E_0^1, E_0^2]



>>> pprint(omega_laser,use_unicode=use_unicode)
[omega^1, omega^2]



We define a few important symbols.

>>> t,hbar,e=symbols("t hbar e",positive=True)
>>> pprint([t,hbar,e],use_unicode=use_unicode)
[t, hbar, e]



We write two electric fields propagating trough the $\\hat{x}$ direction polarized in the $\\hat{z}$ direction. First the wave vectors:

>>> phi=0; theta=pi/2; alpha=pi/2; beta=0
    
>>> k1=Matrix([cos(phi)*sin(theta),sin(phi)*sin(theta),cos(theta)])
>>> k2=Matrix([cos(phi)*sin(theta),sin(phi)*sin(theta),cos(theta)])
    
>>> pprint([k1,k2],use_unicode=use_unicode)
[[1], [1]]
 [ ]  [ ] 
 [0]  [0] 
 [ ]  [ ] 
 [0]  [0] 



The polarization vectors.

>>> ep1=polarization_vector(phi,theta,alpha,beta, 1)
>>> ep2=polarization_vector(phi,theta,alpha,beta, 1)
    
>>> em1=ep1.conjugate()
>>> em2=ep2.conjugate()
    
>>> ep=[ep1,ep2]
>>> em=[em1,em2]
    
>>> pprint([ep,em],use_unicode=use_unicode)
[[[0], [0]], [[0], [0]]]
  [ ]  [ ]    [ ]  [ ]  
  [0]  [0]    [0]  [0]  
  [ ]  [ ]    [ ]  [ ]  
  [1]  [1]    [1]  [1]  



The electric field (evaluated in $\\vec{R}=0$).

>>> zero_vect=Matrix([0,0,0])
>>> E_cartesian = [(E0[l]*ep[l]*exp(-I*omega_laser[l]*t) + E0[l].conjugate()*em[l]*exp( I*omega_laser[l]*t))/2 
...                     for l in range(Nl)]
    
>>> pprint(E_cartesian,use_unicode=use_unicode)
[[                   0                    ], [                   0            
 [                                        ]  [                                
 [                   0                    ]  [                   0            
 [                                        ]  [                                
 [       -I*omega^1*t    I*omega^1*t _____]  [       -I*omega^2*t    I*omega^2
 [E_0^1*e               e           *E_0^1]  [E_0^2*e               e         
 [------------------- + ------------------]  [------------------- + ----------
 [         2                    2         ]  [         2                    2 
<BLANKLINE>
        ]]
        ] 
        ] 
        ] 
*t _____] 
  *E_0^2] 
--------] 
        ] 



We write the electric field in the helicity basis.

>>> E=[cartesian_to_helicity(E_cartesian[l]) for l in range(Nl)]
>>> pprint(E,use_unicode=use_unicode)
[[                   0                    ], [                   0            
 [                                        ]  [                                
 [       -I*omega^1*t    I*omega^1*t _____]  [       -I*omega^2*t    I*omega^2
 [E_0^1*e               e           *E_0^1]  [E_0^2*e               e         
 [------------------- + ------------------]  [------------------- + ----------
 [         2                    2         ]  [         2                    2 
 [                                        ]  [                                
 [                   0                    ]  [                   0            
<BLANKLINE>
        ]]
        ] 
*t _____] 
  *E_0^2] 
--------] 
        ] 
        ] 
        ] 



We define the position operator.

>>> r=define_r_components(Ne,helicity=True,explicitly_hermitian=True)
>>> #Ladder means that r_{p;31}=0
>>> r=[ri.subs({r[0][2,0]:0,r[1][2,0]:0,r[2][2,0]:0}) for ri in r]
>>> pprint(r,num_columns=150,use_unicode=use_unicode)
[[    0      -r_{+1;21}      0     ], [   0      r_{0;21}     0    ], [    0      -r_{-1;21}      0     ]]
 [                                 ]  [                            ]  [                                 ] 
 [r_{-1;21}      0       -r_{+1;32}]  [r_{0;21}     0      r_{0;32}]  [r_{+1;21}      0       -r_{-1;32}] 
 [                                 ]  [                            ]  [                                 ] 
 [    0      r_{-1;32}       0     ]  [   0      r_{0;32}     0    ]  [    0      r_{+1;32}       0     ] 



The frequencies of the energy levels, the resonant frequencies, and the decay frequencies.

>>> omega_level,omega,gamma=define_frequencies(Ne,explicitly_antisymmetric=True)
>>> #Ladder means gamma31=0
>>> gamma=gamma.subs({gamma[2,0]:0})
    
>>> pprint(omega_level,use_unicode=use_unicode)
[omega_1, omega_2, omega_3]



>>> pprint(omega,use_unicode=use_unicode)
[   0      -omega_21  -omega_31]
[                              ]
[omega_21      0      -omega_32]
[                              ]
[omega_31  omega_32       0    ]



>>> pprint(gamma,use_unicode=use_unicode)
[   0      -gamma_21      0    ]
[                              ]
[gamma_21      0      -gamma_32]
[                              ]
[   0      gamma_32       0    ]



The atomic hamiltonian is

>>> H0=Matrix([[hbar*omega_level[i]*KroneckerDelta(i,j) for j in range(Ne)] for i in range(Ne)])
>>> pprint(H0,use_unicode=use_unicode)
[hbar*omega_1       0             0      ]
[                                        ]
[     0        hbar*omega_2       0      ]
[                                        ]
[     0             0        hbar*omega_3]



The interaction hamiltonian is

>>> zero_matrix=zeros(Ne,Ne)
>>> H1=sum([ e*helicity_dot_product(E[l],r) for l in range(Nl)],zero_matrix)

and the complete hamiltonian is

>>> H=H0+H1

# Rotating wave approximation
Notice that the electric field can be separated by terms with positive and negative frequency:

>>> E_cartesian_p=[E0[l]/2*ep[l]*exp(-I*omega_laser[l]*t) for l in range(Nl)]
>>> E_cartesian_m=[E0[l].conjugate()/2*em[l]*exp( I*omega_laser[l]*t) for l in range(Nl)]
    
>>> E_p=[cartesian_to_helicity(E_cartesian_p[l]) for l in range(Nl)]
>>> E_m=[cartesian_to_helicity(E_cartesian_m[l]) for l in range(Nl)]
    
>>> pprint([E_p,E_m],use_unicode=use_unicode)
[[[         0         ], [         0         ]], [[        0         ], [     
  [                   ]  [                   ]    [                  ]  [     
  [       -I*omega^1*t]  [       -I*omega^2*t]    [ I*omega^1*t _____]  [ I*om
  [E_0^1*e            ]  [E_0^2*e            ]    [e           *E_0^1]  [e    
  [-------------------]  [-------------------]    [------------------]  [-----
  [         2         ]  [         2         ]    [        2         ]  [     
  [                   ]  [                   ]    [                  ]  [     
  [         0         ]  [         0         ]    [        0         ]  [     
<BLANKLINE>
   0         ]]]
             ]  
ega^2*t _____]  
       *E_0^2]  
-------------]  
   2         ]  
             ]  
   0         ]  



>>> pprint( simplify(sum([E[l] for l in range(Nl)],zero_vect)-(sum([E_p[l]+E_m[l] for l in range(Nl)],zero_vect) )),use_unicode=use_unicode)
[0]
[ ]
[0]
[ ]
[0]



The position operator can also be separated in this way. We go to the interaction picture (with $\\hat{H}_0$ as the undisturbed hamiltonian)

>>> r_I=[ Matrix([[exp(I*omega[i,j]*t)*r[p][i,j] for j in range(Ne)] for i in range(Ne)]) for p in range(3)]
>>> pprint(r_I[0],use_unicode=use_unicode)
[                                     -I*omega_21*t                           
[           0             -r_{+1;21}*e                           0            
[                                                                             
[           I*omega_21*t                                         -I*omega_32*t
[r_{-1;21}*e                          0              -r_{+1;32}*e             
[                                                                             
[                                     I*omega_32*t                            
[           0              r_{-1;32}*e                           0            
<BLANKLINE>
]
]
]
]
]
]
]
]



>>> pprint(r_I[1],use_unicode=use_unicode)
[                                  -I*omega_21*t                         ]
[          0             r_{0;21}*e                          0           ]
[                                                                        ]
[          I*omega_21*t                                     -I*omega_32*t]
[r_{0;21}*e                         0             r_{0;32}*e             ]
[                                                                        ]
[                                  I*omega_32*t                          ]
[          0             r_{0;32}*e                          0           ]



>>> pprint(r_I[2],use_unicode=use_unicode)
[                                     -I*omega_21*t                           
[           0             -r_{-1;21}*e                           0            
[                                                                             
[           I*omega_21*t                                         -I*omega_32*t
[r_{+1;21}*e                          0              -r_{-1;32}*e             
[                                                                             
[                                     I*omega_32*t                            
[           0              r_{+1;32}*e                           0            
<BLANKLINE>
]
]
]
]
]
]
]
]



Which can be decomposed in positive and negative frequencies as

>>> r_I_p=[ Matrix([[ delta_greater(j,i)*exp(-I*omega[j,i]*t)*r[p][i,j] for j in range(Ne)]for i in range(Ne)]) for p in range(3)]
>>> pprint(r_I_p[0],use_unicode=use_unicode)
[               -I*omega_21*t                           ]
[0  -r_{+1;21}*e                           0            ]
[                                                       ]
[                                          -I*omega_32*t]
[0              0              -r_{+1;32}*e             ]
[                                                       ]
[0              0                          0            ]



>>> pprint(r_I_p[1],use_unicode=use_unicode)
[             -I*omega_21*t                         ]
[0  r_{0;21}*e                          0           ]
[                                                   ]
[                                      -I*omega_32*t]
[0             0             r_{0;32}*e             ]
[                                                   ]
[0             0                        0           ]



>>> pprint(r_I_p[2],use_unicode=use_unicode)
[               -I*omega_21*t                           ]
[0  -r_{-1;21}*e                           0            ]
[                                                       ]
[                                          -I*omega_32*t]
[0              0              -r_{-1;32}*e             ]
[                                                       ]
[0              0                          0            ]



>>> r_I_m=[ Matrix([[ delta_lesser( j,i)*exp( I*omega[i,j]*t)*r[p][i,j] for j in range(Ne)]for i in range(Ne)]) for p in range(3)]
>>> pprint(r_I_m[0],use_unicode=use_unicode)
[           0                        0             0]
[                                                   ]
[           I*omega_21*t                            ]
[r_{-1;21}*e                         0             0]
[                                                   ]
[                                    I*omega_32*t   ]
[           0             r_{-1;32}*e              0]



>>> pprint(r_I_m[1],use_unicode=use_unicode)
[          0                       0             0]
[                                                 ]
[          I*omega_21*t                           ]
[r_{0;21}*e                        0             0]
[                                                 ]
[                                  I*omega_32*t   ]
[          0             r_{0;32}*e              0]



>>> pprint(r_I_m[2],use_unicode=use_unicode)
[           0                        0             0]
[                                                   ]
[           I*omega_21*t                            ]
[r_{+1;21}*e                         0             0]
[                                                   ]
[                                    I*omega_32*t   ]
[           0             r_{+1;32}*e              0]



that summed equal $\\vec{\\hat{r}}_I$

>>> pprint( [r_I[p]-(r_I_p[p]+r_I_m[p]) for p in range(3)] ,use_unicode=use_unicode)
[[0  0  0], [0  0  0], [0  0  0]]
 [       ]  [       ]  [       ] 
 [0  0  0]  [0  0  0]  [0  0  0] 
 [       ]  [       ]  [       ] 
 [0  0  0]  [0  0  0]  [0  0  0] 



Thus the interaction hamiltonian in the interaciton picture is
\\begin{equation}
    \\hat{H}_{1I}=e\\vec{E}\\cdot \\vec{\\hat{r}}_I= e(\\vec{E}^{(+)}\\cdot \\vec{\\hat{r}}^{(+)}_I + \\vec{E}^{(+)}\\cdot \\vec{\\hat{r}}^{(-)}_I + \\vec{E}^{(-)}\\cdot \\vec{\\hat{r}}^{(+)}_I + \\vec{E}^{(-)}\\cdot \\vec{\\hat{r}}^{(-)}_I)
\\end{equation}

>>> H1I=sum([ e*helicity_dot_product(E[l],r_I) for l in range(Nl)],zero_matrix)

Since both $\\omega^l$ and $\\omega_{ij}$ are in the order of THz, the terms that have frequencies with the same sign are summed, and thus also of the order of THz. The frequencies in the terms with oposite signs however, are detunings of the order of MHz. Since we are only interested in the coarse-grained evolution of the density matrix, we may omit the fast terms and approximate

\\begin{equation}
    \\hat{H}_{1I} \\simeq \\hat{H}_{1I,RWA}= e( \\vec{E}^{(+)}\\cdot \\vec{\\hat{r}}^{(-)}_I + \\vec{E}^{(-)}\\cdot \\vec{\\hat{r}}^{(+)}_I )
\\end{equation}

That is known as the rotating wave approximation (RWA).

>>> H1IRWA=sum( [ (e*(helicity_dot_product(E_p[l],r_I_m)+helicity_dot_product(E_m[l],r_I_p))) for l in range(Nl)],zero_matrix)

For instance, the $(\\hat{H}_{1I,RWA})_{21}$ element is

>>> pprint(H1IRWA[1,0].expand(),use_unicode=use_unicode)
                  -I*omega^1*t  I*omega_21*t                     -I*omega^2*t 
E_0^1*e*r_{0;21}*e            *e               E_0^2*e*r_{0;21}*e            *
-------------------------------------------- + -------------------------------
                     2                                              2         
<BLANKLINE>
 I*omega_21*t
e            
-------------
<BLANKLINE>



But if the detuning $\\omega_{21}-\\omega^1 \\ll \\omega_{21}-\\omega^2$ (the second field is far detuned from the $1 \\rightarrow 2$ transition), then $\\omega_{21}-\\omega^2$ may be also considered too high a frequency to be relevant to coarse-grained evolution. So we might neclect that term in $(\\hat{H}_{1I,RWA})_{21}$ and similarly neglect the $\\omega_{32}-\\omega^1$ for term in $(\\hat{H}_{1I,RWA})_{32}$:

>>> pprint(H1IRWA[2,1].expand(),use_unicode=use_unicode)
                  -I*omega^1*t  I*omega_32*t                     -I*omega^2*t 
E_0^1*e*r_{0;32}*e            *e               E_0^2*e*r_{0;32}*e            *
-------------------------------------------- + -------------------------------
                     2                                              2         
<BLANKLINE>
 I*omega_32*t
e            
-------------
<BLANKLINE>



In other words, if the detunings in our experiments allow the approximmation, we might choose which frequency components $\\omega^l$ excite which transitions. Let us say that $L_{ij}$ is the set of $l$ such that $\\omega^l$ excites the transition $i\\rightarrow j$

>>> Lij=[[1,2,[1]],[2,3,[2]]]
>>> Lij=formatLij(Lij,Ne)
>>> print array(Lij)
[[[] [1] []]
 [[1] [] [2]]
 [[] [2] []]]



Thus the interacion hamiltonian in the interaction picture can be approximated as

>>> H1IRWA =sum([ e*( helicity_dot_product( E_p[l],vector_element(r_I_m,i,j)) ) * ket(i+1,Ne)*bra(j+1,Ne) 
...             for l in range(Nl) for j in range(Ne) for i in range(Ne) if l+1 in Lij[i][j] ],zero_matrix)
>>> H1IRWA+=sum([ e*( helicity_dot_product( E_m[l],vector_element(r_I_p,i,j)) ) * ket(i+1,Ne)*bra(j+1,Ne) 
...             for l in range(Nl) for j in range(Ne) for i in range(Ne) if l+1 in Lij[i][j] ],zero_matrix)
    
>>> pprint(H1IRWA,num_columns=150,use_unicode=use_unicode)
[                                                          I*omega^1*t  -I*omega_21*t _____                                              ]
[                                              e*r_{0;21}*e           *e             *E_0^1                                              ]
[                     0                        --------------------------------------------                       0                      ]
[                                                                   2                                                                    ]
[                                                                                                                                        ]
[                  -I*omega^1*t  I*omega_21*t                                                            I*omega^2*t  -I*omega_32*t _____]
[E_0^1*e*r_{0;21}*e            *e                                                            e*r_{0;32}*e           *e             *E_0^2]
[--------------------------------------------                       0                        --------------------------------------------]
[                     2                                                                                           2                      ]
[                                                                                                                                        ]
[                                                                -I*omega^2*t  I*omega_32*t                                              ]
[                                              E_0^2*e*r_{0;32}*e            *e                                                          ]
[                     0                        --------------------------------------------                       0                      ]
[                                                                   2                                                                    ]



Returning to the Schrödinger picture we have.

>>> r_p=[ Matrix([[ delta_greater(j,i)*r[p][i,j] for j in range(Ne)]for i in range(Ne)]) for p in range(3)]
>>> pprint(r_p,num_columns=150,use_unicode=use_unicode)
[[0  -r_{+1;21}      0     ], [0  r_{0;21}     0    ], [0  -r_{-1;21}      0     ]]
 [                         ]  [                     ]  [                         ] 
 [0      0       -r_{+1;32}]  [0     0      r_{0;32}]  [0      0       -r_{-1;32}] 
 [                         ]  [                     ]  [                         ] 
 [0      0           0     ]  [0     0         0    ]  [0      0           0     ] 



>>> r_m=[ Matrix([[ delta_lesser( j,i)*r[p][i,j] for j in range(Ne)]for i in range(Ne)]) for p in range(3)]
>>> pprint(r_m,num_columns=150,use_unicode=use_unicode)
[[    0          0      0], [   0         0      0], [    0          0      0]]
 [                       ]  [                     ]  [                       ] 
 [r_{-1;21}      0      0]  [r_{0;21}     0      0]  [r_{+1;21}      0      0] 
 [                       ]  [                     ]  [                       ] 
 [    0      r_{-1;32}  0]  [   0      r_{0;32}  0]  [    0      r_{+1;32}  0] 



>>> pprint( [r[p]-(r_p[p]+r_m[p]) for p in range(3)] ,use_unicode=use_unicode)
[[0  0  0], [0  0  0], [0  0  0]]
 [       ]  [       ]  [       ] 
 [0  0  0]  [0  0  0]  [0  0  0] 
 [       ]  [       ]  [       ] 
 [0  0  0]  [0  0  0]  [0  0  0] 



Thus the interaction hamiltonian in the Schrödinger picture in the rotating wave approximation is

>>> H1RWA =sum([ e*( helicity_dot_product( E_p[l],vector_element(r_m,i,j)) ) * ket(i+1,Ne)*bra(j+1,Ne) 
...             for l in range(Nl) for j in range(Ne) for i in range(Ne) if l+1 in Lij[i][j] ],zero_matrix)
>>> H1RWA+=sum([ e*( helicity_dot_product( E_m[l],vector_element(r_p,i,j)) ) * ket(i+1,Ne)*bra(j+1,Ne) 
...             for l in range(Nl) for j in range(Ne) for i in range(Ne) if l+1 in Lij[i][j] ],zero_matrix)
    
>>> pprint(H1RWA,num_columns=150,use_unicode=use_unicode)
[                                            I*omega^1*t _____                                ]
[                                e*r_{0;21}*e           *E_0^1                                ]
[              0                 -----------------------------                 0              ]
[                                              2                                              ]
[                                                                                             ]
[                  -I*omega^1*t                                              I*omega^2*t _____]
[E_0^1*e*r_{0;21}*e                                              e*r_{0;32}*e           *E_0^2]
[------------------------------                0                 -----------------------------]
[              2                                                               2              ]
[                                                                                             ]
[                                                  -I*omega^2*t                               ]
[                                E_0^2*e*r_{0;32}*e                                           ]
[              0                 ------------------------------                0              ]
[                                              2                                              ]



And the complete hamiltonian in the Schrödinger picture in the rotating wave approximation is

>>> HRWA=H0+H1RWA
>>> pprint(HRWA,use_unicode=use_unicode)
[                                            I*omega^1*t _____                
[                                e*r_{0;21}*e           *E_0^1                
[         hbar*omega_1           -----------------------------                
[                                              2                              
[                                                                             
[                  -I*omega^1*t                                              I
[E_0^1*e*r_{0;21}*e                                              e*r_{0;32}*e 
[------------------------------           hbar*omega_2           -------------
[              2                                                              
[                                                                             
[                                                  -I*omega^2*t               
[                                E_0^2*e*r_{0;32}*e                           
[              0                 ------------------------------          hbar*
[                                              2                              
<BLANKLINE>
                ]
                ]
 0              ]
                ]
                ]
*omega^2*t _____]
          *E_0^2]
----------------]
 2              ]
                ]
                ]
                ]
omega_3         ]
                ]



# Rotating Frame
Next we will make a phase transformation in order to eliminate the explicit time dependance of the equations.

>>> c,ctilde,phase=define_psi_coefficients(Ne)
>>> pprint([c,ctilde,phase],use_unicode=use_unicode)
[[c1(t)], [\tilde{c}_{1}(t)], [theta1]]
 [     ]  [                ]  [      ] 
 [c2(t)]  [\tilde{c}_{2}(t)]  [theta2] 
 [     ]  [                ]  [      ] 
 [c3(t)]  [\tilde{c}_{3}(t)]  [theta3] 



>>> psi=Matrix([ exp(I*phase[i]*t)*ctilde[i] for i in range(Ne)])
>>> pprint(psi,use_unicode=use_unicode)
[                  I*t*theta1]
[\tilde{c}_{1}(t)*e          ]
[                            ]
[                  I*t*theta2]
[\tilde{c}_{2}(t)*e          ]
[                            ]
[                  I*t*theta3]
[\tilde{c}_{3}(t)*e          ]



The Schrödinger equation $i\\hbar \\partial_t |\\psi\\rangle=\\hat{H}_{RWA}$ is

>>> lhs=Matrix([(I*hbar*Derivative(psi[i],t).doit()).expand() for i in range(Ne)])
>>> pprint(lhs,use_unicode=use_unicode)
[                               I*t*theta1]
[-hbar*theta1*\tilde{c}_{1}(t)*e          ]
[                                         ]
[                               I*t*theta2]
[-hbar*theta2*\tilde{c}_{2}(t)*e          ]
[                                         ]
[                               I*t*theta3]
[-hbar*theta3*\tilde{c}_{3}(t)*e          ]



>>> rhs=HRWA*psi

We multiply each of these equations by $e^{-i \\theta_i t}$ and substracting $i \\theta_i \\tilde{c}_i$

>>> lhs_new=Matrix([simplify(  lhs[i]*exp(-I*phase[i]*t) +hbar*phase[i]*ctilde[i] ) for i in range(Ne)])
>>> pprint(lhs_new,use_unicode=use_unicode)
[0]
[ ]
[0]
[ ]
[0]



>>> rhs_new=Matrix([simplify(  rhs[i]*exp(-I*phase[i]*t) +hbar*phase[i]*ctilde[i] ) for i in range(Ne)])
>>> pprint(rhs_new,num_columns=120,use_unicode=use_unicode)
[                                                                  I*omega^1*t  -I*t*theta1  I*t*theta2 _____         
[                                     e*r_{0;21}*\tilde{c}_{2}(t)*e           *e           *e          *E_0^1         
[                                     ----------------------------------------------------------------------- + hbar*o
[                                                                        2                                            
[                                                                                                                     
[                                   -I*omega^1*t  I*t*theta1  -I*t*theta2                                I*omega^2*t  
[E_0^1*e*r_{0;21}*\tilde{c}_{1}(t)*e            *e          *e              e*r_{0;32}*\tilde{c}_{3}(t)*e           *e
[------------------------------------------------------------------------ + ------------------------------------------
[                                   2                                                                          2      
[                                                                                                                     
[                                                                        -I*omega^2*t  I*t*theta2  -I*t*theta3        
[                                     E_0^2*e*r_{0;32}*\tilde{c}_{2}(t)*e            *e          *e                   
[                                     ------------------------------------------------------------------------ + hbar*
[                                                                        2                                            
<BLANKLINE>
                                                                                            ]
                                                                                            ]
mega_1*\tilde{c}_{1}(t) + hbar*theta1*\tilde{c}_{1}(t)                                      ]
                                                                                            ]
                                                                                            ]
-I*t*theta2  I*t*theta3 _____                                                               ]
           *e          *E_0^2                                                               ]
----------------------------- + hbar*omega_2*\tilde{c}_{2}(t) + hbar*theta2*\tilde{c}_{2}(t)]
                                                                                            ]
                                                                                            ]
                                                                                            ]
                                                                                            ]
omega_3*\tilde{c}_{3}(t) + hbar*theta3*\tilde{c}_{3}(t)                                     ]
                                                                                            ]



It can be seen that the equations loose their explicit time dependance only if $\\omega^{1} - \\theta_{1} + \\theta_{2}=0$ and $\\omega^{2} - \\theta_{2} + \\theta_{3}=0$. Which is satisfied if

>>> phase_transformation=solve([omega_laser[0]-phase[0]+phase[1],omega_laser[1]-phase[1]+phase[2]],[phase[1],phase[2]],
...                            dict=True)[0]
>>> pprint(phase_transformation,use_unicode=use_unicode)
{theta2: -omega^1 + theta1, theta3: -omega^1 - omega^2 + theta1}



There is a free parameter $\\theta_1$, which is to be expected, since state vetors $|\\psi\\rangle$ always have a global phase invariance

>>> pprint(psi.subs(phase_transformation),use_unicode=use_unicode)
[                             I*t*theta1            ]
[           \tilde{c}_{1}(t)*e                      ]
[                                                   ]
[                       I*t*(-omega^1 + theta1)     ]
[     \tilde{c}_{2}(t)*e                            ]
[                                                   ]
[                  I*t*(-omega^1 - omega^2 + theta1)]
[\tilde{c}_{3}(t)*e                                 ]



Thus the equations become

>>> pprint(lhs_new,use_unicode=use_unicode)
[0]
[ ]
[0]
[ ]
[0]



>>> rhs_new=simplify(rhs_new.subs(phase_transformation)).expand()

It can be seen that this is the Schrödinger equation derived from an effective hamiltonian $\\tilde{H}$

>>> Htilde=Matrix([ [Derivative(rhs_new[i],ctilde[j]).doit() for j in range(Ne)] for i in range(Ne)])
>>> pprint(Htilde,use_unicode=use_unicode)
[                                                    _____                    
[                                         e*r_{0;21}*E_0^1                    
[hbar*omega_1 + hbar*theta1               ----------------                    
[                                                2                            
[                                                                             
[                                                                             
[     E_0^1*e*r_{0;21}                                                        
[     ----------------       -hbar*omega^1 + hbar*omega_2 + hbar*theta1       
[            2                                                                
[                                                                             
[                                         E_0^2*e*r_{0;32}                    
[            0                            ----------------               -hbar
[                                                2                            
<BLANKLINE>
                                                    ]
                                                    ]
                       0                            ]
                                                    ]
                                                    ]
                          _____                     ]
               e*r_{0;32}*E_0^2                     ]
               ----------------                     ]
                      2                             ]
                                                    ]
                                                    ]
*omega^1 - hbar*omega^2 + hbar*omega_3 + hbar*theta1]
                                                    ]



We can see that it is convenient to choose $\\theta_1=-\\omega_1$ to simplify the hamiltonian. Also, we can recognize $\\omega^1-\\omega_2+\\omega_1=\\delta^1$ as the detuning of the first field relative to the atomic transition $\\omega_{21}=\\omega_2-\\omega_1$, and the same for $\\omega^2-\\omega_3+\\omega_2=\\delta^2$.

>>> delta1,delta2=symbols("delta1 delta2",real=True)
>>> Htilde=Htilde.subs({phase[0]:-omega_level[0]})
>>> Htilde=Htilde.subs({omega_laser[0]:delta1+omega_level[1]-omega_level[0]})
>>> Htilde=Htilde.subs({omega_laser[1]:delta2+omega_level[2]-omega_level[1]})
    
>>> Htilde=Htilde.expand()
    
>>> pprint(Htilde,use_unicode=use_unicode)
[                             _____                            ]
[                  e*r_{0;21}*E_0^1                            ]
[       0          ----------------              0             ]
[                         2                                    ]
[                                                              ]
[                                                    _____     ]
[E_0^1*e*r_{0;21}                         e*r_{0;32}*E_0^2     ]
[----------------    -delta1*hbar         ----------------     ]
[       2                                        2             ]
[                                                              ]
[                  E_0^2*e*r_{0;32}                            ]
[       0          ----------------  -delta1*hbar - delta2*hbar]
[                         2                                    ]



If we define the Rabi frequencies $\\Omega_1 =e E_0^1 r_{0;21}/\\hbar$ and $\\Omega_2 =e E_0^2 r_{0;32}/\\hbar$

>>> Omega1,Omega2=symbols("Omega1 Omega2",real=True)
>>> Omega1,Omega2=symbols("Omega1 Omega2")
>>> Htilde=Htilde.subs({E0[0]:Omega1*hbar/r[1][1,0]/e})
>>> Htilde=Htilde.subs({E0[1]:Omega2*hbar/r[1][2,1]/e})
    
>>> pprint(Htilde,use_unicode=use_unicode)
[                  ______                             ]
[             hbar*Omega1                             ]
[     0       -----------               0             ]
[                  2                                  ]
[                                                     ]
[                                       ______        ]
[Omega1*hbar                       hbar*Omega2        ]
[-----------  -delta1*hbar         -----------        ]
[     2                                 2             ]
[                                                     ]
[             Omega2*hbar                             ]
[     0       -----------   -delta1*hbar - delta2*hbar]
[                  2                                  ]



We define the density matrix.

>>> rho=define_density_matrix(Ne)
>>> pprint( rho ,use_unicode=use_unicode)
[rho11  rho12  rho13]
[                   ]
[rho21  rho22  rho23]
[                   ]
[rho31  rho32  rho33]



The hamiltonian part of the equations is
\\begin{equation}
    \\dot{\\hat{\\rho}}=\\frac{i}{\\hbar}[\\hat{\\rho}, \\hat{\\tilde{H}}]
\\end{equation}

>>> hamiltonian_terms=(I/hbar*(rho*Htilde-Htilde*rho)).expand()

There are two Lindblad operators, since there are two spontaneous decay channels.

>>> lindblad_terms =gamma[1,0]*lindblad_operator(ket(1,Ne)*bra(2,Ne),rho)
>>> lindblad_terms+=gamma[2,1]*lindblad_operator(ket(2,Ne)*bra(3,Ne),rho)
    
>>> pprint(lindblad_terms, num_columns=120,use_unicode=use_unicode)
[                          -gamma_21*rho12                    -gamma_32*rho13          ]
[ gamma_21*rho22           ----------------                   ----------------         ]
[                                 2                                  2                 ]
[                                                                                      ]
[-gamma_21*rho21                                        gamma_21*rho23   gamma_32*rho23]
[----------------  -gamma_21*rho22 + gamma_32*rho33   - -------------- - --------------]
[       2                                                     2                2       ]
[                                                                                      ]
[-gamma_32*rho31     gamma_21*rho32   gamma_32*rho32                                   ]
[----------------  - -------------- - --------------           -gamma_32*rho33         ]
[       2                  2                2                                          ]



# Optical Bloch Equations
The Optical Bloch equations are thus.

>>> eqs=hamiltonian_terms + lindblad_terms

>>> for mu in range(0,Ne**2-1 -(Ne**2 - Ne)/2+1):
...     ii,jj,s=IJ(mu,Ne)
...     i=ii-1; j=jj-1
...     pprint([Derivative(rho[i,j],t),"=",eqs[i,j]],use_unicode=use_unicode)
... 
                                                         ______ 
 d             I*Omega1*rho12                    I*rho21*Omega1 
[--(rho11), =, -------------- + gamma_21*rho22 - --------------]
 dt                  2                                 2        
<BLANKLINE>
 d               I*Omega1*rho12   I*Omega2*rho23                              
[--(rho22), =, - -------------- + -------------- - gamma_21*rho22 + gamma_32*r
 dt                    2                2                                     
<BLANKLINE>
               ______           ______ 
       I*rho21*Omega1   I*rho32*Omega2 
ho33 + -------------- - --------------]
             2                2        
                                                           ______ 
 d               I*Omega2*rho23                    I*rho32*Omega2 
[--(rho33), =, - -------------- - gamma_32*rho33 + --------------]
 dt                    2                                 2        
<BLANKLINE>
 d               I*Omega1*rho11   I*Omega1*rho22                    gamma_21*r
[--(rho21), =, - -------------- + -------------- + I*delta1*rho21 - ----------
 dt                    2                2                                 2   
<BLANKLINE>
               ______ 
ho21   I*rho31*Omega2 
---- - --------------]
             2        
 d             I*Omega1*rho32   I*Omega2*rho21                                
[--(rho31), =, -------------- - -------------- + I*delta1*rho31 + I*delta2*rho
 dt                  2                2                                       
<BLANKLINE>
     gamma_32*rho31 
31 - --------------]
           2        
<BLANKLINE>
 d               I*Omega2*rho22   I*Omega2*rho33                    gamma_21*r
[--(rho32), =, - -------------- + -------------- + I*delta2*rho32 - ----------
 dt                    2                2                                 2   
<BLANKLINE>
                                ______ 
ho32   gamma_32*rho32   I*rho31*Omega1 
---- - -------------- + --------------]
             2                2        



which is how most literature will show the equations. However, a more convenient way to express this equations is to explicitly asume a normalized and hermitian density matrix

>>> rho=define_density_matrix(Ne,explicitly_hermitian=True,normalized=True)
>>> pprint( rho ,use_unicode=use_unicode)
[                    _____  _____]
[-rho22 - rho33 + 1  rho21  rho31]
[                                ]
[                           _____]
[      rho21         rho22  rho32]
[                                ]
[      rho31         rho32  rho33]



>>> hamiltonian_terms = (I/hbar*(rho*Htilde-Htilde*rho)).expand()
>>> lindblad_terms    =gamma[1,0]*lindblad_operator(ket(1,Ne)*bra(2,Ne),rho)
>>> lindblad_terms   +=gamma[2,1]*lindblad_operator(ket(2,Ne)*bra(3,Ne),rho)
    
>>> eqs=hamiltonian_terms + lindblad_terms
>>> pprint(eqs,num_columns=120,use_unicode=use_unicode)
[                               _____                            ______                                      _____    
[                      I*Omega1*rho21                    I*rho21*Omega1                             I*Omega2*rho31    
[                      -------------- + gamma_21*rho22 - --------------                             -------------- - I
[                            2                                 2                                          2           
[                                                                                                                     
[                                                                                       ______             _____      
[                 I*Omega1*rho33   I*Omega1                    gamma_21*rho21   I*rho31*Omega2    I*Omega1*rho21   I*O
[I*Omega1*rho22 + -------------- - -------- + I*delta1*rho21 - -------------- - --------------  - -------------- + ---
[                       2             2                              2                2                 2             
[                                                                                                                     
[                                                                                                                     
[     I*Omega1*rho32   I*Omega2*rho21                                     gamma_32*rho31          I*Omega2*rho22   I*O
[     -------------- - -------------- + I*delta1*rho31 + I*delta2*rho31 - --------------        - -------------- + ---
[           2                2                                                  2                       2             
<BLANKLINE>
                         _____                            ______     ______                                           
        _____   gamma_21*rho21           ______   I*rho33*Omega1   I*Omega1                         _____            _
*delta1*rho21 - -------------- - I*rho22*Omega1 - -------------- + --------              - I*delta1*rho31 - I*delta2*r
                      2                                 2             2                                               
<BLANKLINE>
      _____                                             ______           ______             _____                     
mega2*rho32                                     I*rho21*Omega1   I*rho32*Omega2    I*Omega1*rho31            _____   g
----------- - gamma_21*rho22 + gamma_32*rho33 + -------------- - --------------  - -------------- - I*delta2*rho32 - -
   2                                                  2                2                 2                            
<BLANKLINE>
                                                                         ______                                      _
mega2*rho33                    gamma_21*rho32   gamma_32*rho32   I*rho31*Omega1                             I*Omega2*r
----------- + I*delta2*rho32 - -------------- - -------------- + --------------                           - ----------
   2                                 2                2                2                                          2   
<BLANKLINE>
                _____     ______ _____     ______ _____         ]
____   gamma_32*rho31   I*Omega1*rho32   I*Omega2*rho21         ]
ho31 - -------------- - -------------- + --------------         ]
             2                2                2                ]
                                                                ]
        _____            _____           ______           ______]
amma_21*rho32   gamma_32*rho32   I*rho22*Omega2   I*rho33*Omega2]
------------- - -------------- + -------------- - --------------]
     2                2                2                2       ]
                                                                ]
____                            ______                          ]
ho32                    I*rho32*Omega2                          ]
---- - gamma_32*rho33 + --------------                          ]
                              2                                 ]



and only consider the equations for the populations $\\rho_{ii}$ for $i>1$ and the real and imaginary parts of the coherences below the diagonal.

>>> pprint( re(eqs[1,1]) ,use_unicode=use_unicode)
                                     /       _____\     /       _____\     /  
                                   im\Omega1*rho21/   im\Omega2*rho32/   im\rh
-gamma_21*rho22 + gamma_32*rho33 + ---------------- - ---------------- - -----
                                          2                  2                
<BLANKLINE>
    ______\     /      ______\
o21*Omega1/   im\rho32*Omega2/
----------- + ----------------
  2                  2        



>>> pprint( re(eqs[1,0]) ,use_unicode=use_unicode)
<BLANKLINE>
                    gamma_21*re(rho21)                      rho33*im(Omega1)  
-delta1*im(rho21) - ------------------ - rho22*im(Omega1) - ---------------- +
                            2                                      2          
<BLANKLINE>
                /      ______\
 im(Omega1)   im\rho31*Omega2/
 ---------- + ----------------
     2               2        



>>> pprint( im(eqs[1,0]) ,use_unicode=use_unicode)
<BLANKLINE>
                   gamma_21*im(rho21)                      rho33*re(Omega1)   
delta1*re(rho21) - ------------------ + rho22*re(Omega1) + ---------------- - 
                           2                                      2           
<BLANKLINE>
               /      ______\
re(Omega1)   re\rho31*Omega2/
---------- - ----------------
    2               2        



If the density matrix is represented as a vector whose components are the these independent components of the density matrix

>>> rho_vect=define_rho_vector(rho,Ne)
>>> pprint(rho_vect,use_unicode=use_unicode)
[  rho22  ]
[         ]
[  rho33  ]
[         ]
[re(rho21)]
[         ]
[re(rho31)]
[         ]
[re(rho32)]
[         ]
[im(rho21)]
[         ]
[im(rho31)]
[         ]
[im(rho32)]



Then the equations can be re-written as linear combinations of these components plus an independent term.
\\begin{equation}
    \\dot{\\vec{\\rho}} = \\hat{A} \\vec{\\rho} + \\vec{b}
\\end{equation}
with $\\hat{A}$ a linear operator acting in this vector space and $\\vec{b}$ the vector of independent terms.

>>> A,b=calculate_A_b(eqs,rho,Ne)
>>> pprint([A,b],num_columns=120,use_unicode=use_unicode)
[[ -gamma_21      gamma_32     im(Omega1)          0              -im(Omega2)       -re(Omega1)         0             
 [                                                                                                                    
 [     0         -gamma_32         0               0              im(Omega2)             0              0             
 [                                                                                                                    
 [              -im(Omega1)    -gamma_21     -im(Omega2)                                            re(Omega2)        
 [-im(Omega1)   ------------   ----------    ------------              0              -delta1       ----------        
 [                   2             2              2                                                     2             
 [                                                                                                                    
 [                             im(Omega2)     -gamma_32          -im(Omega1)        re(Omega2)                        
 [     0             0         ----------     ----------         ------------       ----------   -delta1 - delta2     
 [                                 2              2                   2                 2                             
 [                                                                                                                    
 [ im(Omega2)   -im(Omega2)                   im(Omega1)       gamma_21   gamma_32                 -re(Omega1)        
 [ ----------   ------------       0          ----------     - -------- - --------       0         ------------       
 [     2             2                            2               2          2                          2             
 [                                                                                                                    
 [               re(Omega1)                  -re(Omega2)                            -gamma_21      -im(Omega2)        
 [ re(Omega1)    ----------      delta1      ------------              0            ----------     ------------       
 [                   2                            2                                     2               2             
 [                                                                                                                    
 [                            -re(Omega2)                         re(Omega1)        im(Omega2)      -gamma_32         
 [     0             0        ------------  delta1 + delta2       ----------        ----------      ----------        
 [                                 2                                  2                 2               2             
 [                                                                                                                    
 [-re(Omega2)    re(Omega2)                   re(Omega1)                                            im(Omega1)       g
 [------------   ----------        0          ----------            delta2               0          ----------     - -
 [     2             2                            2                                                     2             
<BLANKLINE>
  re(Omega2)      ], [     0      ]]
                  ]  [            ] 
  -re(Omega2)     ]  [     0      ] 
                  ]  [            ] 
                  ]  [-im(Omega1) ] 
       0          ]  [------------] 
                  ]  [     2      ] 
                  ]  [            ] 
 -re(Omega1)      ]  [     0      ] 
 ------------     ]  [            ] 
      2           ]  [     0      ] 
                  ]  [            ] 
                  ]  [ re(Omega1) ] 
    -delta2       ]  [ ---------- ] 
                  ]  [     2      ] 
                  ]  [            ] 
                  ]  [     0      ] 
       0          ]  [            ] 
                  ]  [     0      ] 
                  ]                 
 -im(Omega1)      ]                 
 ------------     ]                 
      2           ]                 
                  ]                 
amma_21   gamma_32]                 
------- - --------]                 
  2          2    ]                 



Explicitly, this is

>>> eqs_new=A*rho_vect - b
>>> pprint(eqs_new,num_columns=120,use_unicode=use_unicode)
[         -gamma_21*rho22 + gamma_32*rho33 - re(Omega1)*im(rho21) + re(Omega2)*im(rho32) + re(rho21)*im(Omega1) - re(r
[                                                                                                                     
[                                         -gamma_32*rho33 - re(Omega2)*im(rho32) + re(rho32)*im(Omega2)               
[                                                                                                                     
[                        gamma_21*re(rho21)                      rho33*im(Omega1)   re(Omega2)*im(rho31)   re(rho31)*i
[    -delta1*im(rho21) - ------------------ - rho22*im(Omega1) - ---------------- + -------------------- - -----------
[                                2                                      2                    2                      2 
[                                                                                                                     
[  gamma_32*re(rho31)                                  re(Omega1)*im(rho32)   re(Omega2)*im(rho21)   re(rho21)*im(Omeg
[- ------------------ + (-delta1 - delta2)*im(rho31) - -------------------- + -------------------- + -----------------
[          2                                                    2                      2                      2       
[                                                                                                                     
[                       rho22*im(Omega2)   rho33*im(Omega2)   /  gamma_21   gamma_32\             re(Omega1)*im(rho31)
[   -delta2*im(rho32) + ---------------- - ---------------- + |- -------- - --------|*re(rho32) - --------------------
[                              2                  2           \     2          2    /                      2          
[                                                                                                                     
[                       gamma_21*im(rho21)                      rho33*re(Omega1)   re(Omega1)   re(Omega2)*re(rho31)  
[    delta1*re(rho21) - ------------------ + rho22*re(Omega1) + ---------------- - ---------- - -------------------- -
[                               2                                      2               2                 2            
[                                                                                                                     
[  gamma_32*im(rho31)                                 re(Omega1)*re(rho32)   re(Omega2)*re(rho21)   im(Omega1)*im(rho3
[- ------------------ + (delta1 + delta2)*re(rho31) + -------------------- - -------------------- - ------------------
[          2                                                   2                      2                      2        
[                                                                                                                     
[                      rho22*re(Omega2)   rho33*re(Omega2)   /  gamma_21   gamma_32\             re(Omega1)*re(rho31) 
[   delta2*re(rho32) - ---------------- + ---------------- + |- -------- - --------|*im(rho32) + -------------------- 
[                             2                  2           \     2          2    /                      2           
<BLANKLINE>
ho32)*im(Omega2)          ]
                          ]
                          ]
                          ]
m(Omega2)   im(Omega1)    ]
--------- + ----------    ]
                2         ]
                          ]
a2)   re(rho32)*im(Omega1)]
--- - --------------------]
               2          ]
                          ]
   re(rho31)*im(Omega1)   ]
 + --------------------   ]
            2             ]
                          ]
 im(Omega2)*im(rho31)     ]
 --------------------     ]
          2               ]
                          ]
2)   im(Omega2)*im(rho21) ]
-- + -------------------- ]
              2           ]
                          ]
  im(Omega1)*im(rho31)    ]
+ --------------------    ]
           2              ]



Which is the same as the equations in the previous form.

>>> ss_comp={ rho[i,j]:re(rho[i,j])+I*im(rho[i,j]) for j in range(Ne) for i in range(Ne)}
    
>>> test= (eqs_new - Matrix([re(eqs[1,1]),re(eqs[2,2]),
...                           re(eqs[1,0]),re(eqs[2,0]),re(eqs[2,1]),
...                           im(eqs[1,0]),im(eqs[2,0]),im(eqs[2,1])])).expand()
    
>>> test=test.subs(ss_comp)
>>> pprint(test,use_unicode=use_unicode)
[0]
[ ]
[0]
[ ]
[0]
[ ]
[0]
[ ]
[0]
[ ]
[0]
[ ]
[0]
[ ]
[0]



[1]  H.J. Metcalf and P. van der Straten. Laser Cooling and Trapping. Graduate Texts in Contempo-
rary Physics. Springer New York, 2001.

[2] Daniel Adam Steck. Quantum and Atom Optics. Oregon Center for Optics and Department of Physics, University of Oregon Copyright © 200

[]

"""
__doc__=__doc__.replace("+IGNORE_PLOT_STEP1", "+ELLIPSIS\n[<...>]")
__doc__=__doc__.replace("+IGNORE_PLOT_STEP2", "+ELLIPSIS\n<...>")
__doc__=__doc__.replace("+IGNORE_PLOT_STEP3", "+ELLIPSIS\n(...)")
__doc__=__doc__.replace("+IGNORE_PLOT_STEP4", "\n")
