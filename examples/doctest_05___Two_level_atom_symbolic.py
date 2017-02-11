# -*- coding: utf-8 -*-
# Copyright (C) 2017 Oscar Gerardo Lazo Arjona
# mailto: oscar.lazoarjona@physics.ox.ac.uk

__doc__ = r"""

>>> from fast import *
>>> init_printing()
>>> use_unicode=True; use_unicode=False

>>> from sympy import sin,cos,exp,sqrt,pi,zeros

We define the number of states and of radiation fields.

>>> Ne=2
>>> Nl=1

We define the variables related to the laser field.

>>> E0,omega_laser=define_laser_variables(Nl)
>>> pprint(E0,use_unicode=use_unicode)
[E_0^1]



>>> pprint(omega_laser,use_unicode=use_unicode)
[omega^1]



We define a few important symbols.

>>> t,hbar,e=symbols("t hbar e",positive=True)
>>> pprint([t,hbar,e],use_unicode=use_unicode)
[t, hbar, e]



We write an electric field propagating trough the $\\hat{x}$ direction polarized in the $\\hat{z}$ direction. First the wave vector:

>>> phi=0; theta=pi/2; alpha=pi/2; beta=0
    
>>> k=Matrix([cos(phi)*sin(theta),sin(phi)*sin(theta),cos(theta)])
>>> pprint(k,use_unicode=use_unicode)
[1]
[ ]
[0]
[ ]
[0]



The polarization vectors.

>>> ep=polarization_vector(phi,theta,alpha,beta, 1)
>>> em=polarization_vector(phi,theta,alpha,beta,-1)
>>> pprint([ep,em],use_unicode=use_unicode)
[[0], [0]]
 [ ]  [ ] 
 [0]  [0] 
 [ ]  [ ] 
 [1]  [1] 



The electric field (evaluated in $\\vec{R}=0$).

>>> E_cartesian=(E0[0]/2*ep*exp(-I*omega_laser[0]*t) + E0[0].conjugate()/2*em*exp( I*omega_laser[0]*t))
>>> pprint(E_cartesian,use_unicode=use_unicode)
[                   0                    ]
[                                        ]
[                   0                    ]
[                                        ]
[       -I*omega^1*t    I*omega^1*t _____]
[E_0^1*e               e           *E_0^1]
[------------------- + ------------------]
[         2                    2         ]



We write the electric field in the helicity basis.

>>> E=cartesian_to_helicity(E_cartesian)
>>> pprint(E,use_unicode=use_unicode)
[                   0                    ]
[                                        ]
[       -I*omega^1*t    I*omega^1*t _____]
[E_0^1*e               e           *E_0^1]
[------------------- + ------------------]
[         2                    2         ]
[                                        ]
[                   0                    ]



We define the position operator.

>>> r=define_r_components(Ne,helicity=True,explicitly_hermitian=True)
>>> pprint(r,use_unicode=use_unicode)
[[    0      -r_{+1;21}], [   0      r_{0;21}], [    0      -r_{-1;21}]]
 [                     ]  [                  ]  [                     ] 
 [r_{-1;21}      0     ]  [r_{0;21}     0    ]  [r_{+1;21}      0     ] 



The frequencies of the energy levels, the resonant frequencies, and the decay frequencies.

>>> omega_level,omega,gamma=define_frequencies(Ne,explicitly_antisymmetric=True)
>>> pprint(omega_level,use_unicode=use_unicode)
[omega_1, omega_2]



>>> pprint(omega,use_unicode=use_unicode)
[   0      -omega_21]
[                   ]
[omega_21      0    ]



>>> pprint(gamma,use_unicode=use_unicode)
[   0      -gamma_21]
[                   ]
[gamma_21      0    ]



The atomic hamiltonian is

>>> H0=Matrix([[hbar*omega_level[i]*KroneckerDelta(i,j) for j in range(Ne)] for i in range(Ne)])
>>> pprint(H0,use_unicode=use_unicode)
[hbar*omega_1       0      ]
[                          ]
[     0        hbar*omega_2]



The interaction hamiltonian is

>>> H1=e*helicity_dot_product(E,r)
>>> pprint(H1,num_columns=120,use_unicode=use_unicode)
[                                                                  /       -I*omega^1*t    I*omega^1*t _____\]
[                                                                  |E_0^1*e               e           *E_0^1|]
[                          0                            e*r_{0;21}*|------------------- + ------------------|]
[                                                                  \         2                    2         /]
[                                                                                                            ]
[           /       -I*omega^1*t    I*omega^1*t _____\                                                       ]
[           |E_0^1*e               e           *E_0^1|                                                       ]
[e*r_{0;21}*|------------------- + ------------------|                            0                          ]
[           \         2                    2         /                                                       ]



and the complete hamiltonian is

>>> H=H0+H1
>>> pprint(H,num_columns=120,use_unicode=use_unicode)
[                                                                  /       -I*omega^1*t    I*omega^1*t _____\]
[                                                                  |E_0^1*e               e           *E_0^1|]
[                    hbar*omega_1                       e*r_{0;21}*|------------------- + ------------------|]
[                                                                  \         2                    2         /]
[                                                                                                            ]
[           /       -I*omega^1*t    I*omega^1*t _____\                                                       ]
[           |E_0^1*e               e           *E_0^1|                                                       ]
[e*r_{0;21}*|------------------- + ------------------|                      hbar*omega_2                     ]
[           \         2                    2         /                                                       ]



# Rotating wave approximation
Notice that the electric field can be separated by terms with positive and negative frequency:

>>> E_cartesian_p=E0[0]            /2*ep*exp(-I*omega_laser[0]*t)
>>> E_cartesian_m=E0[0].conjugate()/2*em*exp( I*omega_laser[0]*t)
    
>>> E_p=cartesian_to_helicity(E_cartesian_p)
>>> E_m=cartesian_to_helicity(E_cartesian_m)
    
>>> pprint([E_p,E_m],use_unicode=use_unicode)
[[         0         ], [        0         ]]
 [                   ]  [                  ] 
 [       -I*omega^1*t]  [ I*omega^1*t _____] 
 [E_0^1*e            ]  [e           *E_0^1] 
 [-------------------]  [------------------] 
 [         2         ]  [        2         ] 
 [                   ]  [                  ] 
 [         0         ]  [        0         ] 



>>> pprint( simplify(E-(E_p+E_m)) ,use_unicode=use_unicode)
[0]
[ ]
[0]
[ ]
[0]



The position operator can also be separated in this way. We go to the interaction picture (with $\\hat{H}_0$ as the undisturbed hamiltonian)

>>> r_I=[ Matrix([[exp(I*omega[i,j]*t)*r[p][i,j] for j in range(Ne)] for i in range(Ne)]) for p in range(3)]
>>> pprint(r_I[0],use_unicode=use_unicode)
[                                     -I*omega_21*t]
[           0             -r_{+1;21}*e             ]
[                                                  ]
[           I*omega_21*t                           ]
[r_{-1;21}*e                          0            ]



>>> pprint(r_I[1],use_unicode=use_unicode)
[                                  -I*omega_21*t]
[          0             r_{0;21}*e             ]
[                                               ]
[          I*omega_21*t                         ]
[r_{0;21}*e                         0           ]



>>> pprint(r_I[2],use_unicode=use_unicode)
[                                     -I*omega_21*t]
[           0             -r_{-1;21}*e             ]
[                                                  ]
[           I*omega_21*t                           ]
[r_{+1;21}*e                          0            ]



Which can be decomposed as

>>> r_I_p=[ Matrix([[ delta_greater(j,i)*exp(-I*omega[j,i]*t)*r[p][i,j] for j in range(Ne)]for i in range(Ne)]) for p in range(3)]
>>> pprint(r_I_p[0],use_unicode=use_unicode)
[               -I*omega_21*t]
[0  -r_{+1;21}*e             ]
[                            ]
[0              0            ]



>>> pprint(r_I_p[1],use_unicode=use_unicode)
[             -I*omega_21*t]
[0  r_{0;21}*e             ]
[                          ]
[0             0           ]



>>> pprint(r_I_p[2],use_unicode=use_unicode)
[               -I*omega_21*t]
[0  -r_{-1;21}*e             ]
[                            ]
[0              0            ]



>>> r_I_m=[ Matrix([[ delta_lesser( j,i)*exp( I*omega[i,j]*t)*r[p][i,j] for j in range(Ne)]for i in range(Ne)]) for p in range(3)]
>>> pprint(r_I_m[0],use_unicode=use_unicode)
[           0             0]
[                          ]
[           I*omega_21*t   ]
[r_{-1;21}*e              0]



>>> pprint(r_I_m[1],use_unicode=use_unicode)
[          0             0]
[                         ]
[          I*omega_21*t   ]
[r_{0;21}*e              0]



>>> pprint(r_I_m[2],use_unicode=use_unicode)
[           0             0]
[                          ]
[           I*omega_21*t   ]
[r_{+1;21}*e              0]



that summed equal $\\vec{\\hat{r}}_I$

>>> pprint( [r_I[p]-(r_I_p[p]+r_I_m[p]) for p in range(3)] ,use_unicode=use_unicode)
[[0  0], [0  0], [0  0]]
 [    ]  [    ]  [    ] 
 [0  0]  [0  0]  [0  0] 



Thus the interaction hamiltonian in the interaciton picture is
\\begin{equation}
    \\hat{H}_{1I}=e\\vec{E}\\cdot \\vec{\\hat{r}}_I= e(\\vec{E}^{(+)}\\cdot \\vec{\\hat{r}}^{(+)}_I + \\vec{E}^{(+)}\\cdot \\vec{\\hat{r}}^{(-)}_I + \\vec{E}^{(-)}\\cdot \\vec{\\hat{r}}^{(+)}_I + \\vec{E}^{(-)}\\cdot \\vec{\\hat{r}}^{(-)}_I)
\\end{equation}

>>> H1I=e*helicity_dot_product(E,r_I)
>>> pprint(H1I,num_columns=120,use_unicode=use_unicode)
[                                                                                /       -I*omega^1*t    I*omega^1*t _
[                                                                                |E_0^1*e               e           *E
[                                 0                                   e*r_{0;21}*|------------------- + --------------
[                                                                                \         2                    2     
[                                                                                                                     
[           /       -I*omega^1*t    I*omega^1*t _____\                                                                
[           |E_0^1*e               e           *E_0^1|  I*omega_21*t                                                  
[e*r_{0;21}*|------------------- + ------------------|*e                                               0              
[           \         2                    2         /                                                                
<BLANKLINE>
____\               ]
_0^1|  -I*omega_21*t]
----|*e             ]
    /               ]
                    ]
                    ]
                    ]
                    ]
                    ]



Since both $\\omega^l$ and $\\omega_{ij}$ are in the order of THz, the terms that have frequencies with the same sign are summed, and thus also of the order of THz. The frequencies in the terms with oposite signs however, are detunings of the order of MHz. Since we are only interested in the coarse-grained evolution of the density matrix, we may omit the fast terms and approximate

\\begin{equation}
    \\hat{H}_{1I} \\simeq \\hat{H}_{1I,RWA}= e( \\vec{E}^{(+)}\\cdot \\vec{\\hat{r}}^{(-)}_I + \\vec{E}^{(-)}\\cdot \\vec{\\hat{r}}^{(+)}_I )
\\end{equation}

That is known as the rotating wave approximation (RWA).

>>> H1IRWA=e*(helicity_dot_product(E_p,r_I_m)+helicity_dot_product(E_m,r_I_p))
>>> pprint(H1IRWA,use_unicode=use_unicode)
[                                                          I*omega^1*t  -I*ome
[                                              e*r_{0;21}*e           *e      
[                     0                        -------------------------------
[                                                                   2         
[                                                                             
[                  -I*omega^1*t  I*omega_21*t                                 
[E_0^1*e*r_{0;21}*e            *e                                             
[--------------------------------------------                       0         
[                     2                                                       
<BLANKLINE>
ga_21*t _____]
       *E_0^1]
-------------]
             ]
             ]
             ]
             ]
             ]
             ]



 Returning to the Schrödinger picture we have.

>>> r_p=[ Matrix([[ delta_greater(j,i)*r[p][i,j] for j in range(Ne)]for i in range(Ne)]) for p in range(3)]
>>> pprint(r_p,use_unicode=use_unicode)
[[0  -r_{+1;21}], [0  r_{0;21}], [0  -r_{-1;21}]]
 [             ]  [           ]  [             ] 
 [0      0     ]  [0     0    ]  [0      0     ] 



>>> r_m=[ Matrix([[ delta_lesser( j,i)*r[p][i,j] for j in range(Ne)]for i in range(Ne)]) for p in range(3)]
>>> pprint(r_m,use_unicode=use_unicode)
[[    0      0], [   0      0], [    0      0]]
 [            ]  [           ]  [            ] 
 [r_{-1;21}  0]  [r_{0;21}  0]  [r_{+1;21}  0] 



>>> pprint( [r[p]-(r_p[p]+r_m[p]) for p in range(3)] ,use_unicode=use_unicode)
[[0  0], [0  0], [0  0]]
 [    ]  [    ]  [    ] 
 [0  0]  [0  0]  [0  0] 



Thus the interaction hamiltonian in the Schrödinger picture in the rotating wave approximation is

>>> H1RWA=e*(helicity_dot_product(E_p,r_m)+helicity_dot_product(E_m,r_p))
>>> pprint(H1RWA,use_unicode=use_unicode)
[                                            I*omega^1*t _____]
[                                e*r_{0;21}*e           *E_0^1]
[              0                 -----------------------------]
[                                              2              ]
[                                                             ]
[                  -I*omega^1*t                               ]
[E_0^1*e*r_{0;21}*e                                           ]
[------------------------------                0              ]
[              2                                              ]



And the complete hamiltonian in the Schrödinger picture in the rotating wave approximation is

>>> HRWA=H0+H1RWA
>>> pprint(HRWA,use_unicode=use_unicode)
[                                            I*omega^1*t _____]
[                                e*r_{0;21}*e           *E_0^1]
[         hbar*omega_1           -----------------------------]
[                                              2              ]
[                                                             ]
[                  -I*omega^1*t                               ]
[E_0^1*e*r_{0;21}*e                                           ]
[------------------------------          hbar*omega_2         ]
[              2                                              ]



# Rotating Frame
Next we will make a phase transformation in order to eliminate the explicit time dependance of the equations.

>>> c,ctilde,phase=define_psi_coefficients(Ne)
>>> pprint([c,ctilde,phase],use_unicode=use_unicode)
[[c1(t)], [\tilde{c}_{1}(t)], [theta1]]
 [     ]  [                ]  [      ] 
 [c2(t)]  [\tilde{c}_{2}(t)]  [theta2] 



>>> psi=Matrix([ exp(I*phase[i]*t)*ctilde[i] for i in range(Ne)])
>>> pprint(psi,use_unicode=use_unicode)
[                  I*t*theta1]
[\tilde{c}_{1}(t)*e          ]
[                            ]
[                  I*t*theta2]
[\tilde{c}_{2}(t)*e          ]



The Schrödinger equation $i\\hbar \\partial_t |\\psi\\rangle=\\hat{H}_{RWA}$ is

>>> lhs=Matrix([(I*hbar*Derivative(psi[i],t).doit()).expand() for i in range(Ne)])
>>> pprint(lhs,use_unicode=use_unicode)
[                               I*t*theta1]
[-hbar*theta1*\tilde{c}_{1}(t)*e          ]
[                                         ]
[                               I*t*theta2]
[-hbar*theta2*\tilde{c}_{2}(t)*e          ]



>>> rhs=HRWA*psi
>>> pprint(rhs,num_columns=120,use_unicode=use_unicode)
[                             I*omega^1*t  I*t*theta2 _____                                             ]
[e*r_{0;21}*\tilde{c}_{2}(t)*e           *e          *E_0^1                                  I*t*theta1 ]
[---------------------------------------------------------- + hbar*omega_1*\tilde{c}_{1}(t)*e           ]
[                            2                                                                          ]
[                                                                                                       ]
[                                   -I*omega^1*t  I*t*theta1                                            ]
[E_0^1*e*r_{0;21}*\tilde{c}_{1}(t)*e            *e                                            I*t*theta2]
[----------------------------------------------------------- + hbar*omega_2*\tilde{c}_{2}(t)*e          ]
[                             2                                                                         ]



We multiply each of these equations by $e^{-i \\theta_i t}$ and substracting $i \\theta_i \\tilde{c}_i$

>>> lhs_new=Matrix([simplify(  lhs[i]*exp(-I*phase[i]*t) +hbar*phase[i]*ctilde[i] ) for i in range(Ne)])
>>> pprint(lhs_new,use_unicode=use_unicode)
[0]
[ ]
[0]



>>> rhs_new=Matrix([simplify(  rhs[i]*exp(-I*phase[i]*t) +hbar*phase[i]*ctilde[i] ) for i in range(Ne)])
>>> pprint(rhs_new,num_columns=120,use_unicode=use_unicode)
[                             I*omega^1*t  -I*t*theta1  I*t*theta2 _____                                              
[e*r_{0;21}*\tilde{c}_{2}(t)*e           *e           *e          *E_0^1                                              
[----------------------------------------------------------------------- + hbar*omega_1*\tilde{c}_{1}(t) + hbar*theta1
[                                   2                                                                                 
[                                                                                                                     
[                                   -I*omega^1*t  I*t*theta1  -I*t*theta2                                             
[E_0^1*e*r_{0;21}*\tilde{c}_{1}(t)*e            *e          *e                                                        
[------------------------------------------------------------------------ + hbar*omega_2*\tilde{c}_{2}(t) + hbar*theta
[                                   2                                                                                 
<BLANKLINE>
                  ]
                  ]
*\tilde{c}_{1}(t) ]
                  ]
                  ]
                  ]
                  ]
2*\tilde{c}_{2}(t)]
                  ]



It can be seen that the equations loose their explicit time dependance only if $\\omega^{1} - \\theta_{1} + \\theta_{2}=0$. Which is satisfied if

>>> phase_transformation=solve(omega_laser[0]+phase[1]-phase[0],phase[1],dict=True)[0]
>>> pprint(phase_transformation,use_unicode=use_unicode)
{theta2: -omega^1 + theta1}



There is a free parameter $\\theta_1$, which is to be expected, since state vetors $|\\psi\\rangle$ always have a global phase invariance

>>> pprint(psi.subs(phase_transformation),use_unicode=use_unicode)
[                        I*t*theta1       ]
[      \tilde{c}_{1}(t)*e                 ]
[                                         ]
[                  I*t*(-omega^1 + theta1)]
[\tilde{c}_{2}(t)*e                       ]



Thus the equations become

>>> pprint(lhs_new,use_unicode=use_unicode)
[0]
[ ]
[0]



>>> rhs_new=simplify(rhs_new.subs(phase_transformation))
>>> pprint(rhs_new,use_unicode=use_unicode)
[                                 _____                                       
[     e*r_{0;21}*\tilde{c}_{2}(t)*E_0^1                                       
[     --------------------------------- + hbar*(omega_1 + theta1)*\tilde{c}_{1
[                     2                                                       
[                                                                             
[E_0^1*e*r_{0;21}*\tilde{c}_{1}(t)                                            
[--------------------------------- + hbar*(-omega^1 + omega_2 + theta1)*\tilde
[                2                                                            
<BLANKLINE>
          ]
          ]
}(t)      ]
          ]
          ]
          ]
{c}_{2}(t)]
          ]



It can be seen that this is the Schrödinger equation derived from an effective hamiltonian $\\tilde{H}$

>>> Htilde=Matrix([ [Derivative(rhs_new[i],ctilde[j]).doit() for j in range(Ne)] for i in range(Ne)])
>>> pprint(Htilde,use_unicode=use_unicode)
[                                             _____         ]
[                                  e*r_{0;21}*E_0^1         ]
[hbar*(omega_1 + theta1)           ----------------         ]
[                                         2                 ]
[                                                           ]
[   E_0^1*e*r_{0;21}                                        ]
[   ----------------      hbar*(-omega^1 + omega_2 + theta1)]
[          2                                                ]



We can see that it is convenient to choose $\\theta_1=-\\omega_1$ to simplify the hamiltonian. Also, we can recognize $\\omega^1-\\omega_2+\\omega_1=\\delta$ as the detuning of the laser field relative to the atomic transition $\\omega_{21}=\\omega_2-\\omega_1$.

>>> delta=Symbol("delta",real=True)
>>> Htilde=Htilde.subs({phase[0]:-omega_level[0]}).subs({omega_laser[0]:delta+omega_level[1]-omega_level[0]})
>>> pprint(Htilde,use_unicode=use_unicode)
[                             _____]
[                  e*r_{0;21}*E_0^1]
[       0          ----------------]
[                         2        ]
[                                  ]
[E_0^1*e*r_{0;21}                  ]
[----------------    -delta*hbar   ]
[       2                          ]



If we define the Rabi frequency $\\Omega =e E_0^1 r_{0;21}/\\hbar$

>>> Omega=Symbol("Omega",real=False)
>>> Htilde=Htilde.subs({E0[0]:Omega*hbar/r[1][1,0]/e})
>>> pprint(Htilde,use_unicode=use_unicode)
[                 _____ ]
[            hbar*Omega ]
[    0       ---------- ]
[                2      ]
[                       ]
[Omega*hbar             ]
[----------  -delta*hbar]
[    2                  ]



We define the density matrix.

>>> rho=define_density_matrix(Ne)
>>> pprint( rho ,use_unicode=use_unicode)
[rho11  rho12]
[            ]
[rho21  rho22]



The hamiltonian part of the equations is
\\begin{equation}
    \\dot{\\hat{\\rho}}=\\frac{i}{\\hbar}[\\hat{\\rho}, \\hat{\\tilde{H}}]
\\end{equation}

>>> hamiltonian_terms=(I/hbar*(rho*Htilde-Htilde*rho)).expand()
>>> pprint(hamiltonian_terms,use_unicode=use_unicode)
[                                 _____                                    ___
[         I*Omega*rho12   I*rho21*Omega                            I*rho11*Ome
[         ------------- - -------------           -I*delta*rho12 + -----------
[               2               2                                        2    
[                                                                             
[                                                                             
[  I*Omega*rho11   I*Omega*rho22                           I*Omega*rho12   I*r
[- ------------- + ------------- + I*delta*rho21         - ------------- + ---
[        2               2                                       2            
<BLANKLINE>
__           _____]
ga   I*rho22*Omega]
-- - -------------]
           2      ]
                  ]
     _____        ]
ho21*Omega        ]
----------        ]
   2              ]



There is only one Lindblad operator, since there is only one spontaneous decay channel.

>>> lindblad_terms=gamma[1,0]*lindblad_operator(ket(1,Ne)*bra(2,Ne),rho)
>>> pprint(lindblad_terms, num_columns=120,use_unicode=use_unicode)
[                  -gamma_21*rho12 ]
[ gamma_21*rho22   ----------------]
[                         2        ]
[                                  ]
[-gamma_21*rho21                   ]
[----------------  -gamma_21*rho22 ]
[       2                          ]



# Optical Bloch Equations
The Optical Bloch equations are thus.

>>> eqs=hamiltonian_terms + lindblad_terms
>>> pprint(eqs,num_columns=120,use_unicode=use_unicode)
[                                                  _____                                                     _____    
[         I*Omega*rho12                    I*rho21*Omega                            gamma_21*rho12   I*rho11*Omega   I
[         ------------- + gamma_21*rho22 - -------------           -I*delta*rho12 - -------------- + ------------- - -
[               2                                2                                        2                2          
[                                                                                                                     
[                                                                                                                    _
[  I*Omega*rho11   I*Omega*rho22                   gamma_21*rho21           I*Omega*rho12                    I*rho21*O
[- ------------- + ------------- + I*delta*rho21 - --------------         - ------------- - gamma_21*rho22 + ---------
[        2               2                               2                        2                                2  
<BLANKLINE>
       _____]
*rho22*Omega]
------------]
     2      ]
            ]
____        ]
mega        ]
----        ]
            ]



which is how most literature will show the equations. However, a more convenient way to express this equations is to explicitly asume a normalized and hermitian density matrix

>>> rho=define_density_matrix(Ne,explicitly_hermitian=True,normalized=True)
>>> pprint( rho ,use_unicode=use_unicode)
[            _____]
[-rho22 + 1  rho21]
[                 ]
[  rho21     rho22]



>>> hamiltonian_terms = (I/hbar*(rho*Htilde-Htilde*rho)).expand()
>>> lindblad_terms    =gamma[1,0]*lindblad_operator(ket(1,Ne)*bra(2,Ne),rho)
>>> eqs=hamiltonian_terms + lindblad_terms
>>> pprint(eqs,num_columns=120,use_unicode=use_unicode)
[             _____                            _____                                  _____                     _____]
[     I*Omega*rho21                    I*rho21*Omega                 _____   gamma_21*rho21           _____   I*Omega]
[     ------------- + gamma_21*rho22 - -------------       - I*delta*rho21 - -------------- - I*rho22*Omega + -------]
[           2                                2                                     2                             2   ]
[                                                                                                                    ]
[                                                                         _____                            _____     ]
[                I*Omega                   gamma_21*rho21         I*Omega*rho21                    I*rho21*Omega     ]
[I*Omega*rho22 - ------- + I*delta*rho21 - --------------       - ------------- - gamma_21*rho22 + -------------     ]
[                   2                            2                      2                                2           ]



and only consider the equations for the populations $\\rho_{ii}$ for $i>1$ and the real and imaginary parts of the coherences below the diagonal.

>>> ss_comp={ rho[i,j]:re(rho[i,j])+I*im(rho[i,j]) for j in range(Ne) for i in range(Ne)}
>>> pprint( re(eqs[1,1].subs(ss_comp)) ,use_unicode=use_unicode)
-gamma_21*rho22 - re(Omega)*im(rho21) + re(rho21)*im(Omega)



>>> pprint( re(eqs[1,0].subs(ss_comp)) ,use_unicode=use_unicode)
                   gamma_21*re(rho21)                     im(Omega)
-delta*im(rho21) - ------------------ - rho22*im(Omega) + ---------
                           2                                  2    



>>> pprint( im(eqs[1,0].subs(ss_comp)) ,use_unicode=use_unicode)
                  gamma_21*im(rho21)                     re(Omega)
delta*re(rho21) - ------------------ + rho22*re(Omega) - ---------
                          2                                  2    



If the density matrix is represented as a vector whose components are the these independent components of the density matrix

>>> rho_vect=define_rho_vector(rho,Ne)
>>> pprint(rho_vect,use_unicode=use_unicode)
[  rho22  ]
[         ]
[re(rho21)]
[         ]
[im(rho21)]



Then the equations can be re-written as linear combinations of these components plus an independent term.
\\begin{equation}
    \\dot{\\vec{\\rho}} = \\hat{A} \\vec{\\rho} + \\vec{b}
\\end{equation}
with $\\hat{A}$ a linear operator acting in this vector space and $\\vec{b}$ the vector of independent terms.

>>> A,b=calculate_A_b(eqs,rho,Ne)
>>> pprint([A,b],use_unicode=use_unicode)
[[-gamma_21   im(Omega)   -re(Omega)], [     0     ]]
 [                                  ]  [           ] 
 [            -gamma_21             ]  [-im(Omega) ] 
 [-im(Omega)  ----------    -delta  ]  [-----------] 
 [                2                 ]  [     2     ] 
 [                                  ]  [           ] 
 [                        -gamma_21 ]  [ re(Omega) ] 
 [re(Omega)     delta     ----------]  [ --------- ] 
 [                            2     ]  [     2     ] 



Explicitly, this is

>>> eqs_new=A*rho_vect - b
>>> pprint(eqs_new,use_unicode=use_unicode)
[    -gamma_21*rho22 - re(Omega)*im(rho21) + re(rho21)*im(Omega)    ]
[                                                                   ]
[                   gamma_21*re(rho21)                     im(Omega)]
[-delta*im(rho21) - ------------------ - rho22*im(Omega) + ---------]
[                           2                                  2    ]
[                                                                   ]
[                  gamma_21*im(rho21)                     re(Omega) ]
[delta*re(rho21) - ------------------ + rho22*re(Omega) - --------- ]
[                          2                                  2     ]



Which is the same as the equations in the previous form.

>>> pprint( eqs_new - Matrix([re(eqs[1,1]),re(eqs[1,0]),im(eqs[1,0])]).subs(ss_comp) ,use_unicode=use_unicode)
[0]
[ ]
[0]
[ ]
[0]



The steady state solution of this equations is

>>> sol=solve(list(eqs_new),list(rho_vect))
>>> for mu in range(3):
...     pprint( {rho_vect[mu]:sol[rho_vect[mu]]} ,num_columns=120,use_unicode=use_unicode)
... 
                       2            2                      
                     re (Omega) + im (Omega)               
{rho22: --------------------------------------------------}
               2           2       2              2        
        4*delta  + gamma_21  + 2*re (Omega) + 2*im (Omega) 
                  2*delta*re(Omega) + gamma_21*im(Omega)       
{re(rho21): --------------------------------------------------}
                   2           2       2              2        
            4*delta  + gamma_21  + 2*re (Omega) + 2*im (Omega) 
                  2*delta*im(Omega) - gamma_21*re(Omega)       
{im(rho21): --------------------------------------------------}
                   2           2       2              2        
            4*delta  + gamma_21  + 2*re (Omega) + 2*im (Omega) 



According to literature [1], the solution should be

>>> s0=2*(re(Omega)**2+im(Omega)**2)/gamma[1,0]**2
    
>>> s=s0/(1+(2*delta/gamma[1,0])**2)
    
    
>>> rho21=-I*Omega/(2*(gamma[1,0]/2-I*delta)*(1+s))
    
>>> rerho22=( s/(2*(1+s)) ).simplify()
>>> rerho21=re(rho21).simplify()
>>> imrho21=im(rho21).simplify()
    
>>> test=[ sol[rho[1,1]]-rerho22, sol[re(rho[1,0])]-rerho21, sol[im(rho[1,0])]-imrho21 ]
    
>>> pprint( [testi.subs({Omega:re(Omega)+I*im(Omega)}).factor() for testi in test] ,use_unicode=use_unicode)
[0, 0, 0]



So our development produces the same results as the literature.

The saturation intensity is defined as the intensity needed to accumulate $\\frac{1}{4}$ of the population in the excited state when the field is in resonance ($\\delta=0$).

>>> saturation_eq=sol[rho[1,1]].subs({delta:0})-1/Integer(4)
>>> pprint( saturation_eq ,use_unicode=use_unicode)
          2            2                   
        re (Omega) + im (Omega)           1
--------------------------------------- - -
        2       2              2          4
gamma_21  + 2*re (Omega) + 2*im (Omega)    



>>> Omega_amp,alpha=symbols("\Omega_a alpha",real=True)
>>> ss={Omega:Omega_amp*cos(alpha)+I*Omega_amp*sin(alpha)}
>>> saturation_eq= saturation_eq.subs(ss).factor().simplify()
>>> pprint(saturation_eq,use_unicode=use_unicode)
            2           2  
  2*\Omega_a  - gamma_21   
---------------------------
  /          2           2\
4*\2*\Omega_a  + gamma_21 /



>>> Omega_sat=solve( saturation_eq ,Omega_amp)[1]
>>> pprint(Omega_sat,use_unicode=use_unicode)
  ___         
\/ 2 *gamma_21
--------------
      2       



Since $\\Omega =e E_0^1 r_{0;21}/\\hbar$ it follows that

>>> E0_sat=Omega_sat*hbar/e/r[1][1,0]
>>> pprint(E0_sat,use_unicode=use_unicode)
  ___              
\/ 2 *gamma_21*hbar
-------------------
    2*e*r_{0;21}   



The full width at half maximum of $\\rho_{22}$ is

>>> hm1,hm2=solve(sol[rho[1,1]]-sol[rho[1,1]].subs({delta:0})/2,delta)
>>> FWHM=hm2-hm1
>>> FWHM=FWHM.subs(ss).simplify()
>>> pprint(FWHM,use_unicode=use_unicode)
   _________________________
  /           2           2 
\/  2*\Omega_a  + gamma_21  



[1]  H.J. Metcalf and P. van der Straten. Laser Cooling and Trapping. Graduate Texts in Contempo-
rary Physics. Springer New York, 2001.

[]

"""
__doc__=__doc__.replace("+IGNORE_PLOT_STEP1", "+ELLIPSIS\n[<...>]")
__doc__=__doc__.replace("+IGNORE_PLOT_STEP2", "+ELLIPSIS\n<...>")
__doc__=__doc__.replace("+IGNORE_PLOT_STEP3", "+ELLIPSIS\n(...)")
__doc__=__doc__.replace("+IGNORE_PLOT_STEP4", "\n")
