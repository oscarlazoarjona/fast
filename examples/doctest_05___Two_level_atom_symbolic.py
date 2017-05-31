# -*- coding: utf-8 -*-
# Copyright (C) 2017 Oscar Gerardo Lazo Arjona
# mailto: oscar.lazoarjona@physics.ox.ac.uk

__doc__ = r"""

>>> from fast import *
>>> from matplotlib import pyplot
>>> from sympy import sin,cos,exp,sqrt,pi,zeros,I

>>> init_printing()
>>> print_ascii=True
>>> #print_ascii=False

>>> path="folder_05___Two_level_atom_symbolic/" 
>>> name='suite'

We will be deriving the optical Bloch equations for a two level system as that in the figure.

>>> fig=pyplot.figure(); ax=fig.add_subplot(111,aspect="equal")
    
>>> p1=[0.5,1]; p2=[1.5,3]
>>> draw_state(ax,p1,text=r"$|1\rangle$",l=1.0,alignment='right',label_displacement=0.05,fontsize=25,linewidth=4.0)
>>> draw_state(ax,p2,text=r"$|2\rangle$",l=1.0,alignment='right',label_displacement=0.05,fontsize=25,linewidth=4.0)
    
>>> excitation(ax,[p1[0]+0.25,p1[1]],[p2[0]+0.25,p2[1]], fc="red", ec="red",width=0.01, head_width=0.1, head_length=0.1) # doctest: +IGNORE_PLOT_STEP4
>>> decay(     ax,[p1[0]-0.25,p1[1]],[p2[0]-0.25,p2[1]], 0.05,10.0,color="red",linewidth=1.0) # doctest: +IGNORE_PLOT_STEP4
    
>>> pyplot.axis('off') # doctest: +IGNORE_PLOT_STEP3
>>> pyplot.savefig(path+name+'_diagram.png',bbox_inches="tight") # doctest: +IGNORE_PLOT_STEP4
<matplotlib.figure.Figure at 0x7fe853346c50>



We define the number of states and of radiation fields.

>>> Ne=2
>>> Nl=1

We define the variables related to the laser field.

>>> E0,omega_laser=define_laser_variables(Nl)
>>> fprint(E0,print_ascii=print_ascii)
[E_0^1]



>>> fprint(omega_laser,print_ascii=print_ascii)
[varpi_1]



We define a few important symbols.

>>> t,hbar,e=symbols("t hbar e",positive=True)
>>> fprint([t,hbar,e],print_ascii=print_ascii)
[t, hbar, e]



We write an electric field propagating trough the $\\hat{x}$ direction polarized in the $\\hat{z}$ direction. First the wave vector:

>>> phi=0; theta=pi/2; alpha=pi/2; beta=0
    
>>> k=Matrix([cos(phi)*sin(theta),sin(phi)*sin(theta),cos(theta)])
>>> fprint(k,print_ascii=print_ascii)
[1]
[ ]
[0]
[ ]
[0]



The polarization vectors.

>>> ep=polarization_vector(phi,theta,alpha,beta, 1)
>>> em=polarization_vector(phi,theta,alpha,beta,-1)
>>> fprint([ep,em],print_ascii=print_ascii)
[[0], [0]]
 [ ]  [ ] 
 [0]  [0] 
 [ ]  [ ] 
 [1]  [1] 



The electric field (evaluated in $\\vec{R}=0$).

>>> E_cartesian=(E0[0]/2*ep*exp(-I*omega_laser[0]*t) + E0[0].conjugate()/2*em*exp( I*omega_laser[0]*t))
>>> fprint(E_cartesian,print_ascii=print_ascii)
[                   0                    ]
[                                        ]
[                   0                    ]
[                                        ]
[       -I*t*varpi_1    I*t*varpi_1 _____]
[E_0^1*e               e           *E_0^1]
[------------------- + ------------------]
[         2                    2         ]



We draw this electric field.

>>> l1=PlaneWave(phi,theta,0,0,color="red")
>>> laseres=[l1]
>>> Nl=len(laseres)
    
>>> fig = pyplot.figure(); ax = fig.gca(projection='3d')
>>> draw_lasers_3d(ax,laseres,path+'lasers.png') # doctest: +IGNORE_PLOT_STEP4
<matplotlib.figure.Figure at 0x7fe8524d54d0>



We write the electric field in the helicity basis (see notebook "Vectors in the helicity basis and the electric field").

>>> E=cartesian_to_helicity(E_cartesian)
>>> fprint(E,print_ascii=print_ascii)
[                   0                    ]
[                                        ]
[       -I*t*varpi_1    I*t*varpi_1 _____]
[E_0^1*e               e           *E_0^1]
[------------------- + ------------------]
[         2                    2         ]
[                                        ]
[                   0                    ]



We define the position operator in the helicity basis.

>>> r=define_r_components(Ne,helicity=True,explicitly_hermitian=True)
>>> fprint(r,print_ascii=print_ascii)
[[    0      -r_{+1;21}], [   0      r_{0;21}], [    0      -r_{-1;21}]]
 [                     ]  [                  ]  [                     ] 
 [r_{-1;21}      0     ]  [r_{0;21}     0    ]  [r_{+1;21}      0     ] 



The frequencies of the energy levels, the resonant frequencies, and the decay frequencies.

>>> omega_level,omega,gamma=define_frequencies(Ne,explicitly_antisymmetric=True)
>>> fprint(omega_level,print_ascii=print_ascii)
[omega_1, omega_2]



>>> fprint(omega,print_ascii=print_ascii)
[   0      -omega_21]
[                   ]
[omega_21      0    ]



>>> fprint(gamma,print_ascii=print_ascii)
[   0      -gamma_21]
[                   ]
[gamma_21      0    ]



The atomic hamiltonian is

>>> H0=Matrix([[hbar*omega_level[i]*KroneckerDelta(i,j) for j in range(Ne)] for i in range(Ne)])
>>> fprint(H0,print_ascii=print_ascii)
[hbar*omega_1       0      ]
[                          ]
[     0        hbar*omega_2]



The interaction hamiltonian is

>>> H1=e*helicity_dot_product(E,r)
>>> fprint(H1,print_ascii=print_ascii)
[                                                                  /       -I*t*varpi_1    I*t*varpi_1 _____\]
[                                                                  |E_0^1*e               e           *E_0^1|]
[                          0                            e*r_{0;21}*|------------------- + ------------------|]
[                                                                  \         2                    2         /]
[                                                                                                            ]
[           /       -I*t*varpi_1    I*t*varpi_1 _____\                                                       ]
[           |E_0^1*e               e           *E_0^1|                                                       ]
[e*r_{0;21}*|------------------- + ------------------|                            0                          ]
[           \         2                    2         /                                                       ]



and the complete hamiltonian is

>>> H=H0+H1
>>> fprint(H,print_ascii=print_ascii)
[                                                                  /       -I*t*varpi_1    I*t*varpi_1 _____\]
[                                                                  |E_0^1*e               e           *E_0^1|]
[                    hbar*omega_1                       e*r_{0;21}*|------------------- + ------------------|]
[                                                                  \         2                    2         /]
[                                                                                                            ]
[           /       -I*t*varpi_1    I*t*varpi_1 _____\                                                       ]
[           |E_0^1*e               e           *E_0^1|                                                       ]
[e*r_{0;21}*|------------------- + ------------------|                      hbar*omega_2                     ]
[           \         2                    2         /                                                       ]



# Rotating wave approximation
Notice that the electric field can be separated by terms with positive and negative frequency:

>>> E_cartesian_p=E0[0]            /2*ep*exp(-I*omega_laser[0]*t)
>>> E_cartesian_m=E0[0].conjugate()/2*em*exp( I*omega_laser[0]*t)
    
>>> E_p=cartesian_to_helicity(E_cartesian_p)
>>> E_m=cartesian_to_helicity(E_cartesian_m)
    
>>> fprint([E_p,E_m],print_ascii=print_ascii)
[[         0         ], [        0         ]]
 [                   ]  [                  ] 
 [       -I*t*varpi_1]  [ I*t*varpi_1 _____] 
 [E_0^1*e            ]  [e           *E_0^1] 
 [-------------------]  [------------------] 
 [         2         ]  [        2         ] 
 [                   ]  [                  ] 
 [         0         ]  [        0         ] 



We check that this decomposition actually equalls the field

>>> fprint( simplify(E-(E_p+E_m)) ,print_ascii=print_ascii)
[0]
[ ]
[0]
[ ]
[0]



The position operator can also be separated in this way. We go to the interaction picture (with $\\hat{H}_0$ as the undisturbed hamiltonian)

>>> r_I=[ Matrix([[exp(I*omega[i,j]*t)*r[p][i,j] for j in range(Ne)] for i in range(Ne)]) for p in range(3)]
>>> fprint(r_I[0],print_ascii=print_ascii)
[                                     -I*omega_21*t]
[           0             -r_{+1;21}*e             ]
[                                                  ]
[           I*omega_21*t                           ]
[r_{-1;21}*e                          0            ]



>>> fprint(r_I[1],print_ascii=print_ascii)
[                                  -I*omega_21*t]
[          0             r_{0;21}*e             ]
[                                               ]
[          I*omega_21*t                         ]
[r_{0;21}*e                         0           ]



>>> fprint(r_I[2],print_ascii=print_ascii)
[                                     -I*omega_21*t]
[           0             -r_{-1;21}*e             ]
[                                                  ]
[           I*omega_21*t                           ]
[r_{+1;21}*e                          0            ]



Which can be decomposed as

>>> r_I_p=[ Matrix([[ delta_greater(j,i)*exp(-I*omega[j,i]*t)*r[p][i,j] for j in range(Ne)]for i in range(Ne)]) for p in range(3)]
>>> fprint(r_I_p[0],print_ascii=print_ascii)
[               -I*omega_21*t]
[0  -r_{+1;21}*e             ]
[                            ]
[0              0            ]



>>> fprint(r_I_p[1],print_ascii=print_ascii)
[             -I*omega_21*t]
[0  r_{0;21}*e             ]
[                          ]
[0             0           ]



>>> fprint(r_I_p[2],print_ascii=print_ascii)
[               -I*omega_21*t]
[0  -r_{-1;21}*e             ]
[                            ]
[0              0            ]



>>> r_I_m=[ Matrix([[ delta_lesser( j,i)*exp( I*omega[i,j]*t)*r[p][i,j] for j in range(Ne)]for i in range(Ne)]) for p in range(3)]
>>> fprint(r_I_m[0],print_ascii=print_ascii)
[           0             0]
[                          ]
[           I*omega_21*t   ]
[r_{-1;21}*e              0]



>>> fprint(r_I_m[1],print_ascii=print_ascii)
[          0             0]
[                         ]
[          I*omega_21*t   ]
[r_{0;21}*e              0]



>>> fprint(r_I_m[2],print_ascii=print_ascii)
[           0             0]
[                          ]
[           I*omega_21*t   ]
[r_{+1;21}*e              0]



that summed equal $\\vec{\\hat{r}}_I$

>>> fprint( [r_I[p]-(r_I_p[p]+r_I_m[p]) for p in range(3)] ,print_ascii=print_ascii)
[[0  0], [0  0], [0  0]]
 [    ]  [    ]  [    ] 
 [0  0]  [0  0]  [0  0] 



Thus the interaction hamiltonian in the interaciton picture is
\\begin{equation}
    \\hat{H}_{1I}=e\\vec{E}\\cdot \\vec{\\hat{r}}_I= e(\\vec{E}^{(+)}\\cdot \\vec{\\hat{r}}^{(+)}_I + \\vec{E}^{(+)}\\cdot \\vec{\\hat{r}}^{(-)}_I + \\vec{E}^{(-)}\\cdot \\vec{\\hat{r}}^{(+)}_I + \\vec{E}^{(-)}\\cdot \\vec{\\hat{r}}^{(-)}_I)
\\end{equation}

>>> H1I=e*helicity_dot_product(E,r_I)
>>> fprint(H1I,print_ascii=print_ascii)
[                                                                                /       -I*t*varpi_1    I*t*varpi_1 _
[                                                                                |E_0^1*e               e           *E
[                                 0                                   e*r_{0;21}*|------------------- + --------------
[                                                                                \         2                    2     
[                                                                                                                     
[           /       -I*t*varpi_1    I*t*varpi_1 _____\                                                                
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
>>> fprint(H1IRWA,print_ascii=print_ascii)
[                                                          -I*omega_21*t  I*t*varpi_1 _____]
[                                              e*r_{0;21}*e             *e           *E_0^1]
[                     0                        --------------------------------------------]
[                                                                   2                      ]
[                                                                                          ]
[                  I*omega_21*t  -I*t*varpi_1                                              ]
[E_0^1*e*r_{0;21}*e            *e                                                          ]
[--------------------------------------------                       0                      ]
[                     2                                                                    ]



 Returning to the Schrödinger picture we have.

>>> r_p=[ Matrix([[ delta_greater(j,i)*r[p][i,j] for j in range(Ne)]for i in range(Ne)]) for p in range(3)]
>>> fprint(r_p,print_ascii=print_ascii)
[[0  -r_{+1;21}], [0  r_{0;21}], [0  -r_{-1;21}]]
 [             ]  [           ]  [             ] 
 [0      0     ]  [0     0    ]  [0      0     ] 



>>> r_m=[ Matrix([[ delta_lesser( j,i)*r[p][i,j] for j in range(Ne)]for i in range(Ne)]) for p in range(3)]
>>> fprint(r_m,print_ascii=print_ascii)
[[    0      0], [   0      0], [    0      0]]
 [            ]  [           ]  [            ] 
 [r_{-1;21}  0]  [r_{0;21}  0]  [r_{+1;21}  0] 



>>> fprint( [r[p]-(r_p[p]+r_m[p]) for p in range(3)] ,print_ascii=print_ascii)
[[0  0], [0  0], [0  0]]
 [    ]  [    ]  [    ] 
 [0  0]  [0  0]  [0  0] 



Thus the interaction hamiltonian in the Schrödinger picture in the rotating wave approximation is

>>> H1RWA=e*(helicity_dot_product(E_p,r_m)+helicity_dot_product(E_m,r_p))
>>> fprint(H1RWA,print_ascii=print_ascii)
[                                            I*t*varpi_1 _____]
[                                e*r_{0;21}*e           *E_0^1]
[              0                 -----------------------------]
[                                              2              ]
[                                                             ]
[                  -I*t*varpi_1                               ]
[E_0^1*e*r_{0;21}*e                                           ]
[------------------------------                0              ]
[              2                                              ]



And the complete hamiltonian in the Schrödinger picture in the rotating wave approximation is

>>> HRWA=H0+H1RWA
>>> fprint(HRWA,print_ascii=print_ascii)
[                                            I*t*varpi_1 _____]
[                                e*r_{0;21}*e           *E_0^1]
[         hbar*omega_1           -----------------------------]
[                                              2              ]
[                                                             ]
[                  -I*t*varpi_1                               ]
[E_0^1*e*r_{0;21}*e                                           ]
[------------------------------          hbar*omega_2         ]
[              2                                              ]



# Rotating Frame
Next we will make a phase transformation in order to eliminate the explicit time dependance of the equations.

>>> c,ctilde,phase=define_psi_coefficients(Ne)
>>> fprint([c,ctilde,phase], print_ascii=print_ascii)
[[c1(t)], [\tilde{c}_{1}(t)], [theta1]]
 [     ]  [                ]  [      ] 
 [c2(t)]  [\tilde{c}_{2}(t)]  [theta2] 



>>> psi=Matrix([ exp(I*phase[i]*t)*ctilde[i] for i in range(Ne)])
>>> fprint(psi, print_ascii=print_ascii)
[                  I*t*theta1]
[\tilde{c}_{1}(t)*e          ]
[                            ]
[                  I*t*theta2]
[\tilde{c}_{2}(t)*e          ]



The Schrödinger equation $i\\hbar \\partial_t |\\psi\\rangle=\\hat{H}_{RWA}$ is

>>> lhs=Matrix([(I*hbar*Derivative(psi[i],t).doit()).expand() for i in range(Ne)])
>>> fprint(lhs, print_ascii=print_ascii)
[                               I*t*theta1]
[-hbar*theta1*\tilde{c}_{1}(t)*e          ]
[                                         ]
[                               I*t*theta2]
[-hbar*theta2*\tilde{c}_{2}(t)*e          ]



>>> rhs=HRWA*psi
>>> fprint(rhs, print_ascii=print_ascii)
[                             I*t*theta2  I*t*varpi_1 _____                                             ]
[e*r_{0;21}*\tilde{c}_{2}(t)*e          *e           *E_0^1                                  I*t*theta1 ]
[---------------------------------------------------------- + hbar*omega_1*\tilde{c}_{1}(t)*e           ]
[                            2                                                                          ]
[                                                                                                       ]
[                                   I*t*theta1  -I*t*varpi_1                                            ]
[E_0^1*e*r_{0;21}*\tilde{c}_{1}(t)*e          *e                                              I*t*theta2]
[----------------------------------------------------------- + hbar*omega_2*\tilde{c}_{2}(t)*e          ]
[                             2                                                                         ]



We multiply each of these equations by $e^{-i \\theta_i t}$ and substracting $i \\theta_i \\tilde{c}_i$

>>> lhs_new=Matrix([simplify(  lhs[i]*exp(-I*phase[i]*t) +hbar*phase[i]*ctilde[i] ) for i in range(Ne)])
>>> fprint(lhs_new, print_ascii=print_ascii)
[0]
[ ]
[0]



>>> rhs_new=Matrix([simplify(  rhs[i]*exp(-I*phase[i]*t) +hbar*phase[i]*ctilde[i] ) for i in range(Ne)])
>>> fprint(rhs_new, print_ascii=print_ascii)
[                             -I*t*theta1  I*t*theta2  I*t*varpi_1 _____                                              
[e*r_{0;21}*\tilde{c}_{2}(t)*e           *e          *e           *E_0^1                                              
[----------------------------------------------------------------------- + hbar*omega_1*\tilde{c}_{1}(t) + hbar*theta1
[                                   2                                                                                 
[                                                                                                                     
[                                   I*t*theta1  -I*t*theta2  -I*t*varpi_1                                             
[E_0^1*e*r_{0;21}*\tilde{c}_{1}(t)*e          *e           *e                                                         
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
>>> fprint(phase_transformation,print_ascii=print_ascii)
{theta2: theta1 - varpi_1}



There is a free parameter $\\theta_1$, which is to be expected, since state vetors $|\\psi\\rangle$ always have a global phase invariance

>>> fprint(psi.subs(phase_transformation),print_ascii=print_ascii)
[                        I*t*theta1      ]
[      \tilde{c}_{1}(t)*e                ]
[                                        ]
[                  I*t*(theta1 - varpi_1)]
[\tilde{c}_{2}(t)*e                      ]



Thus the equations become

>>> fprint(lhs_new, print_ascii=print_ascii)
[0]
[ ]
[0]



>>> rhs_new=simplify(rhs_new.subs(phase_transformation))
>>> fprint(rhs_new, print_ascii=print_ascii)
[                                 _____                                                ]
[     e*r_{0;21}*\tilde{c}_{2}(t)*E_0^1                                                ]
[     --------------------------------- + hbar*(omega_1 + theta1)*\tilde{c}_{1}(t)     ]
[                     2                                                                ]
[                                                                                      ]
[E_0^1*e*r_{0;21}*\tilde{c}_{1}(t)                                                     ]
[--------------------------------- + hbar*(omega_2 + theta1 - varpi_1)*\tilde{c}_{2}(t)]
[                2                                                                     ]



It can be seen that this is the Schrödinger equation derived from an effective hamiltonian $\\tilde{H}$

>>> Htilde=Matrix([ [Derivative(rhs_new[i],ctilde[j]).doit() for j in range(Ne)] for i in range(Ne)])
>>> fprint(Htilde, print_ascii=print_ascii)
[                                            _____         ]
[                                 e*r_{0;21}*E_0^1         ]
[hbar*(omega_1 + theta1)          ----------------         ]
[                                        2                 ]
[                                                          ]
[   E_0^1*e*r_{0;21}                                       ]
[   ----------------      hbar*(omega_2 + theta1 - varpi_1)]
[          2                                               ]



We can see that it is convenient to choose $\\theta_1=-\\omega_1$ to simplify the hamiltonian. Also, we can recognize $\\omega^1-\\omega_2+\\omega_1=\\delta$ as the detuning of the laser field relative to the atomic transition $\\omega_{21}=\\omega_2-\\omega_1$.

>>> delta=Symbol("delta",real=True)
>>> Htilde=Htilde.subs({phase[0]:-omega_level[0]}).subs({omega_laser[0]:delta+omega_level[1]-omega_level[0]})
>>> fprint(Htilde, print_ascii=print_ascii)
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
>>> fprint(Htilde,print_ascii=print_ascii)
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
>>> fprint( rho , print_ascii=print_ascii)
[rho11  rho12]
[            ]
[rho21  rho22]



The hamiltonian part of the equations is
\\begin{equation}
    \\dot{\\hat{\\rho}}=\\frac{i}{\\hbar}[\\hat{\\rho}, \\hat{\\tilde{H}}]
\\end{equation}

>>> hamiltonian_terms=(I/hbar*(rho*Htilde-Htilde*rho)).expand()
>>> fprint(hamiltonian_terms,print_ascii=print_ascii)
[                                 _____                                    _____           _____]
[         I*Omega*rho12   I*rho21*Omega                            I*rho11*Omega   I*rho22*Omega]
[         ------------- - -------------           -I*delta*rho12 + ------------- - -------------]
[               2               2                                        2               2      ]
[                                                                                               ]
[                                                                                  _____        ]
[  I*Omega*rho11   I*Omega*rho22                           I*Omega*rho12   I*rho21*Omega        ]
[- ------------- + ------------- + I*delta*rho21         - ------------- + -------------        ]
[        2               2                                       2               2              ]



There is only one Lindblad operator, since there is only one spontaneous decay channel.

>>> lindblad_terms=gamma[1,0]*lindblad_operator(ket(1,Ne)*bra(2,Ne),rho)
>>> fprint(lindblad_terms, print_ascii=print_ascii)
[                  -gamma_21*rho12 ]
[ gamma_21*rho22   ----------------]
[                         2        ]
[                                  ]
[-gamma_21*rho21                   ]
[----------------  -gamma_21*rho22 ]
[       2                          ]



# Optical Bloch Equations
$\\textit{The}$ Optical Bloch equations are thus.

>>> eqs=hamiltonian_terms + lindblad_terms
>>> fprint(eqs,print_ascii=print_ascii)
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
>>> fprint( rho ,print_ascii=print_ascii)
[            _____]
[-rho22 + 1  rho21]
[                 ]
[  rho21     rho22]



>>> hamiltonian_terms = (I/hbar*(rho*Htilde-Htilde*rho)).expand()
>>> lindblad_terms    =gamma[1,0]*lindblad_operator(ket(1,Ne)*bra(2,Ne),rho)
>>> eqs=hamiltonian_terms + lindblad_terms
>>> fprint(eqs,print_ascii=print_ascii)
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
>>> fprint( re(eqs[1,1].subs(ss_comp)) ,print_ascii=print_ascii)
-gamma_21*rho22 - re(Omega)*im(rho21) + re(rho21)*im(Omega)



>>> fprint( re(eqs[1,0].subs(ss_comp)) ,print_ascii=print_ascii)
                   gamma_21*re(rho21)                     im(Omega)
-delta*im(rho21) - ------------------ - rho22*im(Omega) + ---------
                           2                                  2    



>>> fprint( im(eqs[1,0].subs(ss_comp)) ,print_ascii=print_ascii)
                  gamma_21*im(rho21)                     re(Omega)
delta*re(rho21) - ------------------ + rho22*re(Omega) - ---------
                          2                                  2    



If the density matrix is represented as a vector whose components are the these independent components of the density matrix

>>> rho_vect=define_rho_vector(rho,Ne)
>>> fprint(rho_vect,print_ascii=print_ascii)
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
>>> fprint([A,b],print_ascii=print_ascii)
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
>>> fprint(eqs_new,print_ascii=print_ascii)
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

>>> fprint( eqs_new - Matrix([re(eqs[1,1]),re(eqs[1,0]),im(eqs[1,0])]).subs(ss_comp) ,print_ascii=print_ascii)
[0]
[ ]
[0]
[ ]
[0]



The steady state solution of this equations is

>>> sol=solve(list(eqs_new),list(rho_vect))
>>> fprint( {rho_vect[0]:sol[rho_vect[0]]} ,print_ascii=print_ascii)
                       2            2                      
                     re (Omega) + im (Omega)               
{rho22: --------------------------------------------------}
               2           2       2              2        
        4*delta  + gamma_21  + 2*re (Omega) + 2*im (Omega) 



>>> fprint( {rho_vect[1]:sol[rho_vect[1]]} ,print_ascii=print_ascii)
                  2*delta*re(Omega) + gamma_21*im(Omega)       
{re(rho21): --------------------------------------------------}
                   2           2       2              2        
            4*delta  + gamma_21  + 2*re (Omega) + 2*im (Omega) 



>>> fprint( {rho_vect[2]:sol[rho_vect[2]]} ,print_ascii=print_ascii)
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
    
>>> fprint( [testi.subs({Omega:re(Omega)+I*im(Omega)}).factor() for testi in test] ,print_ascii=print_ascii)
[0, 0, 0]



So our derivation produces the same results as the literature. The saturation intensity is defined as the intensity needed to accumulate $\\frac{1}{4}$ of the population in the excited state when the field is in resonance ($\\delta=0$).

>>> saturation_eq=sol[rho[1,1]].subs({delta:0})-1/Integer(4)
>>> fprint( saturation_eq ,print_ascii=print_ascii)
          2            2                   
        re (Omega) + im (Omega)           1
--------------------------------------- - -
        2       2              2          4
gamma_21  + 2*re (Omega) + 2*im (Omega)    



>>> Omega_amp,alpha=symbols("\Omega_a alpha",real=True)
>>> ss={Omega:Omega_amp*cos(alpha)+I*Omega_amp*sin(alpha)}
>>> saturation_eq= saturation_eq.subs(ss).factor().simplify()
>>> fprint(saturation_eq,print_ascii=print_ascii)
            2           2  
  2*\Omega_a  - gamma_21   
---------------------------
  /          2           2\
4*\2*\Omega_a  + gamma_21 /



>>> Omega_sat=solve( saturation_eq ,Omega_amp)[1]
>>> fprint(Omega_sat, print_ascii=print_ascii)
  ___         
\/ 2 *gamma_21
--------------
      2       



Since $\\Omega =e E_0^1 r_{0;21}/\\hbar$ it follows that

>>> E0_sat=Omega_sat*hbar/e/r[1][1,0]
>>> fprint(E0_sat, print_ascii=print_ascii)
  ___              
\/ 2 *gamma_21*hbar
-------------------
    2*e*r_{0;21}   



The full width at half maximum of $\\rho_{22}$ is

>>> hm1,hm2=solve(sol[rho[1,1]]-sol[rho[1,1]].subs({delta:0})/2,delta)
>>> FWHM=hm2-hm1
>>> FWHM=FWHM.subs(ss).simplify()
>>> fprint(FWHM, print_ascii=print_ascii)
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
