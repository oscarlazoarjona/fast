# -*- coding: utf-8 -*-
# Copyright (C) 2017 Oscar Gerardo Lazo Arjona
# mailto: oscar.lazo@correo.nucleares.unam.mx

__doc__ = r"""

>>> from fast import *
>>> from matplotlib import pyplot
>>> from sympy import sin,cos,exp,sqrt,pi,zeros,I
>>> from numpy import array

>>> init_printing()
>>> print_ascii=True
>>> # print_ascii = False

>>> path="folder_07___Three_level_atom_Lambda_symbolic/" 
>>> name='suite'

We will be deriving the optical Bloch equations for a three level system in a Lambda configuration as that in the figure.

>>> fig=pyplot.figure(); ax=fig.add_subplot(111,aspect="equal")
    
>>> p1=[0.5,1]; p2=[2.5,1.2]; p3=[1.5,3.5]
>>> draw_state(ax,p1,text=r"$|1\rangle$",l=1.0,alignment='right',label_displacement=0.05,fontsize=25,linewidth=4.0)
>>> draw_state(ax,p2,text=r"$|2\rangle$",l=1.0,alignment='right',label_displacement=0.05,fontsize=25,linewidth=4.0)
>>> draw_state(ax,p3,text=r"$|3\rangle$",l=1.0,alignment='right',label_displacement=0.05,fontsize=25,linewidth=4.0)
    
>>> excitation(ax,[p1[0]+0.25,p1[1]],[p3[0]+0.10,p3[1]], fc="r", ec="r",width=0.01, head_width=0.2, head_length=0.2) # doctest: +IGNORE_PLOT_STEP4
>>> excitation(ax,[p2[0]+0.25,p2[1]],[p3[0]+0.25,p3[1]], fc="b", ec="b",width=0.01, head_width=0.2, head_length=0.2) # doctest: +IGNORE_PLOT_STEP4
    
>>> decay(     ax,[p1[0]-0.25,p1[1]],[p3[0]-0.35,p3[1]], 0.05,10.0,color="r",linewidth=1.0) # doctest: +IGNORE_PLOT_STEP4
>>> decay(     ax,[p2[0]-0.25,p2[1]],[p3[0]-0.25,p3[1]], 0.05,10.0,color="b",linewidth=1.0) # doctest: +IGNORE_PLOT_STEP4
    
>>> pyplot.axis('off') # doctest: +IGNORE_PLOT_STEP3
>>> pyplot.savefig(path+name+'_diagram.png',bbox_inches="tight") # doctest: +IGNORE_PLOT_STEP4
<Figure size 432x288 with 1 Axes>



We define the number of states and of radiation fields.

>>> Ne=3
>>> Nl=2

We define the variables related to the laser field.

>>> E0,omega_laser=define_laser_variables(Nl)
>>> fprint(E0,print_ascii=print_ascii)
[E_{01}, E_{02}]



>>> fprint(omega_laser, print_ascii=print_ascii)
[varpi_1, varpi_2]



We define a few important symbols.

>>> t,hbar,e=symbols("t hbar e",positive=True)
>>> fprint([t,hbar,e], print_ascii=print_ascii)
[t, hbar, e]



We write two electric fields propagating trough the $\\hat{x}$ direction polarized in the $\\hat{z}$ direction. First the wave vectors:

>>> phi1=0 ; theta1=pi/2; alpha1=pi/2; beta1=0
>>> phi2=pi; theta2=pi/2; alpha2=pi/2; beta2=0
    
>>> k1=Matrix([cos(phi1)*sin(theta1),sin(phi1)*sin(theta1),cos(theta1)])
>>> k2=Matrix([cos(phi2)*sin(theta2),sin(phi2)*sin(theta2),cos(theta2)])
    
>>> fprint([k1,k2],print_ascii=print_ascii)
 [1]  [-1] 
 [ ]  [  ] 
[[0], [0 ]]
 [ ]  [  ] 
 [0]  [0 ] 



The polarization vectors.

>>> ep1=polarization_vector(phi1,theta1,alpha1,beta1, 1)
>>> ep2=polarization_vector(phi2,theta2,alpha2,beta2, 1)
    
>>> em1=ep1.conjugate()
>>> em2=ep2.conjugate()
    
>>> ep=[ep1,ep2]
>>> em=[em1,em2]
    
>>> fprint([ep,em], print_ascii=print_ascii)
  [0]  [0]    [0]  [0]  
  [ ]  [ ]    [ ]  [ ]  
[[[0], [0]], [[0], [0]]]
  [ ]  [ ]    [ ]  [ ]  
  [1]  [1]    [1]  [1]  



The electric field (evaluated in $\\vec{R}=0$).

>>> zero_vect=Matrix([0,0,0])
>>> E_cartesian = [(E0[l]*ep[l]*exp(-I*omega_laser[l]*t) + E0[l].conjugate()*em[l]*exp( I*omega_laser[l]*t))/2 
...                     for l in range(Nl)]
    
>>> fprint(E_cartesian, print_ascii=print_ascii)
 [                    0                     ]  [                    0                     ] 
 [                                          ]  [                                          ] 
 [                    0                     ]  [                    0                     ] 
 [                                          ]  [                                          ] 
[[        -I*t*varpi_1    I*t*varpi_1 ______], [        -I*t*varpi_2    I*t*varpi_2 ______]]
 [E_{01}*e               e           *E_{01}]  [E_{02}*e               e           *E_{02}] 
 [-------------------- + -------------------]  [-------------------- + -------------------] 
 [         2                      2         ]  [         2                      2         ] 



>>> l1=PlaneWave(phi1,theta1,alpha1,beta1,color="red")
>>> l2=PlaneWave(phi2,theta2,alpha2,beta2,color="blue")
    
>>> laseres=[l1,l2]
>>> Nl=len(laseres)
    
>>> fig = pyplot.figure(); ax = fig.gca(projection='3d')
>>> draw_lasers_3d(ax,laseres,path+'lasers.png') # doctest: +IGNORE_PLOT_STEP4
<Figure size 432x288 with 1 Axes>



We write the electric field in the helicity basis (see notebook "Vectors in the helicity basis and the electric field").

>>> E=[cartesian_to_helicity(E_cartesian[l]) for l in range(Nl)]
>>> fprint(E, print_ascii=print_ascii)
 [                    0                     ]  [                    0                     ] 
 [                                          ]  [                                          ] 
 [        -I*t*varpi_1    I*t*varpi_1 ______]  [        -I*t*varpi_2    I*t*varpi_2 ______] 
 [E_{01}*e               e           *E_{01}]  [E_{02}*e               e           *E_{02}] 
[[-------------------- + -------------------], [-------------------- + -------------------]]
 [         2                      2         ]  [         2                      2         ] 
 [                                          ]  [                                          ] 
 [                    0                     ]  [                    0                     ] 



We define the position operator.

>>> r=define_r_components(Ne,helicity=True,explicitly_hermitian=True)
>>> #Ladder means that r_{p;31}=0
>>> r=[ri.subs({r[0][1,0]:0,r[1][1,0]:0,r[2][1,0]:0}) for ri in r]
>>> fprint(r, print_ascii=print_ascii)
 [    0          0      -r_{+1;31}]  [   0         0      r_{0;31}]  [    0          0      -r_{-1;31}] 
 [                                ]  [                            ]  [                                ] 
[[    0          0      -r_{+1;32}], [   0         0      r_{0;32}], [    0          0      -r_{-1;32}]]
 [                                ]  [                            ]  [                                ] 
 [r_{-1;31}  r_{-1;32}      0     ]  [r_{0;31}  r_{0;32}     0    ]  [r_{+1;31}  r_{+1;32}      0     ] 



The frequencies of the energy levels, the resonant frequencies, and the decay frequencies.

>>> omega_level,omega,gamma=define_frequencies(Ne,explicitly_antisymmetric=True)
>>> #Ladder means gamma31=0
>>> gamma=gamma.subs({gamma[1,0]:0})
    
>>> fprint(omega_level, print_ascii=print_ascii)
[omega_1, omega_2, omega_3]



>>> fprint(omega, print_ascii=print_ascii)
[   0      -omega_21  -omega_31]
[                              ]
[omega_21      0      -omega_32]
[                              ]
[omega_31  omega_32       0    ]



>>> fprint(gamma, print_ascii=print_ascii)
[   0         0      -gamma_31]
[                             ]
[   0         0      -gamma_32]
[                             ]
[gamma_31  gamma_32      0    ]



The atomic hamiltonian is

>>> H0=Matrix([[hbar*omega_level[i]*KroneckerDelta(i,j) for j in range(Ne)] for i in range(Ne)])
>>> fprint(H0, print_ascii=print_ascii)
[hbar*omega_1       0             0      ]
[                                        ]
[     0        hbar*omega_2       0      ]
[                                        ]
[     0             0        hbar*omega_3]



The interaction hamiltonian is

>>> zero_matrix=zeros(Ne,Ne)
>>> H1=sum([ e*helicity_dot_product(E[l],r) for l in range(Nl)],zero_matrix)
>>> fprint(H1,print_ascii=print_ascii)
[                                                                                                                     
[                                                                                                                     
[                                                        0                                                            
[                                                                                                                     
[                                                                                                                     
[                                                                                                                     
[                                                                                                                     
[                                                        0                                                            
[                                                                                                                     
[                                                                                                                     
[           /        -I*t*varpi_1    I*t*varpi_1 ______\              /        -I*t*varpi_2    I*t*varpi_2 ______\    
[           |E_{01}*e               e           *E_{01}|              |E_{02}*e               e           *E_{02}|    
[e*r_{0;31}*|-------------------- + -------------------| + e*r_{0;31}*|-------------------- + -------------------|  e*
[           \         2                      2         /              \         2                      2         /    
<BLANKLINE>
<BLANKLINE>
<BLANKLINE>
                                                      0                                                          e*r_{
<BLANKLINE>
<BLANKLINE>
<BLANKLINE>
<BLANKLINE>
                                                      0                                                          e*r_{
<BLANKLINE>
<BLANKLINE>
         /        -I*t*varpi_1    I*t*varpi_1 ______\              /        -I*t*varpi_2    I*t*varpi_2 ______\       
         |E_{01}*e               e           *E_{01}|              |E_{02}*e               e           *E_{02}|       
r_{0;32}*|-------------------- + -------------------| + e*r_{0;32}*|-------------------- + -------------------|       
         \         2                      2         /              \         2                      2         /       
<BLANKLINE>
      /        -I*t*varpi_1    I*t*varpi_1 ______\              /        -I*t*varpi_2    I*t*varpi_2 ______\]
      |E_{01}*e               e           *E_{01}|              |E_{02}*e               e           *E_{02}|]
0;31}*|-------------------- + -------------------| + e*r_{0;31}*|-------------------- + -------------------|]
      \         2                      2         /              \         2                      2         /]
                                                                                                            ]
      /        -I*t*varpi_1    I*t*varpi_1 ______\              /        -I*t*varpi_2    I*t*varpi_2 ______\]
      |E_{01}*e               e           *E_{01}|              |E_{02}*e               e           *E_{02}|]
0;32}*|-------------------- + -------------------| + e*r_{0;32}*|-------------------- + -------------------|]
      \         2                      2         /              \         2                      2         /]
                                                                                                            ]
                                                                                                            ]
                                                                                                            ]
                                                   0                                                        ]
                                                                                                            ]



and the complete hamiltonian is

>>> H=H0+H1
>>> fprint(H,print_ascii=print_ascii)
[                                                                                                                     
[                                                                                                                     
[                                                  hbar*omega_1                                                       
[                                                                                                                     
[                                                                                                                     
[                                                                                                                     
[                                                                                                                     
[                                                        0                                                            
[                                                                                                                     
[                                                                                                                     
[           /        -I*t*varpi_1    I*t*varpi_1 ______\              /        -I*t*varpi_2    I*t*varpi_2 ______\    
[           |E_{01}*e               e           *E_{01}|              |E_{02}*e               e           *E_{02}|    
[e*r_{0;31}*|-------------------- + -------------------| + e*r_{0;31}*|-------------------- + -------------------|  e*
[           \         2                      2         /              \         2                      2         /    
<BLANKLINE>
<BLANKLINE>
<BLANKLINE>
                                                      0                                                          e*r_{
<BLANKLINE>
<BLANKLINE>
<BLANKLINE>
<BLANKLINE>
                                                hbar*omega_2                                                     e*r_{
<BLANKLINE>
<BLANKLINE>
         /        -I*t*varpi_1    I*t*varpi_1 ______\              /        -I*t*varpi_2    I*t*varpi_2 ______\       
         |E_{01}*e               e           *E_{01}|              |E_{02}*e               e           *E_{02}|       
r_{0;32}*|-------------------- + -------------------| + e*r_{0;32}*|-------------------- + -------------------|       
         \         2                      2         /              \         2                      2         /       
<BLANKLINE>
      /        -I*t*varpi_1    I*t*varpi_1 ______\              /        -I*t*varpi_2    I*t*varpi_2 ______\]
      |E_{01}*e               e           *E_{01}|              |E_{02}*e               e           *E_{02}|]
0;31}*|-------------------- + -------------------| + e*r_{0;31}*|-------------------- + -------------------|]
      \         2                      2         /              \         2                      2         /]
                                                                                                            ]
      /        -I*t*varpi_1    I*t*varpi_1 ______\              /        -I*t*varpi_2    I*t*varpi_2 ______\]
      |E_{01}*e               e           *E_{01}|              |E_{02}*e               e           *E_{02}|]
0;32}*|-------------------- + -------------------| + e*r_{0;32}*|-------------------- + -------------------|]
      \         2                      2         /              \         2                      2         /]
                                                                                                            ]
                                                                                                            ]
                                                                                                            ]
                                             hbar*omega_3                                                   ]
                                                                                                            ]



# Rotating wave approximation
Notice that the electric field can be separated by terms with positive and negative frequency:

>>> E_cartesian_p=[E0[l]/2*ep[l]*exp(-I*omega_laser[l]*t) for l in range(Nl)]
>>> E_cartesian_m=[E0[l].conjugate()/2*em[l]*exp( I*omega_laser[l]*t) for l in range(Nl)]
    
>>> E_p=[cartesian_to_helicity(E_cartesian_p[l]) for l in range(Nl)]
>>> E_m=[cartesian_to_helicity(E_cartesian_m[l]) for l in range(Nl)]
    
>>> fprint([E_p,E_m], print_ascii=print_ascii)
  [         0          ]  [         0          ]    [         0         ]  [         0         ]  
  [                    ]  [                    ]    [                   ]  [                   ]  
  [        -I*t*varpi_1]  [        -I*t*varpi_2]    [ I*t*varpi_1 ______]  [ I*t*varpi_2 ______]  
  [E_{01}*e            ]  [E_{02}*e            ]    [e           *E_{01}]  [e           *E_{02}]  
[[[--------------------], [--------------------]], [[-------------------], [-------------------]]]
  [         2          ]  [         2          ]    [         2         ]  [         2         ]  
  [                    ]  [                    ]    [                   ]  [                   ]  
  [         0          ]  [         0          ]    [         0         ]  [         0         ]  



>>> fprint( simplify(sum([E[l] for l in range(Nl)],zero_vect)-(sum([E_p[l]+E_m[l] for l in range(Nl)],zero_vect) )), print_ascii=print_ascii)
[0]
[ ]
[0]
[ ]
[0]



The position operator can also be separated in this way. We go to the interaction picture (with $\\hat{H}_0$ as the undisturbed hamiltonian)

>>> r_I=[ Matrix([[exp(I*omega[i,j]*t)*r[p][i,j] for j in range(Ne)] for i in range(Ne)]) for p in range(3)]
>>> fprint(r_I[0], print_ascii=print_ascii)
[                                                              -I*omega_31*t]
[           0                        0             -r_{+1;31}*e             ]
[                                                                           ]
[                                                              -I*omega_32*t]
[           0                        0             -r_{+1;32}*e             ]
[                                                                           ]
[           I*omega_31*t             I*omega_32*t                           ]
[r_{-1;31}*e              r_{-1;32}*e                          0            ]



>>> fprint(r_I[1], print_ascii=print_ascii)
[                                                          -I*omega_31*t]
[          0                       0             r_{0;31}*e             ]
[                                                                       ]
[                                                          -I*omega_32*t]
[          0                       0             r_{0;32}*e             ]
[                                                                       ]
[          I*omega_31*t            I*omega_32*t                         ]
[r_{0;31}*e              r_{0;32}*e                         0           ]



>>> fprint(r_I[2], print_ascii=print_ascii)
[                                                              -I*omega_31*t]
[           0                        0             -r_{-1;31}*e             ]
[                                                                           ]
[                                                              -I*omega_32*t]
[           0                        0             -r_{-1;32}*e             ]
[                                                                           ]
[           I*omega_31*t             I*omega_32*t                           ]
[r_{+1;31}*e              r_{+1;32}*e                          0            ]



Which can be decomposed in positive and negative frequencies as

>>> r_I_p=[ Matrix([[ delta_greater(j,i)*exp(-I*omega[j,i]*t)*r[p][i,j] for j in range(Ne)]for i in range(Ne)]) for p in range(3)]
>>> fprint(r_I_p[0], print_ascii=print_ascii)
[                  -I*omega_31*t]
[0  0  -r_{+1;31}*e             ]
[                               ]
[                  -I*omega_32*t]
[0  0  -r_{+1;32}*e             ]
[                               ]
[0  0              0            ]



>>> fprint(r_I_p[1], print_ascii=print_ascii)
[                -I*omega_31*t]
[0  0  r_{0;31}*e             ]
[                             ]
[                -I*omega_32*t]
[0  0  r_{0;32}*e             ]
[                             ]
[0  0             0           ]



>>> fprint(r_I_p[2], print_ascii=print_ascii)
[                  -I*omega_31*t]
[0  0  -r_{-1;31}*e             ]
[                               ]
[                  -I*omega_32*t]
[0  0  -r_{-1;32}*e             ]
[                               ]
[0  0              0            ]



>>> r_I_m=[ Matrix([[ delta_lesser( j,i)*exp( I*omega[i,j]*t)*r[p][i,j] for j in range(Ne)]for i in range(Ne)]) for p in range(3)]
>>> fprint(r_I_m[0], print_ascii=print_ascii)
[           0                        0             0]
[                                                   ]
[           0                        0             0]
[                                                   ]
[           I*omega_31*t             I*omega_32*t   ]
[r_{-1;31}*e              r_{-1;32}*e              0]



>>> fprint(r_I_m[1], print_ascii=print_ascii)
[          0                       0             0]
[                                                 ]
[          0                       0             0]
[                                                 ]
[          I*omega_31*t            I*omega_32*t   ]
[r_{0;31}*e              r_{0;32}*e              0]



>>> fprint(r_I_m[2], print_ascii=print_ascii)
[           0                        0             0]
[                                                   ]
[           0                        0             0]
[                                                   ]
[           I*omega_31*t             I*omega_32*t   ]
[r_{+1;31}*e              r_{+1;32}*e              0]



that summed equal $\\vec{\\hat{r}}_I$

>>> fprint( [r_I[p]-(r_I_p[p]+r_I_m[p]) for p in range(3)] , print_ascii=print_ascii)
 [0  0  0]  [0  0  0]  [0  0  0] 
 [       ]  [       ]  [       ] 
[[0  0  0], [0  0  0], [0  0  0]]
 [       ]  [       ]  [       ] 
 [0  0  0]  [0  0  0]  [0  0  0] 



Thus the interaction hamiltonian in the interaciton picture is
\\begin{equation}
    \\hat{H}_{1I}=e\\vec{E}\\cdot \\vec{\\hat{r}}_I= e(\\vec{E}^{(+)}\\cdot \\vec{\\hat{r}}^{(+)}_I + \\vec{E}^{(+)}\\cdot \\vec{\\hat{r}}^{(-)}_I + \\vec{E}^{(-)}\\cdot \\vec{\\hat{r}}^{(+)}_I + \\vec{E}^{(-)}\\cdot \\vec{\\hat{r}}^{(-)}_I)
\\end{equation}

>>> H1I=sum([ e*helicity_dot_product(E[l],r_I) for l in range(Nl)],zero_matrix)
>>> fprint(H1I, print_ascii=print_ascii)
[                                                                                                                     
[                                                                                                                     
[                                                                      0                                              
[                                                                                                                     
[                                                                                                                     
[                                                                                                                     
[                                                                                                                     
[                                                                      0                                              
[                                                                                                                     
[                                                                                                                     
[           /        -I*t*varpi_1    I*t*varpi_1 ______\                            /        -I*t*varpi_2    I*t*varpi
[           |E_{01}*e               e           *E_{01}|  I*omega_31*t              |E_{02}*e               e         
[e*r_{0;31}*|-------------------- + -------------------|*e             + e*r_{0;31}*|-------------------- + ----------
[           \         2                      2         /                            \         2                      2
<BLANKLINE>
<BLANKLINE>
<BLANKLINE>
                                                                                                0                     
<BLANKLINE>
<BLANKLINE>
<BLANKLINE>
<BLANKLINE>
                                                                                                0                     
<BLANKLINE>
<BLANKLINE>
_2 ______\                           /        -I*t*varpi_1    I*t*varpi_1 ______\                            /        
  *E_{02}|  I*omega_31*t             |E_{01}*e               e           *E_{01}|  I*omega_32*t              |E_{02}*e
---------|*e              e*r_{0;32}*|-------------------- + -------------------|*e             + e*r_{0;32}*|--------
         /                           \         2                      2         /                            \        
<BLANKLINE>
                                                              /        -I*t*varpi_1    I*t*varpi_1 ______\            
                                                              |E_{01}*e               e           *E_{01}|  -I*omega_3
                                                   e*r_{0;31}*|-------------------- + -------------------|*e          
                                                              \         2                      2         /            
<BLANKLINE>
                                                              /        -I*t*varpi_1    I*t*varpi_1 ______\            
                                                              |E_{01}*e               e           *E_{01}|  -I*omega_3
                                                   e*r_{0;32}*|-------------------- + -------------------|*e          
                                                              \         2                      2         /            
<BLANKLINE>
-I*t*varpi_2    I*t*varpi_2 ______\                                                                                   
               e           *E_{02}|  I*omega_32*t                                                                     
------------ + -------------------|*e                                                                                 
 2                      2         /                                                                                   
<BLANKLINE>
                 /        -I*t*varpi_2    I*t*varpi_2 ______\               ]
1*t              |E_{02}*e               e           *E_{02}|  -I*omega_31*t]
    + e*r_{0;31}*|-------------------- + -------------------|*e             ]
                 \         2                      2         /               ]
                                                                            ]
                 /        -I*t*varpi_2    I*t*varpi_2 ______\               ]
2*t              |E_{02}*e               e           *E_{02}|  -I*omega_32*t]
    + e*r_{0;32}*|-------------------- + -------------------|*e             ]
                 \         2                      2         /               ]
                                                                            ]
                                                                            ]
                                                                            ]
    0                                                                       ]
                                                                            ]



Since both $\\omega^l$ and $\\omega_{ij}$ are in the order of THz, the terms that have frequencies with the same sign are summed, and thus also of the order of THz. The frequencies in the terms with oposite signs however, are detunings of the order of MHz. Since we are only interested in the coarse-grained evolution of the density matrix, we may omit the fast terms and approximate

\\begin{equation}
    \\hat{H}_{1I} \\simeq \\hat{H}_{1I,RWA}= e( \\vec{E}^{(+)}\\cdot \\vec{\\hat{r}}^{(-)}_I + \\vec{E}^{(-)}\\cdot \\vec{\\hat{r}}^{(+)}_I )
\\end{equation}

That is known as the rotating wave approximation (RWA).

>>> H1IRWA=sum( [ (e*(helicity_dot_product(E_p[l],r_I_m)+helicity_dot_product(E_m[l],r_I_p))) for l in range(Nl)],zero_matrix)
>>> fprint(H1IRWA, print_ascii=print_ascii)
[                                                                                                                     
[                                                                                                                     
[                                              0                                                                      
[                                                                                                                     
[                                                                                                                     
[                                                                                                                     
[                                                                                                                     
[                                              0                                                                      
[                                                                                                                     
[                                                                                                                     
[                   I*omega_31*t  -I*t*varpi_1                      I*omega_31*t  -I*t*varpi_2                     I*o
[E_{01}*e*r_{0;31}*e            *e               E_{02}*e*r_{0;31}*e            *e              E_{01}*e*r_{0;32}*e   
[--------------------------------------------- + ---------------------------------------------  ----------------------
[                      2                                               2                                              
<BLANKLINE>
                                                                                     -I*omega_31*t  I*t*varpi_1 ______
                                                                         e*r_{0;31}*e             *e           *E_{01}
                        0                                                ---------------------------------------------
                                                                                               2                      
<BLANKLINE>
                                                                                     -I*omega_32*t  I*t*varpi_1 ______
                                                                         e*r_{0;32}*e             *e           *E_{01}
                        0                                                ---------------------------------------------
                                                                                               2                      
<BLANKLINE>
mega_32*t  -I*t*varpi_1                      I*omega_32*t  -I*t*varpi_2                                               
         *e               E_{02}*e*r_{0;32}*e            *e                                                           
----------------------- + ---------------------------------------------                                               
2                                               2                                                                     
<BLANKLINE>
               -I*omega_31*t  I*t*varpi_2 ______]
   e*r_{0;31}*e             *e           *E_{02}]
 + ---------------------------------------------]
                         2                      ]
                                                ]
               -I*omega_32*t  I*t*varpi_2 ______]
   e*r_{0;32}*e             *e           *E_{02}]
 + ---------------------------------------------]
                         2                      ]
                                                ]
                                                ]
                                                ]
 0                                              ]
                                                ]



The matrix element $(\\hat{H}_{1I,RWA})_{31}$ element is

>>> fprint(H1IRWA[2,0].expand(), print_ascii=print_ascii)
                   I*omega_31*t  -I*t*varpi_1                      I*omega_31*t  -I*t*varpi_2
E_{01}*e*r_{0;31}*e            *e               E_{02}*e*r_{0;31}*e            *e            
--------------------------------------------- + ---------------------------------------------
                      2                                               2                      



But if the detuning $\\omega_{31}-\\omega^1 \\ll \\omega_{31}-\\omega^2$ (the second field is far detuned from the $1 \\rightarrow 3$ transition), then $\\omega_{31}-\\omega^2$ may be also considered too high a frequency to be relevant to coarse-grained evolution. So we might neclect that term in $(\\hat{H}_{1I,RWA})_{31}$ and similarly neglect the $\\omega_{32}-\\omega^1$ for term in $(\\hat{H}_{1I,RWA})_{32}$:

>>> fprint(H1IRWA[2,1].expand(), print_ascii=print_ascii)
                   I*omega_32*t  -I*t*varpi_1                      I*omega_32*t  -I*t*varpi_2
E_{01}*e*r_{0;32}*e            *e               E_{02}*e*r_{0;32}*e            *e            
--------------------------------------------- + ---------------------------------------------
                      2                                               2                      



In other words, if the detunings in our experiments allow the approximmation, we might choose which frequency components $\\omega^l$ excite which transitions. Let us say that $L_{ij}$ is the set of $l$ such that $\\omega^l$ excites the transition $i\\rightarrow j$

>>> Lij=[[1,3,[1]],[2,3,[2]]]
>>> Lij=formatLij(Lij,Ne)
>>> print array(Lij)
[[list([]) list([]) list([1])]
 [list([]) list([]) list([2])]
 [list([1]) list([2]) list([])]]



Thus the interacion hamiltonian in the interaction picture can be approximated as

>>> H1IRWA =sum([ e*( helicity_dot_product( E_p[l],vector_element(r_I_m,i,j)) ) * ket(i+1,Ne)*bra(j+1,Ne) 
...             for l in range(Nl) for j in range(Ne) for i in range(Ne) if l+1 in Lij[i][j] ],zero_matrix)
>>> H1IRWA+=sum([ e*( helicity_dot_product( E_m[l],vector_element(r_I_p,i,j)) ) * ket(i+1,Ne)*bra(j+1,Ne) 
...             for l in range(Nl) for j in range(Ne) for i in range(Ne) if l+1 in Lij[i][j] ],zero_matrix)
    
>>> fprint(H1IRWA, print_ascii=print_ascii)
[                                                                                                          -I*omega_31
[                                                                                              e*r_{0;31}*e           
[                      0                                              0                        -----------------------
[                                                                                                                    2
[                                                                                                                     
[                                                                                                          -I*omega_32
[                                                                                              e*r_{0;32}*e           
[                      0                                              0                        -----------------------
[                                                                                                                    2
[                                                                                                                     
[                   I*omega_31*t  -I*t*varpi_1                     I*omega_32*t  -I*t*varpi_2                         
[E_{01}*e*r_{0;31}*e            *e              E_{02}*e*r_{0;32}*e            *e                                     
[---------------------------------------------  ---------------------------------------------                        0
[                      2                                              2                                               
<BLANKLINE>
*t  I*t*varpi_1 ______]
  *e           *E_{01}]
----------------------]
                      ]
                      ]
*t  I*t*varpi_2 ______]
  *e           *E_{02}]
----------------------]
                      ]
                      ]
                      ]
                      ]
                      ]
                      ]



Returning to the Schrödinger picture we have.

>>> r_p=[ Matrix([[ delta_greater(j,i)*r[p][i,j] for j in range(Ne)]for i in range(Ne)]) for p in range(3)]
>>> fprint(r_p, print_ascii=print_ascii)
 [0  0  -r_{+1;31}]  [0  0  r_{0;31}]  [0  0  -r_{-1;31}] 
 [                ]  [              ]  [                ] 
[[0  0  -r_{+1;32}], [0  0  r_{0;32}], [0  0  -r_{-1;32}]]
 [                ]  [              ]  [                ] 
 [0  0      0     ]  [0  0     0    ]  [0  0      0     ] 



>>> r_m=[ Matrix([[ delta_lesser( j,i)*r[p][i,j] for j in range(Ne)]for i in range(Ne)]) for p in range(3)]
>>> fprint(r_m, print_ascii=print_ascii)
 [    0          0      0]  [   0         0      0]  [    0          0      0] 
 [                       ]  [                     ]  [                       ] 
[[    0          0      0], [   0         0      0], [    0          0      0]]
 [                       ]  [                     ]  [                       ] 
 [r_{-1;31}  r_{-1;32}  0]  [r_{0;31}  r_{0;32}  0]  [r_{+1;31}  r_{+1;32}  0] 



>>> fprint( [r[p]-(r_p[p]+r_m[p]) for p in range(3)] , print_ascii=print_ascii)
 [0  0  0]  [0  0  0]  [0  0  0] 
 [       ]  [       ]  [       ] 
[[0  0  0], [0  0  0], [0  0  0]]
 [       ]  [       ]  [       ] 
 [0  0  0]  [0  0  0]  [0  0  0] 



Thus the interaction hamiltonian in the Schrödinger picture in the rotating wave approximation is

>>> H1RWA =sum([ e*( helicity_dot_product( E_p[l],vector_element(r_m,i,j)) ) * ket(i+1,Ne)*bra(j+1,Ne) 
...             for l in range(Nl) for j in range(Ne) for i in range(Ne) if l+1 in Lij[i][j] ],zero_matrix)
>>> H1RWA+=sum([ e*( helicity_dot_product( E_m[l],vector_element(r_p,i,j)) ) * ket(i+1,Ne)*bra(j+1,Ne) 
...             for l in range(Nl) for j in range(Ne) for i in range(Ne) if l+1 in Lij[i][j] ],zero_matrix)
    
>>> fprint(H1RWA, print_ascii=print_ascii)
[                                                                              I*t*varpi_1 ______]
[                                                                  e*r_{0;31}*e           *E_{01}]
[               0                                0                 ------------------------------]
[                                                                                2               ]
[                                                                                                ]
[                                                                              I*t*varpi_2 ______]
[                                                                  e*r_{0;32}*e           *E_{02}]
[               0                                0                 ------------------------------]
[                                                                                2               ]
[                                                                                                ]
[                   -I*t*varpi_1                     -I*t*varpi_2                                ]
[E_{01}*e*r_{0;31}*e              E_{02}*e*r_{0;32}*e                                            ]
[-------------------------------  -------------------------------                0               ]
[               2                                2                                               ]



And the complete hamiltonian in the Schrödinger picture in the rotating wave approximation is

>>> HRWA=H0+H1RWA
>>> fprint(HRWA, print_ascii=print_ascii)
[                                                                              I*t*varpi_1 ______]
[                                                                  e*r_{0;31}*e           *E_{01}]
[         hbar*omega_1                           0                 ------------------------------]
[                                                                                2               ]
[                                                                                                ]
[                                                                              I*t*varpi_2 ______]
[                                                                  e*r_{0;32}*e           *E_{02}]
[               0                          hbar*omega_2            ------------------------------]
[                                                                                2               ]
[                                                                                                ]
[                   -I*t*varpi_1                     -I*t*varpi_2                                ]
[E_{01}*e*r_{0;31}*e              E_{02}*e*r_{0;32}*e                                            ]
[-------------------------------  -------------------------------           hbar*omega_3         ]
[               2                                2                                               ]



# Rotating Frame
Next we will make a phase transformation in order to eliminate the explicit time dependance of the equations.

>>> c,ctilde,phase=define_psi_coefficients(Ne)
>>> fprint([c,ctilde,phase], print_ascii=print_ascii)
 [c1(t)]  [\tilde{c}_{1}(t)]  [theta1] 
 [     ]  [                ]  [      ] 
[[c2(t)], [\tilde{c}_{2}(t)], [theta2]]
 [     ]  [                ]  [      ] 
 [c3(t)]  [\tilde{c}_{3}(t)]  [theta3] 



>>> psi=Matrix([ exp(I*phase[i]*t)*ctilde[i] for i in range(Ne)])
>>> fprint(psi, print_ascii=print_ascii)
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
>>> fprint(lhs, print_ascii=print_ascii)
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
>>> fprint(lhs_new, print_ascii=print_ascii)
[0]
[ ]
[0]
[ ]
[0]



>>> rhs_new=Matrix([simplify(  rhs[i]*exp(-I*phase[i]*t) +hbar*phase[i]*ctilde[i] ) for i in range(Ne)])
>>> fprint(rhs_new, print_ascii=print_ascii)
[                                                                  -I*t*theta1 + I*t*theta3 + I*t*varpi_1 ______      
[                                     e*r_{0;31}*\tilde{c}_{3}(t)*e                                      *E_{01}      
[                                     -------------------------------------------------------------------------- + hba
[                                                                         2                                           
[                                                                                                                     
[                                                                  -I*t*theta2 + I*t*theta3 + I*t*varpi_2 ______      
[                                     e*r_{0;32}*\tilde{c}_{3}(t)*e                                      *E_{02}      
[                                     -------------------------------------------------------------------------- + hba
[                                                                         2                                           
[                                                                                                                     
[                                    I*t*theta1 - I*t*theta3 - I*t*varpi_1                                       I*t*t
[E_{01}*e*r_{0;31}*\tilde{c}_{1}(t)*e                                        E_{02}*e*r_{0;32}*\tilde{c}_{2}(t)*e     
[------------------------------------------------------------------------- + -----------------------------------------
[                                    2                                                                           2    
<BLANKLINE>
                                                                                               ]
                                                                                               ]
r*omega_1*\tilde{c}_{1}(t) + hbar*theta1*\tilde{c}_{1}(t)                                      ]
                                                                                               ]
                                                                                               ]
                                                                                               ]
                                                                                               ]
r*omega_2*\tilde{c}_{2}(t) + hbar*theta2*\tilde{c}_{2}(t)                                      ]
                                                                                               ]
                                                                                               ]
heta2 - I*t*theta3 - I*t*varpi_2                                                               ]
                                                                                               ]
-------------------------------- + hbar*omega_3*\tilde{c}_{3}(t) + hbar*theta3*\tilde{c}_{3}(t)]
                                                                                               ]



It can be seen that the equations loose their explicit time dependance only if $\\omega^{1} - \\theta_{1} + \\theta_{3}=0$ and $\\omega^{2} - \\theta_{2} + \\theta_{3}=0$. Which is satisfied if

>>> phase_transformation=solve([omega_laser[0]-phase[0]+phase[2],omega_laser[1]-phase[1]+phase[2]],[phase[1],phase[2]],
...                            dict=True)[0]
>>> fprint(phase_transformation, print_ascii=print_ascii)
{theta2: theta1 - varpi_1 + varpi_2, theta3: theta1 - varpi_1}



There is a free parameter $\\theta_1$, which is to be expected, since state vetors $|\\psi\\rangle$ always have a global phase invariance

>>> fprint(psi.subs(phase_transformation), print_ascii=print_ascii)
[                             I*t*theta1           ]
[           \tilde{c}_{1}(t)*e                     ]
[                                                  ]
[                  I*t*(theta1 - varpi_1 + varpi_2)]
[\tilde{c}_{2}(t)*e                                ]
[                                                  ]
[                       I*t*(theta1 - varpi_1)     ]
[     \tilde{c}_{3}(t)*e                           ]



Thus the equations become

>>> fprint(lhs_new, print_ascii=print_ascii)
[0]
[ ]
[0]
[ ]
[0]



>>> rhs_new=simplify(rhs_new.subs(phase_transformation)).expand()
>>> fprint(rhs_new, print_ascii=print_ascii)
[                                                              ______                                                 
[                                  e*r_{0;31}*\tilde{c}_{3}(t)*E_{01}                                                 
[                                  ---------------------------------- + hbar*omega_1*\tilde{c}_{1}(t) + hbar*theta1*\t
[                                                  2                                                                  
[                                                                                                                     
[                              ______                                                                                 
[  e*r_{0;32}*\tilde{c}_{3}(t)*E_{02}                                                                                 
[  ---------------------------------- + hbar*omega_2*\tilde{c}_{2}(t) + hbar*theta1*\tilde{c}_{2}(t) - hbar*varpi_1*\t
[                  2                                                                                                  
[                                                                                                                     
[E_{01}*e*r_{0;31}*\tilde{c}_{1}(t)   E_{02}*e*r_{0;32}*\tilde{c}_{2}(t)                                              
[---------------------------------- + ---------------------------------- + hbar*omega_3*\tilde{c}_{3}(t) + hbar*theta1
[                2                                    2                                                               
<BLANKLINE>
                                                 ]
                                                 ]
ilde{c}_{1}(t)                                   ]
                                                 ]
                                                 ]
                                                 ]
                                                 ]
ilde{c}_{2}(t) + hbar*varpi_2*\tilde{c}_{2}(t)   ]
                                                 ]
                                                 ]
                                                 ]
*\tilde{c}_{3}(t) - hbar*varpi_1*\tilde{c}_{3}(t)]
                                                 ]



It can be seen that this is the Schrödinger equation derived from an effective hamiltonian $\\tilde{H}$

>>> Htilde=Matrix([ [Derivative(rhs_new[i],ctilde[j]).doit() for j in range(Ne)] for i in range(Ne)])
>>> fprint(Htilde, print_ascii=print_ascii)
[                                                                                                             ______  
[                                                                                                  e*r_{0;31}*E_{01}  
[hbar*omega_1 + hbar*theta1                             0                                          -----------------  
[                                                                                                          2          
[                                                                                                                     
[                                                                                                             ______  
[                                                                                                  e*r_{0;32}*E_{02}  
[            0               hbar*omega_2 + hbar*theta1 - hbar*varpi_1 + hbar*varpi_2              -----------------  
[                                                                                                          2          
[                                                                                                                     
[    E_{01}*e*r_{0;31}                          E_{02}*e*r_{0;32}                                                     
[    -----------------                          -----------------                      hbar*omega_3 + hbar*theta1 - hb
[            2                                          2                                                             
<BLANKLINE>
          ]
          ]
          ]
          ]
          ]
          ]
          ]
          ]
          ]
          ]
          ]
ar*varpi_1]
          ]



We can see that it is convenient to choose $\\theta_1=-\\omega_1$ to simplify the hamiltonian. Also, we can recognize $\\omega^1-\\omega_3+\\omega_1=\\delta^1$ as the detuning of the first field relative to the atomic transition $\\omega_{21}=\\omega_2-\\omega_1$, and the same for $\\omega^2-\\omega_3+\\omega_2=\\delta^2$.

>>> delta1,delta2=symbols("delta1 delta2",real=True)
>>> Htilde=Htilde.subs({phase[0]:-omega_level[0]})
>>> Htilde=Htilde.subs({omega_laser[0]:delta1+omega_level[2]-omega_level[0]})
>>> Htilde=Htilde.subs({omega_laser[1]:delta2+omega_level[2]-omega_level[1]})
    
>>> Htilde=Htilde.expand()
    
>>> fprint(Htilde, print_ascii=print_ascii)
[                                                          ______]
[                                               e*r_{0;31}*E_{01}]
[        0                      0               -----------------]
[                                                       2        ]
[                                                                ]
[                                                          ______]
[                                               e*r_{0;32}*E_{02}]
[        0          -delta1*hbar + delta2*hbar  -----------------]
[                                                       2        ]
[                                                                ]
[E_{01}*e*r_{0;31}      E_{02}*e*r_{0;32}                        ]
[-----------------      -----------------         -delta1*hbar   ]
[        2                      2                                ]



If we define the Rabi frequencies $\\Omega_1 =e E_0^1 r_{0;21}/\\hbar$ and $\\Omega_2 =e E_0^2 r_{0;32}/\\hbar$

>>> Omega1,Omega2=symbols("Omega1 Omega2",real=True)
>>> Omega1,Omega2=symbols("Omega1 Omega2")
>>> Htilde=Htilde.subs({E0[0]:Omega1*hbar/r[1][2,0]/e})
>>> Htilde=Htilde.subs({E0[1]:Omega2*hbar/r[1][2,1]/e})
    
>>> fprint(Htilde, print_ascii=print_ascii)
[                                              ______ ]
[                                         hbar*Omega1 ]
[     0                   0               ----------- ]
[                                              2      ]
[                                                     ]
[                                              ______ ]
[                                         hbar*Omega2 ]
[     0       -delta1*hbar + delta2*hbar  ----------- ]
[                                              2      ]
[                                                     ]
[Omega1*hbar         Omega2*hbar                      ]
[-----------         -----------          -delta1*hbar]
[     2                   2                           ]



We define the density matrix.

>>> rho=define_density_matrix(Ne)
>>> fprint( rho , print_ascii=print_ascii)
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
>>> fprint(hamiltonian_terms, print_ascii=print_ascii)
[                                           ______                                                                    
[                  I*Omega1*rho13   I*rho31*Omega1                     I*Omega2*rho13                                 
[                  -------------- - --------------                     -------------- - I*delta1*rho12 + I*delta2*rho1
[                        2                2                                  2                                        
[                                                                                                                     
[                                                            ______                                              _____
[ I*Omega1*rho23                                     I*rho31*Omega2                     I*Omega2*rho23   I*rho32*Omega
[ -------------- + I*delta1*rho21 - I*delta2*rho21 - --------------                     -------------- - -------------
[       2                                                  2                                  2                2      
[                                                                                                                     
[                                                                                                                     
[  I*Omega1*rho11   I*Omega1*rho33   I*Omega2*rho21                     I*Omega1*rho12   I*Omega2*rho22   I*Omega2*rho
[- -------------- + -------------- - -------------- + I*delta1*rho31  - -------------- - -------------- + ------------
[        2                2                2                                  2                2                2     
<BLANKLINE>
            ______                             ______           ______           ______ ]
    I*rho32*Omega1                     I*rho11*Omega1   I*rho12*Omega2   I*rho33*Omega1 ]
2 - --------------   -I*delta1*rho13 + -------------- + -------------- - -------------- ]
          2                                  2                2                2        ]
                                                                                        ]
_                                              ______           ______           ______ ]
2                                      I*rho21*Omega1   I*rho22*Omega2   I*rho33*Omega2 ]
-                    -I*delta2*rho23 + -------------- + -------------- - -------------- ]
                                             2                2                2        ]
                                                                                        ]
                                                                 ______           ______]
33                     I*Omega1*rho13   I*Omega2*rho23   I*rho31*Omega1   I*rho32*Omega2]
-- + I*delta2*rho32  - -------------- - -------------- + -------------- + --------------]
                             2                2                2                2       ]



There are two Lindblad operators, since there are two spontaneous decay channels.

>>> lindblad_terms =gamma[2,0]*lindblad_operator(ket(1,Ne)*bra(3,Ne),rho)
>>> lindblad_terms+=gamma[2,1]*lindblad_operator(ket(2,Ne)*bra(3,Ne),rho)
    
>>> fprint(lindblad_terms, print_ascii=print_ascii)
[                                                                        gamma_31*rho13   gamma_32*rho13]
[         gamma_31*rho33                            0                  - -------------- - --------------]
[                                                                              2                2       ]
[                                                                                                       ]
[                                                                        gamma_31*rho23   gamma_32*rho23]
[                0                           gamma_32*rho33            - -------------- - --------------]
[                                                                              2                2       ]
[                                                                                                       ]
[  gamma_31*rho31   gamma_32*rho31    gamma_31*rho32   gamma_32*rho32                                   ]
[- -------------- - --------------  - -------------- - --------------  -gamma_31*rho33 - gamma_32*rho33 ]
[        2                2                 2                2                                          ]



# Optical Bloch Equations
The Optical Bloch equations are thus.

>>> eqs=hamiltonian_terms + lindblad_terms

>>> eqsign=symbols("=")
>>> eqs_list=[]
>>> for mu in range(0,Ne**2-1 -(Ne**2 - Ne)/2+1):
...     ii,jj,s=IJ(mu,Ne)
...     i=ii-1; j=jj-1
...     eqs_list+=[[Derivative(rho[i,j],t),eqsign,eqs[i,j]]]
>>> eqs_list=Matrix(eqs_list)
>>> fprint(eqs_list, print_ascii=print_ascii)
[                                                                                  ______                           ]
[d                                       I*Omega1*rho13                    I*rho31*Omega1                           ]
[--(rho11)  =                            -------------- + gamma_31*rho33 - --------------                           ]
[dt                                            2                                 2                                  ]
[                                                                                                                   ]
[                                                                                  ______                           ]
[d                                       I*Omega2*rho23                    I*rho32*Omega2                           ]
[--(rho22)  =                            -------------- + gamma_32*rho33 - --------------                           ]
[dt                                            2                                 2                                  ]
[                                                                                                                   ]
[                                                                                            ______           ______]
[d               I*Omega1*rho13   I*Omega2*rho23                                     I*rho31*Omega1   I*rho32*Omega2]
[--(rho33)  =  - -------------- - -------------- - gamma_31*rho33 - gamma_32*rho33 + -------------- + --------------]
[dt                    2                2                                                  2                2       ]
[                                                                                                                   ]
[                                                                                           ______                  ]
[d                               I*Omega1*rho23                                     I*rho31*Omega2                  ]
[--(rho21)  =                    -------------- + I*delta1*rho21 - I*delta2*rho21 - --------------                  ]
[dt                                    2                                                  2                         ]
[                                                                                                                   ]
[d               I*Omega1*rho11   I*Omega1*rho33   I*Omega2*rho21                    gamma_31*rho31   gamma_32*rho31]
[--(rho31)  =  - -------------- + -------------- - -------------- + I*delta1*rho31 - -------------- - --------------]
[dt                    2                2                2                                 2                2       ]
[                                                                                                                   ]
[d               I*Omega1*rho12   I*Omega2*rho22   I*Omega2*rho33                    gamma_31*rho32   gamma_32*rho32]
[--(rho32)  =  - -------------- - -------------- + -------------- + I*delta2*rho32 - -------------- - --------------]
[dt                    2                2                2                                 2                2       ]



which is how most literature will show the equations. However, a more convenient way to express this equations is to explicitly asume a normalized and hermitian density matrix

>>> rho=define_density_matrix(Ne,explicitly_hermitian=True,normalized=True)
>>> fprint( rho , print_ascii=print_ascii)
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

and only consider the equations for the populations $\\rho_{ii}$ for $i>1$ and the real and imaginary parts of the coherences below the diagonal.

If the density matrix is represented as a vector whose components are the these independent components of the density matrix

>>> unfolding = Unfolding(Ne, real=True, normalized=True)
>>> rho_vect = unfolding(rho)
>>> fprint(rho_vect, print_ascii=print_ascii)
[  rho22  ]
[         ]
[  rho33  ]
[         ]
[re(rho21)]
[         ]
[im(rho21)]
[         ]
[re(rho31)]
[         ]
[im(rho31)]
[         ]
[re(rho32)]
[         ]
[im(rho32)]



Then the equations can be re-written as linear combinations of these components plus an independent term.
\\begin{equation}
    \\dot{\\vec{\\rho}} = \\hat{A} \\vec{\\rho} - \\vec{b}
\\end{equation}
with $\\hat{A}$ a linear operator acting in this vector space and $\\vec{b}$ the vector of independent terms.

>>> A,b=calculate_A_b(eqs, unfolding)
>>> fprint([A, b], print_ascii=print_ascii)
 [     0          gamma_32           0                0               0             0        -im(Omega2)   re(Omega2) 
 [                                                                                                                    
 [     0         -gamma_32           0                0           im(Omega1)   -re(Omega1)    im(Omega2)   -re(Omega2)
 [                                                                                                                    
 [                                                               -im(Omega2)    re(Omega2)   -im(Omega1)   re(Omega1) 
 [     0             0               0         -delta1 + delta2  ------------   ----------   ------------  ---------- 
 [                                                                    2             2             2            2      
 [                                                                                                                    
 [                                                               -re(Omega2)   -im(Omega2)    re(Omega1)   im(Omega1) 
 [     0             0        delta1 - delta2         0          ------------  ------------   ----------   ---------- 
 [                                                                    2             2             2            2      
 [                                                                                                                    
 [-im(Omega1)                   im(Omega2)        re(Omega2)      -gamma_32                                           
[[------------  -im(Omega1)     ----------        ----------      ----------     -delta1          0             0     
 [     2                            2                 2               2                                               
 [                                                                                                                    
 [ re(Omega1)                  -re(Omega2)        im(Omega2)                    -gamma_32                             
 [ ----------    re(Omega1)    ------------       ----------        delta1      ----------        0             0     
 [     2                            2                 2                             2                                 
 [                                                                                                                    
 [ im(Omega2)   -im(Omega2)     im(Omega1)       -re(Omega1)                                  -gamma_32               
 [ ----------   ------------    ----------       ------------         0             0         ----------     -delta2  
 [     2             2              2                 2                                           2                   
 [                                                                                                                    
 [-re(Omega2)    re(Omega2)    -re(Omega1)       -im(Omega1)                                               -gamma_32  
 [------------   ----------    ------------      ------------         0             0           delta2     ---------- 
 [     2             2              2                 2                                                        2      
<BLANKLINE>
]                 
]                 
]                 
]                 
]  [     0      ] 
]  [            ] 
]  [     0      ] 
]  [            ] 
]  [     0      ] 
]  [            ] 
]  [     0      ] 
]  [            ] 
]  [-im(Omega1) ] 
], [------------]]
]  [     2      ] 
]  [            ] 
]  [ re(Omega1) ] 
]  [ ---------- ] 
]  [     2      ] 
]  [            ] 
]  [     0      ] 
]  [            ] 
]  [     0      ] 
]                 
]                 
]                 
]                 



Explicitly, this is

>>> eqs_new=A*rho_vect - b
    
>>> eqs_list=[]
>>> for mu in range(0,Ne**2-1 ):
...     eqs_list+=[[Derivative(rho_vect[mu],t),eqsign,eqs_new[mu]]]
>>> eqs_list=Matrix(eqs_list)
>>> fprint(eqs_list, print_ascii=print_ascii)
[  d                                                                                                                  
[  --(rho22)    =                                       gamma_32*rho33 + re(Omega2)*im(rho32) - re(rho32)*im(Omega2)  
[  dt                                                                                                                 
[                                                                                                                     
[  d                                                                                                                  
[  --(rho33)    =                -gamma_32*rho33 - re(Omega1)*im(rho31) - re(Omega2)*im(rho32) + re(rho31)*im(Omega1) 
[  dt                                                                                                                 
[                                                                                                                     
[d                                                       re(Omega1)*im(rho32)   re(Omega2)*im(rho31)   re(rho31)*im(Om
[--(re(rho21))  =         (-delta1 + delta2)*im(rho21) + -------------------- + -------------------- - ---------------
[dt                                                               2                      2                      2     
[                                                                                                                     
[d                                                       re(Omega1)*re(rho32)   re(Omega2)*re(rho31)   im(Omega1)*im(r
[--(im(rho21))  =          (delta1 - delta2)*re(rho21) + -------------------- - -------------------- + ---------------
[dt                                                               2                      2                      2     
[                                                                                                                     
[d                                     gamma_32*re(rho31)   rho22*im(Omega1)                      re(Omega2)*im(rho21)
[--(re(rho31))  =  -delta1*im(rho31) - ------------------ - ---------------- - rho33*im(Omega1) + --------------------
[dt                                            2                   2                                       2          
[                                                                                                                     
[d                                    gamma_32*im(rho31)   rho22*re(Omega1)                      re(Omega1)   re(Omega
[--(im(rho31))  =  delta1*re(rho31) - ------------------ + ---------------- + rho33*re(Omega1) - ---------- - --------
[dt                                           2                   2                                  2                
[                                                                                                                     
[d                                           gamma_32*re(rho32)   rho22*im(Omega2)   rho33*im(Omega2)   re(Omega1)*im(
[--(re(rho32))  =        -delta2*im(rho32) - ------------------ + ---------------- - ---------------- - --------------
[dt                                                  2                   2                  2                    2    
[                                                                                                                     
[d                                           gamma_32*im(rho32)   rho22*re(Omega2)   rho33*re(Omega2)   re(Omega1)*re(
[--(im(rho32))  =         delta2*re(rho32) - ------------------ - ---------------- + ---------------- - --------------
[dt                                                  2                   2                  2                    2    
<BLANKLINE>
                                    ]
                                    ]
                                    ]
                                    ]
                                    ]
+ re(rho32)*im(Omega2)              ]
                                    ]
                                    ]
ega2)   re(rho32)*im(Omega1)        ]
----- - --------------------        ]
                 2                  ]
                                    ]
ho32)   im(Omega2)*im(rho31)        ]
----- - --------------------        ]
                 2                  ]
                                    ]
   re(rho21)*im(Omega2)   im(Omega1)]
 + -------------------- + ----------]
            2                 2     ]
                                    ]
2)*re(rho21)   im(Omega2)*im(rho21) ]
------------ + -------------------- ]
 2                      2           ]
                                    ]
rho21)   re(rho21)*im(Omega1)       ]
------ + --------------------       ]
                  2                 ]
                                    ]
rho21)   im(Omega1)*im(rho21)       ]
------ - --------------------       ]
                  2                 ]



Which is the same as the equations in the previous form.

>>> ss_comp={ rho[i,j]:re(rho[i,j])+I*im(rho[i,j]) for j in range(Ne) for i in range(Ne)}
    
>>> test= (eqs_new - Matrix([re(eqs[1,1]),re(eqs[2,2]),
...                           re(eqs[1,0]),re(eqs[2,0]),re(eqs[2,1]),
...                           im(eqs[1,0]),im(eqs[2,0]),im(eqs[2,1])])).expand()
    
>>> test=test.subs(ss_comp)
>>> fprint(test,print_ascii=print_ascii)
[                                                                                                                     
[                                                                                                                     
[                                                                                                                     
[                                                                                                                     
[                                                                                                                     
[                                                                                                                     
[                                                         gamma_32*re(rho31)   rho22*im(Omega1)                      r
[delta1*re(rho21) + delta1*im(rho31) - delta2*re(rho21) + ------------------ + ---------------- + rho33*im(Omega1) + -
[                                                                 2                   2                               
[                                                                                                                     
[                                          gamma_32*re(rho31)   gamma_32*re(rho32)   rho22*im(Omega1)   rho22*im(Omega
[   -delta1*im(rho31) + delta2*im(rho32) - ------------------ + ------------------ - ---------------- - --------------
[                                                  2                    2                   2                  2      
[                                                                                                                     
[                                                          gamma_32*im(rho31)   rho22*re(Omega1)                      
[-delta1*re(rho21) + delta1*re(rho31) + delta2*re(rho21) - ------------------ + ---------------- + rho33*re(Omega1) - 
[                                                                  2                   2                              
[                                                                                                                     
[                                          gamma_32*re(rho32)   gamma_32*im(rho31)   rho22*re(Omega1)   rho22*im(Omega
[   -delta1*re(rho31) - delta2*im(rho32) - ------------------ + ------------------ - ---------------- + --------------
[                                                  2                    2                   2                  2      
[                                                                                                                     
[                                                                                                                     
<BLANKLINE>
               0                                                                                                      
<BLANKLINE>
               0                                                                                                      
<BLANKLINE>
               0                                                                                                      
<BLANKLINE>
e(Omega1)*re(rho32)   re(Omega2)*re(rho31)   re(Omega2)*im(rho21)   re(rho21)*im(Omega2)   im(Omega1)*im(rho32)   im(O
------------------- - -------------------- - -------------------- - -------------------- + -------------------- - ----
        2                      2                      2                      2                      2                 
<BLANKLINE>
2)                      rho33*im(Omega2)   re(Omega1)*im(rho21)   re(Omega2)*im(rho21)   re(rho21)*im(Omega1)   re(rho
-- - rho33*im(Omega1) + ---------------- + -------------------- + -------------------- - -------------------- + ------
                               2                    2                      2                      2                   
<BLANKLINE>
re(Omega1)*re(rho32)   re(Omega1)   re(Omega2)*re(rho21)   re(Omega2)*re(rho31)   im(Omega1)*im(rho32)   im(Omega2)*im
-------------------- - ---------- - -------------------- + -------------------- - -------------------- + -------------
         2                 2                 2                      2                      2                      2   
<BLANKLINE>
2)                      rho33*im(Omega2)   re(Omega1)*im(rho21)   re(Omega1)   re(Omega2)*re(rho21)   re(rho21)*im(Ome
-- - rho33*re(Omega1) - ---------------- - -------------------- + ---------- + -------------------- + ----------------
                               2                    2                 2                 2                      2      
<BLANKLINE>
               0                                                                                                      
<BLANKLINE>
                              ]
                              ]
                              ]
                              ]
                              ]
                              ]
mega1)   im(Omega2)*im(rho31) ]
------ - -------------------- ]
2                 2           ]
                              ]
21)*im(Omega2)   im(Omega1)   ]
-------------- + ----------   ]
   2                 2        ]
                              ]
(rho21)   im(Omega2)*im(rho31)]
------- + --------------------]
                   2          ]
                              ]
ga1)   im(Omega2)*im(rho21)   ]
---- - --------------------   ]
                2             ]
                              ]
                              ]



[1]  H.J. Metcalf and P. van der Straten. Laser Cooling and Trapping. Graduate Texts in Contempo-
rary Physics. Springer New York, 2001.

[2] Daniel Adam Steck. Quantum and Atom Optics. Oregon Center for Optics and Department of Physics, University of Oregon Copyright © 200

[]

"""
__doc__=__doc__.replace("+IGNORE_PLOT_STEP1", "+ELLIPSIS\n[<...>]")
__doc__=__doc__.replace("+IGNORE_PLOT_STEP2", "+ELLIPSIS\n<...>")
__doc__=__doc__.replace("+IGNORE_PLOT_STEP3", "+ELLIPSIS\n(...)")
__doc__=__doc__.replace("+IGNORE_PLOT_STEP5", "+ELLIPSIS\nText(...)")
__doc__=__doc__.replace("+IGNORE_PLOT_STEP4", "\n")
