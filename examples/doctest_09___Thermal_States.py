# -*- coding: utf-8 -*-
# Copyright (C) 2017 Oscar Gerardo Lazo Arjona
# mailto: oscar.lazoarjona@physics.ox.ac.uk

__doc__ = r"""

Let us calculate the populations $p_i$ for the hyperfine magnetic states of the $5S_{1/2}$ ground state of rubidium in a thermal state. A thermal state has populations

$$p_i=\\frac{ \\exp{(-E_i/k_B T)} }{Z} \\hspace{1cm} Z=\\sum_j \\exp{(-E_j/k_B T)} $$

where $E_i$ are the energies of each state and $k_B$ is Boltzman's constant.

>>> from fast.atomic_structure import State, make_list_of_states
>>> from sympy import Integer
>>> from math import pi,exp,log

>>> element="Rb"; iso=85
>>> ground=State(element,iso,5,0,1/Integer(2))
>>> ground
85Rb 5S_1/2



>>> magnetic_states=make_list_of_states([ground],"magnetic")
>>> for state in magnetic_states:
...     print state
... 
85Rb 5S_1/2^2,-2
85Rb 5S_1/2^2,-1
85Rb 5S_1/2^2,0
85Rb 5S_1/2^2,1
85Rb 5S_1/2^2,2
85Rb 5S_1/2^3,-3
85Rb 5S_1/2^3,-2
85Rb 5S_1/2^3,-1
85Rb 5S_1/2^3,0
85Rb 5S_1/2^3,1
85Rb 5S_1/2^3,2
85Rb 5S_1/2^3,3



>>> e1=magnetic_states[0]

>>> from scipy.constants import physical_constants
>>> hbar=physical_constants["Planck constant over 2 pi"][0] # J s
>>> kB  =physical_constants["Boltzmann constant"][0] #        J/K

>>> E=[hbar*2*pi*e.nu for e in magnetic_states]
>>> for ii in E:
...     print ii
... 
-1.17337391713e-24
-1.17337391713e-24
-1.17337391713e-24
-1.17337391713e-24
-1.17337391713e-24
8.38124226523e-25
8.38124226523e-25
8.38124226523e-25
8.38124226523e-25
8.38124226523e-25
8.38124226523e-25
8.38124226523e-25



>>> T=20#Celsius degrees.
>>> T=T+237.15 #Kelvin
>>> print T
257.15



>>> Z=sum([exp(-E[i]/kB/T) for i in range(len(magnetic_states))])
>>> print Z
12.0000004681



So the populations are:

>>> p=[exp(-E[i]/kB/T)/Z for i in range(len(magnetic_states))]
>>> for ii in p: print ii
0.083360876002
0.083360876002
0.083360876002
0.083360876002
0.083360876002
0.0833136599986
0.0833136599986
0.0833136599986
0.0833136599986
0.0833136599986
0.0833136599986
0.0833136599986



Which is close to assigning populations equally:

>>> p_deg=[1.0/len(magnetic_states) for i in range(len(magnetic_states))]
>>> for ii in p_deg: print ii
0.0833333333333
0.0833333333333
0.0833333333333
0.0833333333333
0.0833333333333
0.0833333333333
0.0833333333333
0.0833333333333
0.0833333333333
0.0833333333333
0.0833333333333
0.0833333333333



Let's make a plot of the variation with temperature.

>>> if iso==85:
...     a=5
... else:
...     a=3
... 

>>> P1=sum(p[:a])
>>> P2=sum(p[a:])
>>> print P1,P2
0.41680438001 0.58319561999



>>> def get_energies(element,iso):
...     if element=="Rb":
...         ground=State(element,iso,5,0,1/Integer(2))
...     else:
...         ground=State(element,iso,6,0,1/Integer(2))
...     magnetic_states=make_list_of_states([ground],"magnetic")
...     E=[hbar*2*pi*e.nu for e in magnetic_states]
...     return E
... 

>>> def pops(T,E,element,iso):
...     if iso ==133:
...         a=7
...     elif iso ==85:
...         a=5
...     else:
...         a=3
...     
...     if element=="Rb":
...         ground=State(element,iso,5,0,1/Integer(2))
...     else:
...         ground=State(element,iso,6,0,1/Integer(2))
...     magnetic_states=make_list_of_states([ground],"magnetic")
...     
...     Z=sum([exp(-E[i]/kB/T) for i in range(len(magnetic_states))])
...     p=[exp(-E[i]/kB/T)/Z for i in range(len(magnetic_states))]
...     return p[0],p[-1],sum(p[:a]),sum(p[a:])
... 

>>> from numpy import logspace,array
>>> T=logspace(-2,3,201)

>>> E85 =get_energies("Rb",85)
>>> E87 =get_energies("Rb",87)
>>> E133=get_energies("Cs",133)

>>> dat=array([pops(Ti,E85,"Rb",85) for Ti in T])
>>> p185=list(dat[:,0])
>>> p285=list(dat[:,1])
>>> P185=list(dat[:,2])
>>> P285=list(dat[:,3])

>>> dat=array([pops(Ti,E87,"Rb",87) for Ti in T])
>>> p187=list(dat[:,0])
>>> p287=list(dat[:,1])
>>> P187=list(dat[:,2])
>>> P287=list(dat[:,3])

>>> dat=array([pops(Ti,E133,"Cs",133) for Ti in T])
>>> p1133=list(dat[:,0])
>>> p2133=list(dat[:,1])
>>> P1133=list(dat[:,2])
>>> P2133=list(dat[:,3])

>>> show_mot_temperature=False
>>> if show_mot_temperature:
...     T_Doppler=145.537e-6
...     T_lim=100e-6
...     T=[T_lim]+list(T)
...     T=array(T)
...     
...     p185=[1/5.0]+p185
...     p285=[0.0  ]+p285
...     P185=[1.0  ]+P185
...     P285=[0.0  ]+P285
...     
...     p187=[1/3.0]+p187
...     p287=[0.0  ]+p287
...     P187=[1.0  ]+P187
...     P287=[0.0  ]+P287
... 

>>> from matplotlib import pyplot


>>> plots_path="folder_09___Thermal_States/" 

>>> pyplot.close("all")
>>> pyplot.semilogx(T,p185,"b",label=r"$\mathrm{lower \ magnetic \ state}$") # doctest: +IGNORE_PLOT_STEP1
>>> pyplot.semilogx(T,p285,"r",label=r"$\mathrm{higher \ magnetic \ state}$") # doctest: +IGNORE_PLOT_STEP1
>>> pyplot.semilogx(T,P185,"m",label=r"$\mathrm{lower  \ multiplet}$") # doctest: +IGNORE_PLOT_STEP1
>>> pyplot.semilogx(T,P285,"g",label=r"$\mathrm{higher \ multiplet}$") # doctest: +IGNORE_PLOT_STEP1
    
>>> pyplot.semilogx(T,p187,"b--") # doctest: +IGNORE_PLOT_STEP1
>>> pyplot.semilogx(T,p287,"r--") # doctest: +IGNORE_PLOT_STEP1
>>> pyplot.semilogx(T,P187,"m--") # doctest: +IGNORE_PLOT_STEP1
>>> pyplot.semilogx(T,P287,"g--") # doctest: +IGNORE_PLOT_STEP1
    
>>> pyplot.semilogx(T,p1133,"b:") # doctest: +IGNORE_PLOT_STEP1
>>> pyplot.semilogx(T,p2133,"r:") # doctest: +IGNORE_PLOT_STEP1
>>> pyplot.semilogx(T,P1133,"m:") # doctest: +IGNORE_PLOT_STEP1
>>> pyplot.semilogx(T,P2133,"g:") # doctest: +IGNORE_PLOT_STEP1
    
    
>>> pyplot.semilogx([273.15,273.15],[0,1],"k") # doctest: +IGNORE_PLOT_STEP1
>>> pyplot.semilogx([273.15+50,273.15+50],[0,1],"k") # doctest: +IGNORE_PLOT_STEP1
    
>>> pyplot.ylabel(r"$\mathrm{population}$",fontsize=15) # doctest: +IGNORE_PLOT_STEP2
>>> pyplot.xlabel(r"$T \ \mathrm{(K)}$",fontsize=15) # doctest: +IGNORE_PLOT_STEP2
>>> pyplot.legend(fontsize=10) # doctest: +IGNORE_PLOT_STEP2
    
>>> pyplot.ylim([0,1]) # doctest: +IGNORE_PLOT_STEP3
>>> pyplot.savefig(plots_path+"/01_populations.png",bbox_inches="tight") # doctest: +IGNORE_PLOT_STEP4
<matplotlib.figure.Figure at 0x7fa52802c850>



Here solid lines show the populations for $^{85}\\mathrm{Rb}$, dashed lines for $^{87}\\mathrm{Rb}$, and dotted lines for $^{133}\\mathrm{Cs}$. We can see that $^{85}\\mathrm{Rb}$ is thermalized at slightly lower temperatures than $^{87}\\mathrm{Rb}$: at about 5 and 10 Kelvin respectively all magnetic states reach the same population.

>>> def entropy(p1,p2,iso):
...     if iso==133:
...         n1=7; n2=9
...     elif iso==85:
...         n1=5; n2=7
...     elif iso==87:
...         n1=3; n2=5
...     return - kB*(n1*p1*log(p1) + n2*p2*log(p2))
... 

>>> S85 =[entropy(p185[i], p285[i], 85) for i in range(len(p185))]
>>> S87 =[entropy(p187[i], p287[i], 87) for i in range(len(p187))]
>>> S133=[entropy(p1133[i],p2133[i],133) for i in range(len(p187))]

>>> pyplot.close("all")
>>> pyplot.semilogx(T,S85, "r",label=r"$^{85}  \mathrm{Rb}$") # doctest: +IGNORE_PLOT_STEP1
>>> pyplot.semilogx(T,S87, "b",label=r"$^{87}  \mathrm{Rb}$") # doctest: +IGNORE_PLOT_STEP1
>>> pyplot.semilogx(T,S133,"g",label=r"$^{133} \mathrm{Cs}$") # doctest: +IGNORE_PLOT_STEP1
    
>>> pyplot.semilogx([273.15,273.15]      ,[1.5e-23,4e-23],"k") # doctest: +IGNORE_PLOT_STEP1
>>> pyplot.semilogx([273.15+50,273.15+50],[1.5e-23,4e-23],"k") # doctest: +IGNORE_PLOT_STEP1
    
>>> pyplot.ylabel(r"$S \ \mathrm{(J K^{-1})}$",fontsize=15) # doctest: +IGNORE_PLOT_STEP2
>>> pyplot.xlabel(r"$T \ \mathrm{(K)}$",fontsize=15) # doctest: +IGNORE_PLOT_STEP2
>>> pyplot.legend(fontsize=15,loc="lower center") # doctest: +IGNORE_PLOT_STEP2
>>> pyplot.savefig(plots_path+"/02_entropy.png",bbox_inches="tight") # doctest: +IGNORE_PLOT_STEP4
<matplotlib.figure.Figure at 0x7fa5272ff290>



>>> def average_energy(p1,p2,iso):
...     if iso==133:
...         n1=7; n2=9; element="Cs" 
...     elif iso==85:
...         n1=5; n2=7; element="Rb" 
...     elif iso==87:
...         n1=3; n2=5; element="Rb" 
...     
...     E=get_energies(element,iso)
...     E1=E[0]; E2=E[-1]
...     return (n1*p1*E1 + n2*p2*E2)*1e-9/hbar/2/pi
... 

>>> E85 =[average_energy(p185[i], p285[i], 85) for i in range(len(p185))]
>>> E87 =[average_energy(p187[i], p287[i], 87) for i in range(len(p187))]
>>> E133=[average_energy(p1133[i],p2133[i],133) for i in range(len(p187))]

>>> pyplot.close("all")
>>> pyplot.semilogx(T,E85, "r",label=r"$^{85}  \mathrm{Rb}$") # doctest: +IGNORE_PLOT_STEP1
>>> pyplot.semilogx(T,E87, "b",label=r"$^{87}  \mathrm{Rb}$") # doctest: +IGNORE_PLOT_STEP1
>>> pyplot.semilogx(T,E133,"g",label=r"$^{133} \mathrm{Cs}$") # doctest: +IGNORE_PLOT_STEP1
    
>>> pyplot.semilogx([273.15,273.15]      ,[-5.1,0],"k") # doctest: +IGNORE_PLOT_STEP1
>>> pyplot.semilogx([273.15+50,273.15+50],[-5.1,0],"k") # doctest: +IGNORE_PLOT_STEP1
    
>>> pyplot.ylabel(r"$E \ \mathrm{(GHz)}$",fontsize=15) # doctest: +IGNORE_PLOT_STEP2
>>> pyplot.xlabel(r"$T \ \mathrm{(K)}$",fontsize=15) # doctest: +IGNORE_PLOT_STEP2
>>> pyplot.legend(fontsize=15,loc="lower center") # doctest: +IGNORE_PLOT_STEP2
    
>>> pyplot.savefig(plots_path+"/03_energy.png",bbox_inches="tight") # doctest: +IGNORE_PLOT_STEP4
<matplotlib.figure.Figure at 0x7fa527144f50>



# An atom in thermal equilibrium
Let us solve the master equation found in [1]. These are simply the usual Bloch equations for a two level system.

>>> from fast import *
>>> from sympy import oo,exp,I
>>> init_printing()
>>> use_unicode=True
>>> use_unicode=False

>>> hbar,k,omega2,T,gamma,n=symbols("hbar k omega2 T gamma n",positive=True)
>>> omega1=symbols("omega1",negative=True)

>>> Omega=symbols("Omega",real=True)
>>> delta=symbols("delta",real=True)

>>> H=hbar*Matrix([[delta,Omega.conjugate()/2],[Omega/2,0]])

>>> rho=define_density_matrix(2,explicitly_hermitian=True,normalized=True)

>>> eqs =I/hbar*(rho*H-H*rho)

However, instead of the usual Lindblad terms, we will use the following

>>> eqs+=(1+n)*gamma*lindblad_operator(ketbra(1,2,2),rho)
>>> eqs+=(  n)*gamma*lindblad_operator(ketbra(2,1,2),rho)

Where $\\bar{n}$ is the "mean photon number" here simply defined as 

>>> ns=exp(-hbar*(omega2-omega1)/k/T)/(1-exp(-hbar*(omega2-omega1)/k/T))
>>> ns=1/(exp(hbar*(omega2-omega1)/k/T)-1)
>>> pprint(ns,use_unicode=use_unicode)
             1              
----------------------------
 hbar*(-omega1 + omega2)    
 -----------------------    
           T*k              
e                        - 1



Where $\\omega_1$ and $\\omega_2$ are the energy frequencies of our two states, and $T$ will become our definition of temperature. We solve these equations for the steady state. Notice that $n\\rightarrow0$ as $T\\rightarrow0$  and $n \\rightarrow\\infty$ as $T\\rightarrow \\infty$.

>>> print ns.limit(T,0), ns.limit(T,oo)
0 oo



>>> eq2=eqs[1,1].subs({rho[1,0]:re(rho[1,0])+I*im(rho[1,0])}).expand()
>>> eq3=re(eqs[1,0].expand())
>>> eq4=im(eqs[0,1].expand())

>>> sol=solve([eq2,eq3,eq4],[rho[1,1],re(rho[1,0]),im(rho[1,0])])

>>> rho11T=(1-sol[rho[1,1]]).expand().factor()
>>> rho22T=sol[rho[1,1]]
>>> rho21T=(sol[re(rho[1,0])]+I*sol[im(rho[1,0])]).factor()

>>> pprint(rho11T,use_unicode=use_unicode,num_columns=150)
       2          2          2            2          2  3          2  2          2          2
2*Omega *n + Omega  + 4*delta *n + 4*delta  + 4*gamma *n  + 8*gamma *n  + 5*gamma *n + gamma 
---------------------------------------------------------------------------------------------
                       /       2          2          2  2          2          2\             
             (2*n + 1)*\2*Omega  + 4*delta  + 4*gamma *n  + 4*gamma *n + gamma /             



>>> pprint(rho22T,use_unicode=use_unicode,num_columns=150)
     2               /       2        2          2\
Omega *(2*n + 1) + n*\4*delta  + gamma *(2*n + 1) /
---------------------------------------------------
          /       2          2        2          2\
(2*n + 1)*\2*Omega  + 4*delta  + gamma *(2*n + 1) /



>>> pprint(rho21T,use_unicode=use_unicode,num_columns=150)
             -Omega*(-2*delta + 2*I*gamma*n + I*gamma)             
-------------------------------------------------------------------
          /       2          2          2  2          2          2\
(2*n + 1)*\2*Omega  + 4*delta  + 4*gamma *n  + 4*gamma *n + gamma /



Obviously, if we take the temperature to zero, we recover the usual solutions to Bloch equations in the steady state.

>>> pprint(rho11T.limit(n,0),use_unicode=use_unicode)
      2          2        2 
 Omega  + 4*delta  + gamma  
----------------------------
       2          2        2
2*Omega  + 4*delta  + gamma 



>>> pprint(rho22T.limit(n,0),use_unicode=use_unicode)
                2           
           Omega            
----------------------------
       2          2        2
2*Omega  + 4*delta  + gamma 



>>> pprint(rho21T.limit(n,0),use_unicode=use_unicode)
-(-2*Omega*delta + I*Omega*gamma) 
----------------------------------
          2          2        2   
   2*Omega  + 4*delta  + gamma    



Now let's explore what happens when the temperature goes to infinity:

>>> print [rho11T.limit(n,oo), rho22T.limit(n,oo), rho21T.limit(n,oo)]
[1/2, 1/2, 0]



Which is exactly what one would expect of thermal states. Furthermore, at finite temperature, and in the abscence of optical fields

>>> vapour={delta:0,Omega:0}
>>> rho11_vapour=rho11T.subs(vapour).factor()
>>> rho22_vapour=rho22T.subs(vapour).factor()
>>> rho21_vapour=rho21T.subs(vapour)
>>> pprint([rho11_vapour,rho22_vapour,rho21_vapour],use_unicode=use_unicode)
  n + 1      n       
[-------, -------, 0]
 2*n + 1  2*n + 1    



Remarkably, the populations do not depend on the decay frequency $\\gamma$! Thus only the temperature determines the populations in the abscence of fields!

Explicitly, the populations are

>>> rho11_vapour=rho11_vapour.subs({n:ns}).expand().factor()
>>> rho22_vapour=rho22_vapour.subs({n:ns}).expand().factor()
>>> pprint([rho11_vapour,rho22_vapour],use_unicode=use_unicode)
          hbar*omega2                  hbar*omega1        
          -----------                  -----------        
              T*k                          T*k            
         e                            e                   
[---------------------------, ---------------------------]
  hbar*omega1    hbar*omega2   hbar*omega1    hbar*omega2 
  -----------    -----------   -----------    ----------- 
      T*k            T*k           T*k            T*k     
 e            + e             e            + e            



We can check that these are in deed thermal states

>>> Z=exp(-hbar*omega1/k/T)+exp(-hbar*omega2/k/T)
>>> rho11_thermal=exp(-hbar*omega1/k/T)/Z
>>> rho22_thermal=exp(-hbar*omega2/k/T)/Z
    
>>> print [(rho11_vapour-rho11_thermal).simplify(),(rho22_vapour-rho22_thermal).simplify()]
[0, 0]



So we may define a thermalization temperature as the temperature required so that in the abscence of fields, the excited state has population 1/4 (in much the same way as the saturation intensity is defined for the field).

>>> Tterm=solve(rho22_vapour-1/Integer(4),T)[0]
>>> pprint(Tterm,use_unicode=use_unicode)
-hbar*(omega1 - omega2) 
------------------------
        k*log(3)        



Let us now apply this to the hyperfine splittings of the ground states of the alkalis that we examined before, which is very questionable.

>>> g1Rb85=State("Rb",85,5,0,1/Integer(2),2)
>>> g2Rb85=State("Rb",85,5,0,1/Integer(2),3)
    
>>> g1Rb87=State("Rb",87,5,0,1/Integer(2),1)
>>> g2Rb87=State("Rb",87,5,0,1/Integer(2),2)
    
>>> g1Cs133=State("Cs",133,6,0,1/Integer(2),3)
>>> g2Cs133=State("Cs",133,6,0,1/Integer(2),4)

>>> from scipy.constants import k as ks
>>> from scipy.constants import hbar as hbars

>>> ns.subs({hbar:hbars,k:ks,omega1:g1Cs133.omega,omega2:g2Cs133.omega,T:293.})
663.632828136148



>>> TtermRb85 =Tterm.subs({hbar:hbars,k:ks,omega1: g1Rb85.omega, omega2:g2Rb85.omega}).n()
>>> TtermRb87 =Tterm.subs({hbar:hbars,k:ks,omega1: g1Rb87.omega, omega2:g2Rb87.omega}).n()
>>> TtermCs133=Tterm.subs({hbar:hbars,k:ks,omega1:g1Cs133.omega,omega2:g2Cs133.omega}).n()
    
>>> print [TtermRb85, TtermRb87, TtermCs133]
[0.132614817483323, 0.298570427225583, 0.401576510739131]



Let us now find the populations for the exited multiplets for a thermal state at these temperatures:

>>> E85 =get_energies("Rb",85)
>>> E87 =get_energies("Rb",87)
>>> E133=get_energies("Cs",133)
>>> print pops(TtermRb85,E85,"Rb",85)[3].subs({hbar:hbars})
0.318181818181818



>>> print pops(TtermRb87,E87,"Rb",87)[3].subs({hbar:hbars})
0.357142857142857



>>> print pops(TtermCs133,E133,"Cs",133)[3].subs({hbar:hbars})
0.300000000000000



If the two level system theory was appropiate for this problem, we should have got 1/4, so we did get pretty good estimates!

>>> pyplot.close("all")

[]

[1] An open systems approach to quantum optics : lectures presented at the Universite Libre de Bruxelles, October 28 to November 4, 1991. Carmichael, Howard.

[]

"""
__doc__=__doc__.replace("+IGNORE_PLOT_STEP1", "+ELLIPSIS\n[<...>]")
__doc__=__doc__.replace("+IGNORE_PLOT_STEP2", "+ELLIPSIS\n<...>")
__doc__=__doc__.replace("+IGNORE_PLOT_STEP3", "+ELLIPSIS\n(...)")
__doc__=__doc__.replace("+IGNORE_PLOT_STEP4", "\n")
