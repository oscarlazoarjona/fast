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


>>> plots_path="folder_09___Thermal_States"

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
<matplotlib.figure.Figure at 0x7ff6d4ee6fd0>



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
<matplotlib.figure.Figure at 0x7ff6d4ee6a50>



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
<matplotlib.figure.Figure at 0x7ff6d4ef5f10>



[]

"""
__doc__=__doc__.replace("+IGNORE_PLOT_STEP1", "+ELLIPSIS\n[<...>]")
__doc__=__doc__.replace("+IGNORE_PLOT_STEP2", "+ELLIPSIS\n<...>")
__doc__=__doc__.replace("+IGNORE_PLOT_STEP3", "+ELLIPSIS\n(...)")
__doc__=__doc__.replace("+IGNORE_PLOT_STEP4", "\n")
