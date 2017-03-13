# -*- coding: utf-8 -*-
# Copyright (C) 2017 Oscar Gerardo Lazo Arjona
# mailto: oscar.lazoarjona@physics.ox.ac.uk

__doc__ = r"""

### Two level atom
We import all the functions of FAST and some other useful stuff.

>>> from fast import *
>>> from math import pi,sqrt

>>> from matplotlib import pyplot
>>> from fast.config import parallel, use_netcdf, fast_path
>>> from numpy import array

We establish the basic characteristics of te experiment. The path where we will work, and the name of the experiment.

>>> path=fast_path[:-5]+"/examples/folder_01___Two_level_atom/" 
>>> name='suite'

The number of states.

>>> Ne=2

The properties of the atom. The frequencies $\\omega_i$ of the energy levels $E_i=\\hbar \\omega_i$.

>>> omega_states=[0.0,100.0]

And from these the matrix of transition frequencies is calculated $\\omega_{ij}=\\omega_i-\\omega_j$.

>>> omega=[[omega_states[i]-omega_states[j] for j in range(Ne)] for i in range(Ne)]
>>> print array(omega)
[[   0. -100.]
 [ 100.    0.]]



The matrix $\\gamma_{ij}$ of decay frequencies for transitions $|i\\rangle \\rightarrow |j\\rangle$.

>>> gamma=1.0
>>> gamma=[[0.,-gamma],[gamma,0.]]
>>> print array(gamma)
[[ 0. -1.]
 [ 1.  0.]]



The matrix form of the position operator $\\hat{\\vec{r}}$ in the helicity basis (in Bohr radii).

>>> r=[ [[0,1],[1,0]] for p in range(3)]
>>> print array(r)
[[[0 1]
  [1 0]]
<BLANKLINE>
 [[0 1]
  [1 0]]
<BLANKLINE>
 [[0 1]
  [1 0]]]



We define the lasers.

>>> l1=PlaneWave(0,pi/2,0,0,color="red")
>>> laseres=[l1]
>>> Nl=len(laseres)
    
>>> fig = pyplot.figure(); ax = fig.gca(projection='3d')
>>> draw_lasers_3d(ax,laseres,path+'lasers.png') # doctest: +IGNORE_PLOT_STEP4
<matplotlib.figure.Figure at 0x7fa6c30e1a50>



We specify the coupling of the lasers. In this case, the transition $|1\\rangle \\rightarrow |2\\rangle$ is coupled by laser $1$.

>>> Lij=[[1,2,[1]]]
>>> Lij=formatLij(Lij,Ne)

We make a diagram level.

>>> fig=pyplot.figure(); ax=fig.add_subplot(111,aspect="equal")
    
>>> p1=[0.5,1]; p2=[1.5,3]
>>> draw_state(ax,p1,text=r"$|1\rangle$",l=1.0,alignment='right',label_displacement=0.05,fontsize=25,linewidth=4.0)
>>> draw_state(ax,p2,text=r"$|2\rangle$",l=1.0,alignment='right',label_displacement=0.05,fontsize=25,linewidth=4.0)
    
>>> excitation(ax,[p1[0]+0.25,p1[1]],[p2[0]+0.25,p2[1]], fc="red", ec="red",width=0.01, head_width=0.1, head_length=0.1) # doctest: +IGNORE_PLOT_STEP4
>>> decay(     ax,[p1[0]-0.25,p1[1]],[p2[0]-0.25,p2[1]], 0.05,10.0,color="red",linewidth=1.0) # doctest: +IGNORE_PLOT_STEP4
    
>>> pyplot.axis('off') # doctest: +IGNORE_PLOT_STEP3
>>> pyplot.savefig(path+name+'_diagram.png',bbox_inches="tight") # doctest: +IGNORE_PLOT_STEP4
<matplotlib.figure.Figure at 0x7fa6bfb14150>



## Time evolution
We write the Fortran code of the experiment.

>>> tw=write_evolution(path,name+"_evolution",laseres,omega,gamma,r,Lij)

We compile it.

>>> tc=compile_code(path,name+"_evolution",lapack=True,parallel=parallel)

We specify the intensity of the laser (as the electric field amplitude). And the detuning.

>>> delta1=-1.0; E00=1.0

We run the time evolution with these parameters.

>>> tr=run_evolution(path,name+"_evolution",[E00],[delta1],  2000,  0.01,Ne,use_netcdf=use_netcdf)

We read the results.

>>> t,rho22,rho21_real,rho21_imag = read_result(path,name+"_evolution",N=Ne,use_netcdf=use_netcdf)

We plot the components we have just read.

>>> fig=pyplot.figure(); ax=fig.add_subplot(111)
>>> ax.plot(t,rho22,'k-'      ,label=r"$\rho_{22}$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(t,rho21_real,'b-' ,label=r"$\mathfrak{Re}\rho_{21}$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(t,rho21_imag,'r-' ,label=r"$\mathfrak{Im}\rho_{21}$") # doctest: +IGNORE_PLOT_STEP1
    
>>> ax.set_ylim([None,0.6]) # doctest: +IGNORE_PLOT_STEP3
>>> ax.legend(loc=0,fontsize=15) # doctest: +IGNORE_PLOT_STEP2
>>> ax.set_xlabel(r"$t$",fontsize=20) # doctest: +IGNORE_PLOT_STEP2
    
>>> pyplot.savefig(path+'evolution_evo.png',bbox_inches='tight',figsize=25) # doctest: +IGNORE_PLOT_STEP4
<matplotlib.figure.Figure at 0x7fa6ec704fd0>



We run the time evolution many times now varying the detuning of the laser.

>>> tr=run_evolution(path,name+"_evolution",[E00],[-20.0],  2000,  0.01,Ne,spectrum_of_laser=1,N_delta=401,frequency_end=20.0,use_netcdf=use_netcdf)
    
>>> delta,rho22,rho21_real,rho21_imag = read_result(path,name+"_evolution",N=Ne,use_netcdf=use_netcdf)

>>> fig=pyplot.figure(); ax=fig.add_subplot(111)
>>> ax.plot(delta,rho22,'k-'         ,label=r"$\rho_{22}$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(delta,rho21_real,'b-'    ,label=r"$\mathfrak{Re}\rho_{21}$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(delta,rho21_imag,'r-'    ,label=r"$\mathfrak{Im}\rho_{21}$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.legend(loc=0,fontsize=15) # doctest: +IGNORE_PLOT_STEP2
>>> ax.set_xlabel(r"$\delta$",fontsize=20) # doctest: +IGNORE_PLOT_STEP2
    
>>> ax.set_xlim([-20,20]) # doctest: +IGNORE_PLOT_STEP3
>>> pyplot.savefig(path+'spectrum_'+name+'_evolution.png',bbox_inches='tight') # doctest: +IGNORE_PLOT_STEP4
<matplotlib.figure.Figure at 0x7fa6bf4fabd0>



## Steady state
We write the Fortran code of the experiment.

>>> tw=write_stationary(path,name+"_steady",laseres,omega,gamma,r,Lij)

We compile it.

>>> tc=compile_code(path,name+"_steady",lapack=True, parallel=parallel)

We explain the initial detuning and the electric field amplitude.

>>> delta=-20; E00=1.0

We run the spectroscopy.

>>> tr=run_stationary(path,name+"_steady",[E00],[delta],spectrum_of_laser=1,N_delta=401,frequency_end=20.0,use_netcdf=use_netcdf)
    
>>> delta,rho22,rho21_real,rho21_imag = read_result(path,name+"_steady",N=Ne,use_netcdf=use_netcdf)

>>> fig=pyplot.figure(); ax = fig.add_subplot(111)
>>> ax.plot(delta,rho22,'k-'         ,label=r"$\rho_{22}$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(delta,rho21_real,'b-'    ,label=r"$\mathfrak{Re}\rho_{21}$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(delta,rho21_imag,'r-'    ,label=r"$\mathfrak{Im}\rho_{21}$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.legend(loc=0,fontsize=15) # doctest: +IGNORE_PLOT_STEP2
>>> ax.set_xlabel(r"$\delta$",fontsize=20) # doctest: +IGNORE_PLOT_STEP2
    
>>> ax.set_xlim([-20,20]) # doctest: +IGNORE_PLOT_STEP3
>>> pyplot.savefig(path+'spectrum_'+name+'.png',bbox_inches='tight') # doctest: +IGNORE_PLOT_STEP4
<matplotlib.figure.Figure at 0x7fa6c4e2e190>



### Power broadening.

>>> if True:
...     delta1=-20.0
    
...     fig=pyplot.figure(); ax = fig.add_subplot(111)
...     E0=[0.25, 0.5, 1.0, 3.0, 10.0]
...     colors=['m','b','g',"orange",'r']
...     for i in range(len(E0)):
...         tr=run_stationary(path,name+"_steady",[E0[i]],[delta1],spectrum_of_laser=1,N_delta=401,frequency_end=20.0,use_netcdf=use_netcdf)
...         delta,rho22,rho21_real,rho21_imag = read_result(path,name+"_steady",N=Ne,use_netcdf=use_netcdf)
...         
...         ax.plot(delta,rho22,'-',color=colors[i],label=r"$E_0="+str(E0[i])+"$") # doctest: +IGNORE_PLOT_STEP1
...         
...         Omega=E0[i]
...         gamma=1.0
...         
...         hwfm=sqrt(gamma**2+ 2*(Omega)**2)/2
...         ax.plot([-hwfm,hwfm],[rho22[200]/2,rho22[200]/2],color=colors[i]) # doctest: +IGNORE_PLOT_STEP1
...         
...         a=0.025
...         ax.plot([-hwfm,-hwfm],[rho22[200]/2+a*rho22[200],rho22[200]/2-a*rho22[200]],color=colors[i]) # doctest: +IGNORE_PLOT_STEP1
...         ax.plot([ hwfm, hwfm],[rho22[200]/2+a*rho22[200],rho22[200]/2-a*rho22[200]],color=colors[i]) # doctest: +IGNORE_PLOT_STEP1
    
...     ax.legend(loc=0,fontsize=12) # doctest: +IGNORE_PLOT_STEP2
...     ax.set_xlim([-20,20]) # doctest: +IGNORE_PLOT_STEP3
...     ax.set_xlabel(r"$\delta$",fontsize=22) # doctest: +IGNORE_PLOT_STEP2
...     ax.set_ylabel(r"$\rho_{22}$",fontsize=22) # doctest: +IGNORE_PLOT_STEP2
    
...     pyplot.savefig(path+name+'_power_broadening.png',bbox_inches='tight') # doctest: +IGNORE_PLOT_STEP4
... 
<matplotlib.figure.Figure at 0x7fa6bfd29d10>



>>> pyplot.close("all")

[]

"""
__doc__=__doc__.replace("+IGNORE_PLOT_STEP1", "+ELLIPSIS\n[<...>]")
__doc__=__doc__.replace("+IGNORE_PLOT_STEP2", "+ELLIPSIS\n<...>")
__doc__=__doc__.replace("+IGNORE_PLOT_STEP3", "+ELLIPSIS\n(...)")
__doc__=__doc__.replace("+IGNORE_PLOT_STEP4", "\n")
