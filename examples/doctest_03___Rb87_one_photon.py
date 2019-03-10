# -*- coding: utf-8 -*-
# Copyright (C) 2017 Oscar Gerardo Lazo Arjona
# mailto: oscar.lazo@correo.nucleares.unam.mx

__doc__ = r"""

>>> from fast import State, Transition, Integer, make_list_of_states, calculate_matrices
>>> from fast import fancy_matrix_plot, fancy_r_plot, PlaneWave
>>> from fast import draw_lasers_3d, formatLij, draw_multiplet, excitation, decay, electric_field_amplitude_intensity
>>> from fast import write_evolution, compile_code, run_evolution, read_result, write_stationary, run_stationary

>>> from math import pi,sqrt
>>> from matplotlib import pyplot
>>> from fast.config import parallel, use_netcdf
>>> from scipy.optimize import curve_fit
>>> from numpy import array
>>> from time import time

>>> path = "folder_03___Rb87_one_photon/" 
>>> name = 'suite' 

We first specify the states that we will use.

>>> element = "Rb" 
>>> isotope = 87
>>> e1 = State(element, isotope, 5, 0, Integer(1)/2)
>>> e3 = State(element, isotope, 5, 1, Integer(3)/2)
>>> fine_states = [e1, e3]
>>> print(fine_states)
[87Rb 5S_1/2, 87Rb 5P_3/2]



We split these states into hyperfine substates.

>>> hyperfine_states = make_list_of_states(fine_states, "hyperfine")
>>> print(hyperfine_states)
[87Rb 5S_1/2^1, 87Rb 5S_1/2^2, 87Rb 5P_3/2^0, 87Rb 5P_3/2^1, 87Rb 5P_3/2^2, 87Rb 5P_3/2^3]



We divide these hypefine states into their magnetic substates.

>>> magnetic_states = make_list_of_states(fine_states, "magnetic")

We define the total number of states.

>>> Ne = len(magnetic_states)

We choose a frequency scale (and therefore a time scale).

>>> Omega = 1e6

We calculate the matrices for all the given states.

>>> omega, gamma, r = calculate_matrices(magnetic_states, Omega)

We plot the resonant frequencies $\\omega_{ij}$.

>>> fig = pyplot.figure(); ax=fig.add_subplot(111)
>>> fancy_matrix_plot(ax, omega, magnetic_states, path, name+'_omega.png', take_abs=True, colorbar=True) # doctest: +IGNORE_PLOT_STEP4
<Figure size 432x288 with 1 Axes>



We plot the matrix of decays $\\gamma_{ij}$.

>>> fig=pyplot.figure(); ax=fig.add_subplot(111)
>>> fancy_matrix_plot(ax,gamma,magnetic_states,path,name+'_gamma.png',take_abs=True,colorbar=True) # doctest: +IGNORE_PLOT_STEP4
<Figure size 432x288 with 1 Axes>



We plot the components of the position operator $\\hat{\\vec{r}}$.

>>> fancy_r_plot(r         ,magnetic_states,path,name+'_r.png', complex_matrix=True) # doctest: +IGNORE_PLOT_STEP4
<Figure size 1080x342 with 3 Axes>



We define the laser we will use.

>>> l1 = PlaneWave(0, 0, 0, pi/8)
>>> lasers = [l1]
>>> Nl = len(lasers)
>>> fig = pyplot.figure(); ax = fig.gca(projection='3d')
>>> draw_lasers_3d(ax,lasers,path+name+'_lasers.png') # doctest: +IGNORE_PLOT_STEP4
<Figure size 432x288 with 1 Axes>



We specify the couplings of the lasers. In this case, all allowed transitions are coupled by laser 1.

>>> Lij=[]
>>> for i in range(Ne):
...     for j in range(i):
...         if r[1][i][j]!=0 or r[0][i][j]!=0 or r[2][i][j]!=0:
...             Lij+=[[j+1,i+1,[1]]]
...             
>>> Lij=formatLij(Lij,Ne)
>>> #fig=pyplot.figure(); ax=fig.add_subplot(111)
>>> #plot_Lij(ax,Lij,Nl,magnetic_states,path,name+'_Lij.png')

We draw a level diagram.

>>> sc=2.0
>>> fig=pyplot.figure(figsize=(sc*2.5,sc*3.0)); ax=fig.add_subplot(111,aspect="equal")
>>> a=25; b=30
    
>>> xs=2.5; xp=22.5
>>> p1 =[xs,  1.0]; w1 =5.0; h1 = 5.0
>>> p2 =[xp, 15.0]; w2 =5.0; h2 = 2.0
    
>>> a=False;a=True
>>> h1=draw_multiplet(ax,e1 ,p1 ,h1 ,w1 ,fside='right',label_fontsize=23,label_separation=0.5,
...                     fsize=15,deltanu_fontsize=12,text=e1._latex_()[16:],text_pos='right',magnetic_lines=a)
>>> h2=draw_multiplet(ax,e3 ,p2 ,h2 ,w2 ,fside='left', label_fontsize=23,label_separation=0.75,
...                     fsize=15,deltanu_fontsize=12,text=e3._latex_()[16:],text_pos='left',magnetic_lines=a)
    
>>> excitation(ax,[h1[1][0]-0.25*5,h1[1][1]],[h2[3][0]-0.25*5,h2[3][1]], fc="red", ec="red",width=0.2, head_width=2, head_length=2) # doctest: +IGNORE_PLOT_STEP4
>>> decay(     ax,[h1[0][0]+0.25*5,h1[0][1]],[h2[2][0]-0.25*5,h2[2][1]], 0.5,10.0,color="red",linewidth=2.0) # doctest: +IGNORE_PLOT_STEP4
>>> ax.text(5,25,r"$^{"+str(isotope)+"}\mathrm{Rb}$",fontsize=45,verticalalignment="center",horizontalalignment="center") # doctest: +IGNORE_PLOT_STEP5
    
>>> pyplot.axis('off') # doctest: +IGNORE_PLOT_STEP3
>>> pyplot.savefig(path+name+'_levels.png',bbox_inches="tight") # doctest: +IGNORE_PLOT_STEP4
>>> pyplot.savefig(path+name+'_levels.pdf',bbox_inches="tight") # doctest: +IGNORE_PLOT_STEP4
<Figure size 360x432 with 1 Axes>



We calculate the amplitude of the electric field corresponding to saturation.

>>> s0 = 1.0e1
>>> E00 = electric_field_amplitude_intensity(s0, Omega=Omega)
>>> E0 = [E00]
>>> print(E0)
[28.50887690495697]



We calculate the transition frequencies of the (in MHz). The transitions $5 S_1/2 F=2 \\rightarrow 5 P_3/2 F=1,2,3$ (the big F's).

>>> nu0=Transition(hyperfine_states[3],hyperfine_states[1]).nu
>>> t1=[(Transition(hyperfine_states[i+3],hyperfine_states[1]).nu -nu0)/Omega for i in range(3)]
>>> print(t1)
[0.0, 156.940152625, 423.5917535625]



The transitions $5 S_1/2 F=1 \\rightarrow 5 P_3/2 F=0,1,2$ (the small F's).

>>> t2 = [(Transition(hyperfine_states[i+2],hyperfine_states[0]).nu -nu0)/Omega for i in range(3)]
>>> print(t2)
[6762.460809, 6834.682610875, 6991.6227635]



# Time evolution
We write the Fortran code of the experiment.

>>> tw = write_evolution(path, name+"_evolution", lasers, omega, gamma, r, Lij, verbose=0)

We compile it.

>>> tc = compile_code(path,name+"_evolution",lapack=True,optimization_flag='',parallel=parallel)

We begin with an initial state where all population is distributed equally between the $5S_{1/2}$ states.

>>> rho0 = [1/8.0 for i in range(7)]+[0.0 for j in range(24-8)]

We se what happens during 1 $\\mu$s of evolution.

>>> dt=1e-3; N_iter=1000
>>> detuning_knob=[t1[2]*2*pi]
    
>>> tr=run_evolution(path,name+"_evolution",E0,detuning_knob, N_iter,dt,Ne,rho0=rho0,use_netcdf=use_netcdf)
    
>>> dat=read_result(path,name+"_evolution",N=Ne,use_netcdf=use_netcdf)
>>> t=dat[0]
>>> popfg1=dat[3-2:3]; pop1= [1-sum([dat[j+1][i] for j in range(23)]) for i in range(len(t))]
>>> popfg1=[pop1]+popfg1
>>> popfg2=dat[8-5:8]
>>> popfmax=dat[24-7:24]

We plot what happens to the states $5S_{1/2}F=1$.

>>> fig=pyplot.figure(); ax=fig.add_subplot(111)
    
>>> ax.plot(t,popfg1[ 2],'-' ,color='r',label=r"$M_F=1$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(t,popfg1[ 1],'-' ,color='b',label=r"$M_F=0$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(t,popfg1[ 0],'--',color='r',label=r"$M_F=-1$") # doctest: +IGNORE_PLOT_STEP1
    
>>> ax.legend(loc=0) # doctest: +IGNORE_PLOT_STEP2
>>> ax.set_xlabel(r"$t\ (\mu\mathrm{s})$",fontsize=20) # doctest: +IGNORE_PLOT_STEP5
>>> ax.set_ylabel(r"$\mathrm{Population}$",fontsize=20) # doctest: +IGNORE_PLOT_STEP5
>>> ax.set_ylim([0,0.14]) # doctest: +IGNORE_PLOT_STEP3
    
>>> pyplot.savefig(path+name+"_tshort_fg1.png",bbox_inches="tight") # doctest: +IGNORE_PLOT_STEP4
>>> pyplot.savefig(path+name+"_tshort_fg1.pdf",bbox_inches="tight") # doctest: +IGNORE_PLOT_STEP4
<Figure size 432x288 with 1 Axes>



We plot what happens to states $5S_{1/2}F=2$.

>>> fig=pyplot.figure(); ax=fig.add_subplot(111)
>>> ax.plot(t,popfg2[ 4],'-'  ,color='r',label=r"$\mathrm{Population} \ M_F=2$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(t,popfg2[ 3],'-'  ,color='g',label=r"$\mathrm{Population} \ M_F=1$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(t,popfg2[ 2],'-'  ,color='b',label=r"$\mathrm{Population} \ M_F=0$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(t,popfg2[ 1],'--' ,color='g',label=r"$\mathrm{Population} \ M_F=-1$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(t,popfg2[ 0],'--' ,color='r',label=r"$\mathrm{Population} \ M_F=-2$") # doctest: +IGNORE_PLOT_STEP1
    
>>> ax.legend(loc=0) # doctest: +IGNORE_PLOT_STEP2
>>> ax.set_xlabel(r"$t\ (\mu\mathrm{s})$",fontsize=20) # doctest: +IGNORE_PLOT_STEP5
>>> ax.set_ylim([0,None]) # doctest: +IGNORE_PLOT_STEP3
>>> pyplot.savefig(path+name+"_tshort_fg2.png",bbox_inches="tight") # doctest: +IGNORE_PLOT_STEP4
>>> pyplot.savefig(path+name+"_tshort_fg2.pdf",bbox_inches="tight") # doctest: +IGNORE_PLOT_STEP4
<Figure size 432x288 with 1 Axes>



We plot what happens to the states $5P_{3/2}F=3$.

>>> fig=pyplot.figure(); ax=fig.add_subplot(111)
>>> ax.plot(t,popfmax[ 6],'-' ,color='r'     ,label=r"$\mathrm{Population} \ M_F=3$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(t,popfmax[ 5],'-' ,color='orange',label=r"$\mathrm{Population} \ M_F=2$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(t,popfmax[ 4],'-' ,color='g'     ,label=r"$\mathrm{Population} \ M_F=1$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(t,popfmax[ 3],'-' ,color='b'     ,label=r"$\mathrm{Population} \ M_F=0$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(t,popfmax[ 2],'--',color='g'     ,label=r"$\mathrm{Population} \ M_F=-1$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(t,popfmax[ 1],'--',color='orange',label=r"$\mathrm{Population} \ M_F=-2$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(t,popfmax[ 0],'--',color='r'     ,label=r"$\mathrm{Population} \ M_F=-3$") # doctest: +IGNORE_PLOT_STEP1
    
>>> ax.legend(loc=0,fontsize=11) # doctest: +IGNORE_PLOT_STEP2
>>> ax.set_xlabel(r"$t\ (\mu\mathrm{s})$",fontsize=18) # doctest: +IGNORE_PLOT_STEP5
>>> ax.set_ylim([0,None]) # doctest: +IGNORE_PLOT_STEP3
>>> pyplot.savefig(path+name+"_tshort_fmax.png",bbox_inches="tight") # doctest: +IGNORE_PLOT_STEP4
>>> pyplot.savefig(path+name+"_tshort_fmax.pdf",bbox_inches="tight") # doctest: +IGNORE_PLOT_STEP4
<Figure size 432x288 with 1 Axes>



>>> dt=1.0e3; N_iter=1000
>>> detuning_knob=[t1[2]*2*pi]
    
>>> tr=run_evolution(path,name+"_evolution",E0,detuning_knob, N_iter,dt,Ne,rho0=rho0,use_netcdf=use_netcdf)
    
>>> dat=read_result(path,name+"_evolution",N=Ne,use_netcdf=use_netcdf)
>>> t=dat[0]
>>> t=[ti/1e6 for ti in t]
>>> popfg1=dat[3-2:3]; pop1= [1-sum([dat[j+1][i] for j in range(23)]) for i in range(len(t))]
>>> popfg1=[pop1]+popfg1
>>> popfg2=dat[8-5:8]
>>> popfmax=dat[24-7:24]

We plot what happens to states $5S_{1/2}F=1$.

>>> fig=pyplot.figure(); ax=fig.add_subplot(111)
    
>>> ax.plot(t,popfg1[ 2],'-' ,color='r',label=r"$\mathrm{Population} \ M_F=1$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(t,popfg1[ 1],'-' ,color='b',label=r"$\mathrm{Population} \ M_F=0$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(t,popfg1[ 0],'--',color='r',label=r"$\mathrm{Population} \ M_F=-1$") # doctest: +IGNORE_PLOT_STEP1
    
>>> ax.legend(loc=0) # doctest: +IGNORE_PLOT_STEP2
>>> ax.set_xlabel(r"$t\ (\mathrm{s})$",fontsize=20) # doctest: +IGNORE_PLOT_STEP5
>>> ax.set_ylim([0,None]) # doctest: +IGNORE_PLOT_STEP3
>>> pyplot.savefig(path+name+"_tlong_fg1.png",bbox_inches="tight") # doctest: +IGNORE_PLOT_STEP4
>>> pyplot.savefig(path+name+"_tlong_fg1.pdf",bbox_inches="tight") # doctest: +IGNORE_PLOT_STEP4
<Figure size 432x288 with 1 Axes>



We plot what happens to states $5S_{1/2}F=2$.

>>> fig=pyplot.figure(); ax=fig.add_subplot(111)
>>> ax.plot(t,popfg2[ 4],'-'  ,color='r',label=r"$\mathrm{Population} \ M_F=2$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(t,popfg2[ 3],'-'  ,color='g',label=r"$\mathrm{Population} \ M_F=1$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(t,popfg2[ 2],'-'  ,color='b',label=r"$\mathrm{Population} \ M_F=0$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(t,popfg2[ 1],'--' ,color='g',label=r"$\mathrm{Population} \ M_F=-1$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(t,popfg2[ 0],'--' ,color='r',label=r"$\mathrm{Population} \ M_F=-2$") # doctest: +IGNORE_PLOT_STEP1
    
>>> ax.legend(loc=0) # doctest: +IGNORE_PLOT_STEP2
>>> ax.set_xlabel(r"$t\ (\mathrm{s})$",fontsize=20) # doctest: +IGNORE_PLOT_STEP5
>>> ax.set_ylim([0,None]) # doctest: +IGNORE_PLOT_STEP3
>>> pyplot.savefig(path+name+"_tlong_fg2.png",bbox_inches="tight") # doctest: +IGNORE_PLOT_STEP4
>>> pyplot.savefig(path+name+"_tlong_fg2.pdf",bbox_inches="tight") # doctest: +IGNORE_PLOT_STEP4
<Figure size 432x288 with 1 Axes>



We plot what happens to states $5P_{3/2}F=3$.

>>> fig=pyplot.figure(); ax=fig.add_subplot(111)
>>> ax.plot(t,popfmax[ 6],'-' ,color='r'        ,label=r"$\mathrm{Population} \ M_F=3$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(t,popfmax[ 5],'-' ,color='orange'   ,label=r"$\mathrm{Population} \ M_F=2$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(t,popfmax[ 4],'-' ,color='g'        ,label=r"$\mathrm{Population} \ M_F=1$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(t,popfmax[ 3],'-' ,color='b'        ,label=r"$\mathrm{Population} \ M_F=0$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(t,popfmax[ 2],'--',color='g'        ,label=r"$\mathrm{Population} \ M_F=-1$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(t,popfmax[ 1],'--',color='orange'   ,label=r"$\mathrm{Population} \ M_F=-2$") # doctest: +IGNORE_PLOT_STEP1
>>> ax.plot(t,popfmax[ 0],'--',color='r'        ,label=r"$\mathrm{Population} \ M_F=-3$") # doctest: +IGNORE_PLOT_STEP1
    
>>> ax.legend(loc=0) # doctest: +IGNORE_PLOT_STEP2
>>> ax.set_xlabel(r"$t\ (\mathrm{s})$",fontsize=20) # doctest: +IGNORE_PLOT_STEP5
>>> ax.set_ylim([0,None]) # doctest: +IGNORE_PLOT_STEP3
>>> pyplot.savefig(path+name+"_tlong_fmax.png",bbox_inches="tight") # doctest: +IGNORE_PLOT_STEP4
>>> pyplot.savefig(path+name+"_tlong_fmax.pdf",bbox_inches="tight") # doctest: +IGNORE_PLOT_STEP4
<Figure size 432x288 with 1 Axes>



## Steady State
We write the Fortran code of the experiment.

>>> tw = write_stationary(path,name+"_steady",lasers,omega,gamma,r,Lij,verbose=0)

We compile it.

>>> tc=compile_code(path,name+"_steady",lapack=True,parallel=parallel)

We now make a test of this stationary two-level state to make sure.

>>> def model_lorentzian(n=1):
...     s='lambda x,'
...     for i in range(1,n+1): s+=' A'+str(i)+','
...     for i in range(1,n+1): s+=' x0'+str(i)+','
...     for i in range(1,n+1): s+=' gamma'+str(i)+','
...     s=s[:-1]
...     s+=' :'
...     for i in range(1,n+1): s+=' A'+str(i)+'*gamma'+str(i)+'**2/(4*(x-x0'+str(i)+')**2+gamma'+str(i)+'**2) +'
...     s=s[:-1]
...     return eval(s)
    
>>> def fit_lorentzians(x,y,p0=None,n=1,return_p0=False):
...     '''Fits a lorentzian curve using p0=[A1,A2,...,x1,x2,...,gamma1,gamma2,...] as an initial guess.
...     It returns a curve with N_points.'''
...     lorentzians=model_lorentzian(n)
...     N=len(x)
...     if return_p0:
...         fit=p0;pcov=None
...     else:
...         fit,pcov = curve_fit(lorentzians,x,y,p0=p0)
...     s=''.join([',fit['+str(i)+']' for i in range(3*n)])
...     s='[lorentzians(xi'+s+') for xi in x]'
...     fitted_curve=eval(s)    
...     for i in range(n):
...         fit[-i-1]=abs(fit[-i-1])
...     return fit,pcov,fitted_curve
    
>>> ########################################################################
>>> pyplot.close("all")
>>> Npow=19; pow_ini=-2.0; pow_end=2.0; pow_step=(pow_end-pow_ini)/(Npow-1)
>>> pows=[10**(pow_ini+i*pow_step) for i in range(Npow)]
>>> detuning_knob=[(t1[2]-30)*2*pi]
>>> frequency_end=(t1[2]+30)*2*pi
>>> #print 111
>>> Gamma=6.065
>>> N_delta=250
    
>>> def aux2(s0):
...     E000=[electric_field_amplitude_intensity(s0,Omega=Omega)]
...     run_stationary(path,name+"_steady",E000,detuning_knob,spectrum_of_laser=1,N_delta=N_delta, frequency_end=frequency_end,use_netcdf=use_netcdf)
...     nu,pop0=read_result(path,name+"_steady",Ne,Ne,1,Ne,use_netcdf=use_netcdf)
...     nu=[nui/2/pi for nui in nu]
...     res=fit_lorentzians(nu,pop0,[pop0[N_delta/2],t1[2],1.04],n=1,return_p0=False)
...     fit,pcov,fitted_curve=res
...     Ai,x0i,gammai=fit
...     return (gammai,Gamma*sqrt(1+2*s0),Ai,2*s0/2/(1+2*s0))
    
>>> def aux():
...     res=[aux2(s0) for s0 in pows]
...     model_widths =[i[0] for i in res]
...     simple_widths=[i[1] for i in res]
...     model_pops   =[i[2] for i in res]
...     simple_pops  =[i[3] for i in res]
...     return model_widths,simple_widths,model_pops,simple_pops
... 

>>> t0 = time()
>>> model_widths,simple_widths,model_pops,simple_pops=aux()
>>> tlong = time() - t0

>>> pyplot.close("all")
>>> Npow=1000; pow_ini=-2.0; pow_end=2.0; pow_step=(pow_end-pow_ini)/(Npow-1)
>>> pows2=[10**(pow_ini+i*pow_step) for i in range(Npow)]
    
>>> #We plot widths
>>> simple_widths =[6.065*sqrt(1+s0) for s0 in pows2]
    
>>> pyplot.semilogx(pows2,simple_widths,'-b',basex=10,label=r'$\Gamma\sqrt{1+I/I_0}$') # doctest: +IGNORE_PLOT_STEP1
>>> pyplot.semilogx(pows,model_widths,'+r',basex=10,label=r'$\mathrm{Modeled \ width}$',markersize=13) # doctest: +IGNORE_PLOT_STEP1
>>> pyplot.plot(pows2,[6.065]*Npow,'g-',label=r'$\mathrm{Natural \ width}$') # doctest: +IGNORE_PLOT_STEP1
    
>>> pyplot.xlabel(r"$I/I_0$",fontsize=18) # doctest: +IGNORE_PLOT_STEP5
>>> pyplot.ylabel(r"$\mathrm{Width \ (MHz)}$",fontsize=18) # doctest: +IGNORE_PLOT_STEP5
    
>>> pyplot.ylim([0,None]) # doctest: +IGNORE_PLOT_STEP3
>>> pyplot.legend(loc=0) # doctest: +IGNORE_PLOT_STEP2
    
>>> pyplot.savefig(path+name+'_2power.png',bbox_inches='tight') # doctest: +IGNORE_PLOT_STEP4
<Figure size 432x288 with 1 Axes>



>>> pyplot.close()
>>> #We plot amplitudes
>>> simple_pops =[ s0/2/(1+ s0) for s0 in pows2]
>>> pyplot.semilogx(pows2,simple_pops,'-b',basex=10,label=r'$\frac{I/I_0}{2(1+I/I_0)}$') # doctest: +IGNORE_PLOT_STEP1
>>> pyplot.semilogx(pows,model_pops,'+r',basex=10,label=r'$\mathrm{Modeled \ \rho_{ee}}$',markersize=13) # doctest: +IGNORE_PLOT_STEP1
>>> pyplot.semilogx([pows2[0],pows2[-1]],[0.25,0.25],'k-',alpha=0.25) # doctest: +IGNORE_PLOT_STEP1
>>> pyplot.semilogx([1,1],[0,0.5],'k-',alpha=0.25) # doctest: +IGNORE_PLOT_STEP1
    
>>> pyplot.xlabel(r"$I/I_0$",fontsize=18) # doctest: +IGNORE_PLOT_STEP5
>>> pyplot.ylabel(r"$\mathrm{Peak\ height}$",fontsize=18) # doctest: +IGNORE_PLOT_STEP5
    
    
>>> pyplot.ylim([0,None]) # doctest: +IGNORE_PLOT_STEP3
>>> pyplot.legend(loc=0) # doctest: +IGNORE_PLOT_STEP2
    
>>> pyplot.ylim([0,None]) # doctest: +IGNORE_PLOT_STEP3
>>> pyplot.legend(loc=0) # doctest: +IGNORE_PLOT_STEP2
    
>>> pyplot.savefig(path+name+'_3power.png',bbox_inches='tight') # doctest: +IGNORE_PLOT_STEP4
<Figure size 432x288 with 1 Axes>



>>> pyplot.close("all")

[]

"""
__doc__=__doc__.replace("+IGNORE_PLOT_STEP1", "+ELLIPSIS\n[<...>]")
__doc__=__doc__.replace("+IGNORE_PLOT_STEP2", "+ELLIPSIS\n<...>")
__doc__=__doc__.replace("+IGNORE_PLOT_STEP3", "+ELLIPSIS\n(...)")
__doc__=__doc__.replace("+IGNORE_PLOT_STEP5", "+ELLIPSIS\nText(...)")
__doc__=__doc__.replace("+IGNORE_PLOT_STEP4", "\n")
