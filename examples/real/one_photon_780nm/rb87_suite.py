# -*- coding: utf-8 -*-

#************************************************************************
#       Copyright (C) 2014 - 2017 Oscar Gerardo Lazo Arjona             *
#              <oscar.lazo@correo.nucleares.unam.mx>                    *
#                                                                       *
#  This file is part of FAST.                                           *
#                                                                       *
#  FAST is free software: you can redistribute it and/or modify         *
#  it under the terms of the GNU General Public License as published by *
#  the Free Software Foundation, either version 3 of the License, or    *
#  (at your option) any later version.                                  *
#                                                                       *
#  FAST is distributed in the hope that it will be useful,              *
#  but WITHOUT ANY WARRANTY; without even the implied warranty of       *
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
#  GNU General Public License for more details.                         *
#                                                                       *
#  You should have received a copy of the GNU General Public License    *
#  along with FAST.  If not, see <http://www.gnu.org/licenses/>.        *
#                                                                       *
#************************************************************************

from fast import *
from matplotlib import pyplot,cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.optimize import curve_fit
from fast.config import parallel, use_netcdf, run_long_tests

pyplot.ioff()

########################################################################
# In this test we calculate the stationary state of Rb 87 using only one
# laser field to drive the levels 5 S_1/2 and 5 P_3/2.
########################################################################

################################################
#We define the basic characteristics of the experiment.
#The path in which files will be placed.
path='./real/one_photon_780nm/'
#The name of te experiments.
name='rb87_suite'

#First we specify the states we will use
element="Rb"
isotope=87
e1=State(element,isotope,5,S,Integer(1)/2)
e3=State(element,isotope,5,P,Integer(3)/2)
fine_states=[e1,e3]

#Next we split these fine states into their hyperfine substates.
hyperfine_states=split_fine_to_hyperfine(fine_states)

#Next we split these hyperfine states into their magnetic substates.
magnetic_states=split_hyperfine_to_magnetic(hyperfine_states)

#The number of states
Ne=len(magnetic_states)
#The number of lasers
Nl=1
#We choose a frequency scale for the spectra
Omega=1e6

#We calculate the matrices for the given states
omega,gamma,r=calculate_matrices(magnetic_states,Omega)

#We plot these matrices.
fig=pyplot.figure(); ax=fig.add_subplot(111)
fancy_matrix_plot(ax,omega,magnetic_states,path,name+'_omega.png',take_abs=True,colorbar=True); fig=pyplot.figure(); ax=fig.add_subplot(111)
fancy_matrix_plot(ax,gamma,magnetic_states,path,name+'_gamma.png',take_abs=True,colorbar=True); fig=pyplot.figure(); ax=fig.add_subplot(111)
fancy_r_plot(        r    ,magnetic_states,path,name+'_r.png'    , complex_matrix=True)

#~ #We check that the sum of all decays is 6.065 MHz.
#~ #print [sum([gamma[i][j] for j in range(i)])/2/pi for i in range(Ne)]
#~ 
#~ #We give the sets of excitations per laser.
Lij=[]
for i in range(Ne):
    for j in range(i):
        if r[1][i][j]!=0 or r[0][i][j]!=0 or r[2][i][j]!=0:
        #if r[1][i][j]!=0:
            Lij+=[[j+1,i+1,[1]]]
            
Lij=formatLij(Lij,Ne)
#And we plot them. In this case this shows only blue squares,
#because all the allowed excitations are driven by the same laser.
#If there were more than one laser it would show different
#excitations driven by different lasers as differently colored squares.
#fig=pyplot.figure(); ax=fig.add_subplot(111)
#plot_Lij(ax,Lij,Nl,magnetic_states,path,name+'_Lij.png')

########################################################################
#We define the lasers.
l1=PlaneWave(0,0,0, pi/8,symbolical=False)#Circular propagating through z
lasers=[l1]
Nl=len(lasers)
fig = pyplot.figure(); ax = fig.gca(projection='3d')
draw_lasers_3d(ax,lasers,path+name+'_lasers.png')

########################################################################
#We draw a level diagram.
fig=pyplot.figure(); ax=fig.add_subplot(111,aspect="equal")
a=25; b=30
#ax.plot([0,a,a,0],[0,0,b,b],'w-',alpha=0.0)

xs=2.5; xp=22.5
p1 =[xs,  1.0]; w1 =5.0; h1 = 5.0
p2 =[xp, 15.0]; w2 =5.0; h2 = 2.0

a=False;a=True
h1=draw_multiplet(ax,e1 ,p1 ,h1 ,w1 ,fside='right',label_fontsize=23,label_separation=0.5,
												  fsize=15,deltanu_fontsize=12,text=e1._latex_()[16:],text_pos='right',magnetic_lines=a)
h2=draw_multiplet(ax,e3 ,p2 ,h2 ,w2 ,fside='left', label_fontsize=23,label_separation=0.75,
												  fsize=15,deltanu_fontsize=12,text=e3._latex_()[16:],text_pos='left',magnetic_lines=a)
excitation(ax,[h1[1][0]-0.25*5,h1[1][1]],[h2[3][0]-0.25*5,h2[3][1]], fc="red", ec="red",width=0.2, head_width=2, head_length=2)
decay(     ax,[h1[0][0]+0.25*5,h1[0][1]],[h2[2][0]-0.25*5,h2[2][1]], 0.5,10.0,color="red",linewidth=2.0)
ax.text(5,25,r"$^{"+str(isotope)+"}\mathrm{Rb}$",fontsize=45,verticalalignment="center",horizontalalignment="center")

pyplot.axis('off')
pyplot.savefig(path+name+'_levels.png',bbox_inches="tight")
pyplot.savefig(path+name+'_levels.pdf',bbox_inches="tight")

########################################################################
#We now calculate the electric field amplitude for a laser beam of
#power 7.25e-3 Watts and a gaussian profile with a 2.5e-3 m standard deviation.
E00=electric_field_amplitude_gaussian(0.005e-3,2.5e-3,Omega=Omega)
#print E00,'E00'
s0=1.0e1
E00=electric_field_amplitude_intensity(s0,Omega=Omega)
#print E00,'E00'
#We set the electric field amplitude.
E0=[E00]

########################################################################
#We calculate the transition frequencies (in MHz).
#The smallest transition frequency.
nu0=Transition(hyperfine_states[3],hyperfine_states[1]).nu
#The first set of frequencies (from 5 S_1/2 F=2 to 5 P_3/2 F=1,2,3).
t1=[(Transition(hyperfine_states[i+3],hyperfine_states[1]).nu -nu0)/Omega for i in range(3)]
#The first set of frequencies (from 5 S_1/2 F=1 to 5 P_3/2 F=0,1,2).
t2=[(Transition(hyperfine_states[i+2],hyperfine_states[0]).nu -nu0)/Omega for i in range(3)]

########################################################################
#                           Time Evolution                             #
########################################################################
#We write the Fortran code of the experiment.
write_evolution(path,name+"_evolution",lasers,omega,gamma,r,Lij,verbose=0)
#We compile it.
compile_code(path,name+"_evolution",lapack=True,optimization_flag='',parallel=parallel)

########################################################################
#We see what happens during 1 us of evolution.
dt=0.001; N_iter=1000
#With the laser tuned to the cyclic upper transition.
detuning_knob=[t1[2]*2*pi]

#The initial state with have the population evenly distributed among the
#magnetic states of 5S1/2.
rho0=[1/8.0 for i in range(7)]+[0.0 for j in range(24-8)]

run_evolution(path,name+"_evolution",E0,detuning_knob, N_iter,dt,Ne,	rho0=rho0,use_netcdf=use_netcdf)
###########################################################
dat=read_result(path,name+"_evolution",N=Ne,use_netcdf=use_netcdf)
t=dat[0]
popfg1=dat[3-2:3]; pop1= [1-sum([dat[j+1][i] for j in range(23)]) for i in range(len(t))]
popfg1=[pop1]+popfg1
popfg2=dat[8-5:8]
popfmax=dat[24-7:24]

###########################################################
fig=pyplot.figure(); ax=fig.add_subplot(111)

ax.plot(t,popfg1[ 2],'-' ,color='r'		,label=r"$\mathrm{Population} \ M_F=1$")
ax.plot(t,popfg1[ 1],'-' ,color='b'		,label=r"$\mathrm{Population} \ M_F=0$")
ax.plot(t,popfg1[ 0],'--',color='r'		,label=r"$\mathrm{Population} \ M_F=-1$")

ax.legend(loc=0)
ax.set_xlabel(r"$t\ (\mu\mathrm{s})$",fontsize=20)
ax.set_ylim([0,0.14])
pyplot.savefig(path+name+"_tshort_fg1.png",bbox_inches="tight")
pyplot.savefig(path+name+"_tshort_fg1.pdf",bbox_inches="tight")

###########################################################
fig=pyplot.figure(); ax=fig.add_subplot(111)
ax.plot(t,popfg2[ 4],'-' ,color='r'		,label=r"$\mathrm{Population} \ M_F=2$")
ax.plot(t,popfg2[ 3],'-' ,color='g'		,label=r"$\mathrm{Population} \ M_F=1$")
ax.plot(t,popfg2[ 2],'-' ,color='b'		,label=r"$\mathrm{Population} \ M_F=0$")
ax.plot(t,popfg2[ 1],'--' ,color='g'		,label=r"$\mathrm{Population} \ M_F=-1$")
ax.plot(t,popfg2[ 0],'--' ,color='r'		,label=r"$\mathrm{Population} \ M_F=-2$")

ax.legend(loc=0)
ax.set_xlabel(r"$t\ (\mu\mathrm{s})$",fontsize=20)
ax.set_ylim([0,None])
pyplot.savefig(path+name+"_tshort_fg2.png",bbox_inches="tight")
pyplot.savefig(path+name+"_tshort_fg2.pdf",bbox_inches="tight")

###########################################################
fig=pyplot.figure(); ax=fig.add_subplot(111)
ax.plot(t,popfmax[ 6],'-' ,color='r'		,label=r"$\mathrm{Population} \ M_F=3$")
ax.plot(t,popfmax[ 5],'-' ,color='orange'	,label=r"$\mathrm{Population} \ M_F=2$")
ax.plot(t,popfmax[ 4],'-' ,color='g'		,label=r"$\mathrm{Population} \ M_F=1$")
ax.plot(t,popfmax[ 3],'-' ,color='b'		,label=r"$\mathrm{Population} \ M_F=0$")
ax.plot(t,popfmax[ 2],'--',color='g'		,label=r"$\mathrm{Population} \ M_F=-1$")
ax.plot(t,popfmax[ 1],'--',color='orange'	,label=r"$\mathrm{Population} \ M_F=-2$")
ax.plot(t,popfmax[ 0],'--',color='r'		,label=r"$\mathrm{Population} \ M_F=-3$")

ax.legend(loc=0)
ax.set_xlabel(r"$t\ (\mu\mathrm{s})$",fontsize=20)
ax.set_ylim([0,None])
pyplot.savefig(path+name+"_tshort_fmax.png",bbox_inches="tight")
pyplot.savefig(path+name+"_tshort_fmax.pdf",bbox_inches="tight")

########################################################################
#We see what happens during 1 s of evolution.
dt=1.0e3; N_iter=1000
#With the laser tuned to the cyclic upper transition.
detuning_knob=[t1[2]*2*pi]

#The initial state with have the population evenly distributed among the
#magnetic states of 5S1/2.
rho0=[1/8.0 for i in range(7)]+[0.0 for j in range(24-8)]

run_evolution(path,name+"_evolution",E0,detuning_knob, N_iter,dt,Ne,	rho0=rho0,use_netcdf=use_netcdf)
###########################################################
dat=read_result(path,name+"_evolution",N=Ne,use_netcdf=use_netcdf)
t=dat[0]
t=[ti/1e6 for ti in t]
popfg1=dat[3-2:3]; pop1= [1-sum([dat[j+1][i] for j in range(23)]) for i in range(len(t))]
popfg1=[pop1]+popfg1
popfg2=dat[8-5:8]
popfmax=dat[24-7:24]

###########################################################
fig=pyplot.figure(); ax=fig.add_subplot(111)

ax.plot(t,popfg1[ 2],'-' ,color='r'		,label=r"$\mathrm{Population} \ M_F=1$")
ax.plot(t,popfg1[ 1],'-' ,color='b'		,label=r"$\mathrm{Population} \ M_F=0$")
ax.plot(t,popfg1[ 0],'--',color='r'		,label=r"$\mathrm{Population} \ M_F=-1$")

ax.legend(loc=0)
ax.set_xlabel(r"$t\ (\mathrm{s})$",fontsize=20)
ax.set_ylim([0,None])
pyplot.savefig(path+name+"_tlong_fg1.png",bbox_inches="tight")
pyplot.savefig(path+name+"_tlong_fg1.pdf",bbox_inches="tight")

###########################################################
fig=pyplot.figure(); ax=fig.add_subplot(111)
ax.plot(t,popfg2[ 4],'-' ,color='r'		,label=r"$\mathrm{Population} \ M_F=2$")
ax.plot(t,popfg2[ 3],'-' ,color='g'		,label=r"$\mathrm{Population} \ M_F=1$")
ax.plot(t,popfg2[ 2],'-' ,color='b'		,label=r"$\mathrm{Population} \ M_F=0$")
ax.plot(t,popfg2[ 1],'--' ,color='g'		,label=r"$\mathrm{Population} \ M_F=-1$")
ax.plot(t,popfg2[ 0],'--' ,color='r'		,label=r"$\mathrm{Population} \ M_F=-2$")

ax.legend(loc=0)
ax.set_xlabel(r"$t\ (\mathrm{s})$",fontsize=20)
ax.set_ylim([0,None])
pyplot.savefig(path+name+"_tlong_fg2.png",bbox_inches="tight")
pyplot.savefig(path+name+"_tlong_fg2.pdf",bbox_inches="tight")

###########################################################
fig=pyplot.figure(); ax=fig.add_subplot(111)
ax.plot(t,popfmax[ 6],'-' ,color='r'		,label=r"$\mathrm{Population} \ M_F=3$")
ax.plot(t,popfmax[ 5],'-' ,color='orange'	,label=r"$\mathrm{Population} \ M_F=2$")
ax.plot(t,popfmax[ 4],'-' ,color='g'		,label=r"$\mathrm{Population} \ M_F=1$")
ax.plot(t,popfmax[ 3],'-' ,color='b'		,label=r"$\mathrm{Population} \ M_F=0$")
ax.plot(t,popfmax[ 2],'--',color='g'		,label=r"$\mathrm{Population} \ M_F=-1$")
ax.plot(t,popfmax[ 1],'--',color='orange'	,label=r"$\mathrm{Population} \ M_F=-2$")
ax.plot(t,popfmax[ 0],'--',color='r'		,label=r"$\mathrm{Population} \ M_F=-3$")

ax.legend(loc=0)
ax.set_xlabel(r"$t\ (\mathrm{s})$",fontsize=20)
ax.set_ylim([0,None])
pyplot.savefig(path+name+"_tlong_fmax.png",bbox_inches="tight")
pyplot.savefig(path+name+"_tlong_fmax.pdf",bbox_inches="tight")
pyplot.close("all")

########################################################################
#                           Steady State                               #
########################################################################
#We write the Fortran code of the experiment.
write_stationary(path,name+"_steady",lasers,omega,gamma,r,Lij,verbose=0)
#We compile it.
compile_code(path,name+"_steady",lapack=True,parallel=parallel)

####################################################################
#We now make a test of this stationary two-level state to make sure.
run_long_tests=1
if run_long_tests>=1:
	def model_lorentzian(n=1):
		s='lambda x,'
		for i in range(1,n+1): s+=' A'+str(i)+','
		for i in range(1,n+1): s+=' x0'+str(i)+','
		for i in range(1,n+1): s+=' gamma'+str(i)+','
		s=s[:-1]
		s+=' :'
		for i in range(1,n+1): s+=' A'+str(i)+'*gamma'+str(i)+'**2/(4*(x-x0'+str(i)+')**2+gamma'+str(i)+'**2) +'
		s=s[:-1]
		return eval(s)

	def fit_lorentzians(x,y,p0=None,n=1,return_p0=False):
		'''Fits a lorentzian curve using p0=[A1,A2,...,x1,x2,...,gamma1,gamma2,...] as an initial guess.
		It returns a curve with N_points.'''
		
		lorentzians=model_lorentzian(n)
		N=len(x)
		if return_p0:
			fit=p0;pcov=None
		else:
			fit,pcov = curve_fit(lorentzians,x,y,p0=p0)
		s=''.join([',fit['+str(i)+']' for i in range(3*n)])
		s='[lorentzians(xi'+s+') for xi in x]'
		fitted_curve=eval(s)
		
		for i in range(n):
			fit[-i-1]=abs(fit[-i-1])
		return fit,pcov,fitted_curve

	########################################################################
	Npow=19; pow_ini=-2.0; pow_end=2.0; pow_step=(pow_end-pow_ini)/(Npow-1)
	pows=[10**(pow_ini+i*pow_step) for i in range(Npow)]
	model_widths=[]; simple_widths=[]
	model_pops=[]; simple_pops=[]
	detuning_knob=[(t1[2]-30)*2*pi]
	frequency_end=(t1[2]+30)*2*pi
	#print 111
	N_delta=250
	for s0 in pows:
		
		E000=[electric_field_amplitude_intensity(s0,Omega=Omega)]
		
		run_stationary(path,name+"_steady",E000,detuning_knob,spectrum_of_laser=1,N_delta=N_delta, frequency_end=frequency_end,use_netcdf=use_netcdf)

		nu,pop0=read_result(path,name+"_steady",Ne,Ne,1,Ne,use_netcdf=use_netcdf)
		nu=[nui/2/pi for nui in nu]
		
		fit,pcov,fitted_curve=fit_lorentzians(nu,pop0,[pop0[N_delta/2],t1[2],1.04],n=1,return_p0=False)
		pyplot.plot(nu,pop0,'r+')
		pyplot.plot(nu,fitted_curve,'b-')
		Ai,x0i,gammai=fit
		
		model_widths+=[gammai]
		Gamma=6.065
		simple_widths+=[Gamma*sqrt(1+2*s0)]
		
		model_pops+=[Ai]
		simple_pops+=[2*s0/2/(1+2*s0)]
		
		print 'I/I0=',s0,'E0=',E000[0]

	pyplot.xlabel(r"$\delta$",fontsize=18)
	pyplot.ylabel(r"$\rho_{22}$",fontsize=18)
	pyplot.savefig(path+name+"_1power.png",bbox_inches="tight")
	pyplot.close("all")

	Npow=1000; pow_ini=-2.0; pow_end=2.0; pow_step=(pow_end-pow_ini)/(Npow-1)
	pows2=[10**(pow_ini+i*pow_step) for i in range(Npow)]

	#We plot widths
	simple_widths =[6.065*sqrt(1+s0) for s0 in pows2]

	pyplot.semilogx(pows2,simple_widths,'-b',basex=10,label=r'$\Gamma\sqrt{1+I/I_0}$')
	pyplot.semilogx(pows,model_widths,'+r',basex=10,label=r'$\mathrm{Modeled \ width}$',markersize=13)
	pyplot.plot(pows2,[6.065]*Npow,'g-',label=r'$\mathrm{Natural \ width}$')

	pyplot.xlabel(r"$I/I_0$",fontsize=18)
	pyplot.ylabel(r"$\mathrm{Width \ (MHz)}$",fontsize=18)

	pyplot.ylim([0,None])
	pyplot.legend(loc=0)

	pyplot.savefig(path+name+'_2power.png',bbox_inches='tight')
	pyplot.close()


	#We plot amplitudes
	simple_pops =[ s0/2/(1+ s0) for s0 in pows2]

	pyplot.semilogx(pows2,simple_pops,'-b',basex=10,label=r'$\frac{I/I_0}{2(1+I/I_0)}$')
	pyplot.semilogx(pows,model_pops,'+r',basex=10,label=r'$\mathrm{Modeled \ \rho_{ee}}$',markersize=13)
	pyplot.xlabel(r"$I/I_0$",fontsize=18)
	pyplot.ylabel(r"$\mathrm{Peak\ height}$",fontsize=18)


	pyplot.ylim([0,None])
	pyplot.legend(loc=0)

	pyplot.ylim([0,None])
	pyplot.legend(loc=0)

	pyplot.savefig(path+name+'_3power.png',bbox_inches='tight')
	pyplot.close()
