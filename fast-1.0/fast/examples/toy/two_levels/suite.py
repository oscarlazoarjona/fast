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
from matplotlib import pyplot
from config import parallel, use_netcdf

################################################
#We define the basic characteristics of the experiment.
#The path in which files will be placed.
path='./examples/toy/two_levels/'
#The name of the experiments.
name='suite'
#The number of states.
Ne=2

################################################
#We define the properties of the atom.
#Initialy the frequencies gamma and omega are given in absolute terms (Hz)
#Later the program will rescale to the given scale.
omega_states=[0,100]
omega=[[omega_states[i]-omega_states[j] for j in range(Ne)] for i in range(Ne)]

gamma=[[0.,-1.],[1.,0.]]
r=[[[0,1],[1,0]] for p in range(3)]

################################################
#We define the lasers.
l1=PlaneWave(0,pi/2,0,0,color="red")
lasers=[l1]
Nl=len(lasers)
fig = pyplot.figure(); ax = fig.gca(projection='3d')
draw_lasers_3d(ax,lasers,name=path+name+'_lasers.png')

################################################
#We specify the couplings between states by. In this case,
#that transition |1> -> |2> is coupled by laser 1.
Lij=[[1,2,[1]]]
Lij=formatLij(Lij,Ne)

########################################################################
#We draw a level diagram.
fig=pyplot.figure(); ax=fig.add_subplot(111,aspect="equal")

p1=[0.5,1]; p2=[1.5,3]
draw_state(ax,p1,text=r"$|1\rangle$",l=1.0,alignment='right',label_displacement=0.05,fontsize=25,linewidth=4.0)
draw_state(ax,p2,text=r"$|2\rangle$",l=1.0,alignment='right',label_displacement=0.05,fontsize=25,linewidth=4.0)

excitation(ax,[p1[0]+0.25,p1[1]],[p2[0]+0.25,p2[1]], fc="red", ec="red",width=0.01, head_width=0.1, head_length=0.1)
decay(     ax,[p1[0]-0.25,p1[1]],[p2[0]-0.25,p2[1]], 0.05,10.0,color="red",linewidth=1.0)

pyplot.axis('off')
pyplot.savefig(path+name+'_diagram.png',bbox_inches="tight")

########################################################################
#                           Time Evolution                             #
########################################################################

#We write the Fortran code of the experiment.
write_evolution(path,name+"_evolution",lasers,omega,gamma,r,Lij)
#We compile it.
compile_code(path,name+"_evolution",lapack=True,parallel=parallel)

################################################
#We give the detuning and electric field amplitude initialy in absolute units.
delta1=-1.0; E00=1.0

################################################
#We run the temporal evolution with the parameters we just calculated.
run_evolution(path,name+"_evolution",[E00],[delta1],  2000,  0.01,Ne,use_netcdf=use_netcdf)
#We read the results.
t,rho22,rho21_real,rho21_imag = read_result(path,name+"_evolution",N=Ne,use_netcdf=use_netcdf)

################################################
#We plot the results and save the results for future reference.
fig=pyplot.figure(); ax=fig.add_subplot(111)
ax.plot(t,rho22,'k-'		,label=r"$\rho_{22}$")
ax.plot(t,rho21_real,'b-'	,label=r"$\mathfrak{Re}\rho_{21}$")
ax.plot(t,rho21_imag,'r-'	,label=r"$\mathfrak{Im}\rho_{21}$")

ax.legend(loc=0,fontsize=17)
ax.set_xlabel(r"$t$",fontsize=20)

pyplot.savefig(path+'evolution_evo.png',bbox_inches='tight')

f=file(path+name+'_evolution1.dat','w')
dat=''.join([str(t[i])+','+str(rho22[i])+','+str(rho21_real[i])+','+str(rho21_imag[i])+'\n' for i in range(len(t))])
f.write(dat)
f.close()

################################################
#Now we will run the experiment many times varying the detuning to obtain
#a spectrum.
run_evolution(path,name+"_evolution",[E00],[-20.0],  2000,  0.01,Ne,spectrum_of_laser=1,N_delta=401,frequency_end=20.0,use_netcdf=use_netcdf)
#We read the results.
delta,rho22,rho21_real,rho21_imag = read_result(path,name+"_evolution",N=Ne,use_netcdf=use_netcdf)
################################################
#We plot and save the results for future reference.
fig=pyplot.figure(); ax=fig.add_subplot(111)
ax.plot(delta,rho22,'k-'		,label=r"$\rho_{22}$")
ax.plot(delta,rho21_real,'b-'	,label=r"$\mathfrak{Re}\rho_{21}$")
ax.plot(delta,rho21_imag,'r-'	,label=r"$\mathfrak{Im}\rho_{21}$")
ax.legend(loc=0,fontsize=17)
ax.set_xlabel(r"$\delta$",fontsize=20)

ax.set_xlim([-20,20])
pyplot.savefig(path+name+'_spectrum_evolution.png',bbox_inches='tight')

f=file(path+name+'_evolution2.dat','w')
dat=''.join([str(t[i])+','+str(rho22[i])+','+str(rho21_real[i])+','+str(rho21_imag[i])+'\n' for i in range(len(delta))])
f.write(dat)
f.close()

########################################################################
#                        Steady State Spectra                          #
########################################################################
#We write the Fortran code of the experiment.
write_stationary(path,name+"_steady",lasers,omega,gamma,r,Lij)
#We compile it.
compile_code(path,name+"_steady",lapack=True, parallel=parallel)

################################################
#We give the detuning and electric field amplitude initialy in absolute units.
delta=-20; E00=1.0

################################################
#We run the spectroscopy with the parameters we just calculated.
run_stationary(path,name+"_steady",[E00],[delta],spectrum_of_laser=1,N_delta=401,frequency_end=20.0,use_netcdf=use_netcdf)
#We read the results.
delta,rho22,rho21_real,rho21_imag = read_result(path,name+"_steady",N=Ne,use_netcdf=use_netcdf)

################################################
#We plot and save the results.
fig=pyplot.figure(); ax = fig.add_subplot(111)
ax.plot(delta,rho22,'k-'		,label=r"$\rho_{22}$")
ax.plot(delta,rho21_real,'b-'	,label=r"$\mathfrak{Re}\rho_{21}$")
ax.plot(delta,rho21_imag,'r-'	,label=r"$\mathfrak{Im}\rho_{21}$")
ax.legend(loc=0,fontsize=17)
ax.set_xlabel(r"$\delta$",fontsize=20)

ax.set_xlim([-20,20])
pyplot.savefig(path+name+'_spectrum_steady.png',bbox_inches='tight')

f=file(path+name+'_stationary2.dat','w')
dat=''.join([str(delta[i])+','+str(rho22[i])+','+str(rho21_real[i])+','+str(rho21_imag[i])+'\n' for i in range(len(delta))])
f.write(dat)
f.close()

########################################################################
# We now demonstrate power broadening.
delta1=-20.0

fig=pyplot.figure(); ax = fig.add_subplot(111)
E0=[0.25, 0.5, 1.0, 3.0, 7.0]
colors=['m','b','g',"orange",'r']
for i in range(len(E0)):
	
	run_stationary(path,name+"_steady",[E0[i]],[delta1],spectrum_of_laser=1,N_delta=401,frequency_end=20.0,use_netcdf=use_netcdf)
	delta,rho22,rho21_real,rho21_imag = read_result(path,name+"_steady",N=Ne,use_netcdf=use_netcdf)

	ax.plot(delta,rho22,'-',color=colors[i],label=r"$E_0="+str(E0[i])+"$")
	
	Omega=E0[i]
	gamma=1.0
	#~ hwfm=sqrt(gamma**2+ (Omega)**2)
	
	#~ ax.plot([-hwfm,hwfm],[rho22[200]/2,rho22[200]/2],color=colors[i])
	#~ a=0.025
	#~ ax.plot([-hwfm,-hwfm],[rho22[200]/2+a*rho22[200],rho22[200]/2-a*rho22[200]],color=colors[i])
	#~ ax.plot([ hwfm, hwfm],[rho22[200]/2+a*rho22[200],rho22[200]/2-a*rho22[200]],color=colors[i])

ax.legend(loc=0,fontsize=18)
ax.set_xlim([-20,20])
ax.set_xlabel(r"$\delta$",fontsize=22)
ax.set_ylabel(r"$\rho_{22}$",fontsize=22)

pyplot.savefig(path+name+'_power_broadening.png',bbox_inches='tight')
