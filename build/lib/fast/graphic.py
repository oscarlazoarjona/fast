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

sage_included = 'sage' in globals().keys()

if not sage_included:
	from math import atan2,sqrt,pi,cos,sin,exp
	from atomic_structure import find_fine_states, split_hyperfine_to_magnetic, make_list_of_states
	from atomic_structure import calculate_boundaries
	from misc import read_result

from math import pi
from colorsys import hls_to_rgb,hsv_to_rgb
from scipy.optimize import curve_fit
from matplotlib import pyplot
import numpy as np
from time import time
import os
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from electric_field import PlaneWave, MotField

def complex_to_color(z):
	if sage_included:
		if imag(z)==0 and real(z)==0:
			return [(0,0,0),0]

		h=(atan2(imag(z),real(z)) + pi)/(2*pi)
		l=sqrt(imag(z)**2+real(z)**2)
	else:
		if z.imag==0 and z.real==0:
			return [(0,0,0),0]

		h=(atan2(z.imag,z.real) + pi)/(2*pi)
		l=sqrt(z.imag**2+z.real**2)
	#print z,h,l
	h=float(h); l=float(l)
	if l> 1.0: return [(1,1,1),l]

	try:
		rgb=hls_to_rgb(h, l/2, 1.0)
	except:
		print h,l,1.0
		print type(h),type(l)
		rgb=hls_to_rgb(h, l, 1.0)
	return [rgb,l]

def complex_matrix_plot(A,logA=False,normalize=False,plot=True,**kwds):
	N=len(A[0])
	if logA:
		Anew=[]
		for i in range(N):
			row=[]
			for j in range(N):
				if A[i][j]!=0:
					row+=[log(log(A[i][j]))]
				else:
					row+=[0.0]
			Anew+=[row]
		A=Anew[:]
		#A=[[log(A[i][j]) for j in range(N)] for i in range(N)]

	if normalize:
		norm=1
		for i in range(N):
			for j in range(N):
				if abs(A[i][j])>norm: norm=abs(A[i][j])

		A=[[A[i][j]/norm for j in range(N)]for i in range(N)]

	#print A
	color_matrix=[]
	lmax=-1
	for i in range(N):
		row=[]
		for j in range(N):
			rgb,l=complex_to_color(A[i][j])
			row+=[rgb]
			if l> lmax:
				lmax=l
		color_matrix+=[row]

	if normalize:
		color_matrix=[[tuple([k/lmax for k in color_matrix[i][j]])  for j in range(N)] for i in range(N)]

	if plot:
		pyplot.imshow(color_matrix,interpolation='none')
		pyplot.savefig('a.png',bbox_inches='tight')
		pyplot.close('all')
	else:
		return color_matrix

def plot_disc(size=100):
	N=float(size)
	mat=[[ j/N + 1j*i/N for j in range(-size,size+1)] for i in reversed(range(-size,size+1))]
	mat=complex_matrix_plot(mat)
	pyplot.close('all')
	pyplot.imshow(mat,extent=[-1,1,-1,1])

	N=101; theta_step=2*pi/(N-1)
	theta=[i*theta_step for i in range(N)]
	x=[cos(theta[i]) for i in range(N)]
	y=[sin(theta[i]) for i in range(N)]

	pyplot.plot(x,y,'k-')
	pyplot.savefig('a.png',bbox_inches='tight')
	pyplot.close()

from matplotlib.colors import LogNorm,Normalize
import matplotlib.ticker as ticker
from matplotlib.colors import LinearSegmentedColormap

aaa=0.15
cdict1 = {'red':   ((0.0, 0.0, 0.0),
                   (0.5, 0.0, aaa),
                   (1.0, 1.0, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.0, 1.0),
                   (0.5, aaa, 0.0),
                   (1.0, 0.0, 0.0))
        }

bbb=0.002

cdict1 = {'blue':   ((0.0    , 0.0, 0.0),
					(0.5-bbb, 0.0, 0.0),
					(0.5+bbb, 0.0, aaa),
					(1.0    , 1.0, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'red':   ((0.0    , 0.0, 1.0),
					(0.5-bbb, aaa, 0.0),
					(0.5+bbb, 0.0, 0.0),
					(1.0    , 0.0, 0.0))
        }

blue_red1 = LinearSegmentedColormap('BlueRed1', cdict1)

aaa=0.15
aaa=0.25+.125
bbb=0.002


cdict2 = {
		'red':     ((0.0    , 1.0, 1.0),
					(0.5-bbb, aaa, 0.0),
					(0.5+bbb, 0.0, aaa),
					(1.0    , 1.0, 1.0)),

		'green':   ((0.0    , 0.0, 0.0),
					(0.5-bbb, 0.0, 0.0),
					(0.5+bbb, 0.0, aaa),
					(1.0    , 1.0, 1.0)),

        'blue':    ((0.0    , 1.0, 1.0),
					(0.5-bbb, aaa, 0.0),
					(0.5+bbb, 0.0, 0.0),
					(1.0    , 0.0, 0.0)),

        }

yellow_purple1 = LinearSegmentedColormap('BlueRed1', cdict2)

def fmt(x, pos):
	if x==0.0: return r'$0$'
	a, b = '{:.1e}'.format(x).split('e')
	b = int(b)
	if b in [0,1]: return r'$'+str(x)+'$'
	return r'${} \times 10^{{{}}}$'.format(a, b)

def fmt_log(x, pos):
	if x==0.0: return r'$0$'
	a, b = '{:.1e}'.format(x).split('e')
	b = int(b)
	#if b in [0,1]: return r'$'+str(x)+'$'
	return r'${} \times 10^{{{}}}$'.format(a, b)

def fancy_matrix_plot(ax,mat,states=None,path='',name='default.png',y_labels=True,
		complex_matrix=False,take_abs=False,take_log=False,colorbar=False,hyperfine_labels=False,center_on_zero=False,**kwds):

	rmat=[mat[i] for i in reversed(range(len(mat[0])))]

	smallest=min([min([j for j in i ]) for i in rmat])
	largest =max([max([j for j in i ]) for i in rmat])

	if take_abs:
		#We make the matrix all positive.
		rmat=[[abs(i) for i in row] for row in rmat]
		if take_log:
			#We make zeros into the smallest.
			smallest=min([min([j for j in i if j != 0.0]) for i in rmat])
			largest =max([max([j for j in i if j != 0.0]) for i in rmat])
			for i in range(len(rmat)):
				for j in range(len(rmat[0])):
					if rmat[i][j]==0.0:
						rmat[i][j]=smallest

	if complex_matrix:
		rmat=complex_matrix_plot(rmat,normalize=True,plot=False)

	#fig, ax = pyplot.subplots()

	if take_log:
		ax.imshow(rmat,interpolation='none', norm=LogNorm(vmin=smallest, vmax=largest) ,**kwds)
	elif center_on_zero:
		largest=max([abs(smallest),largest])
		ax.imshow(rmat,interpolation='none', norm=Normalize(vmin=-largest, vmax=largest),**kwds)
	else:
		ax.imshow(rmat,interpolation='none',**kwds)
	pyplot.axis('off')
	Ne=len(mat[0])
	ax.set_xlim([-0.5,Ne-0.5])
	ax.set_ylim([-0.5,Ne-0.5])

	if states != None:
		fine_states=find_fine_states(states)
		hyperfine_states    =make_list_of_states(fine_states   ,'hyperfine',verbose=0)#split_fine_to_hyperfine(fine_states)
		full_magnetic_states=make_list_of_states(fine_states   ,'magnetic',verbose=0)#split_fine_to_magnetic(fine_states)
		index_list_fine,index_list_hyperfine=calculate_boundaries(fine_states,full_magnetic_states)

		i_fine=[i[0] for i in index_list_fine[1:]]
		i_hyperfine=[i[0] for i in index_list_hyperfine[1:] if i[0] not in i_fine]

		#We place lines.
		for i in i_fine:
			ax.plot([-0.5,Ne-0.5],[Ne-(i+0.5),Ne-(i+0.5)],'r-',linewidth=0.5)
			ax.plot([i-0.5,i-0.5],[-0.5,Ne-0.5],'r-',linewidth=0.5)
		for i in i_hyperfine:
			ax.plot([-0.5,Ne-0.5],[Ne-(i+0.5),Ne-(i+0.5)],'b-',linewidth=0.5)
			ax.plot([i-0.5,i-0.5],[-0.5,Ne-0.5],'b-',linewidth=0.5)

		#We place the axis labels.
		if hyperfine_labels:
			for i in range(len(hyperfine_states)):
				a=hyperfine_states[i]._latex_()[18:]

				y=Ne- (index_list_hyperfine[i][1]+index_list_hyperfine[i][0])/2 -1.2
				x=    (index_list_hyperfine[i][1]+index_list_hyperfine[i][0])/2 -0.5
				if y_labels: ax.text(-0.75,y,'$'+a+'$',horizontalalignment='right')
				ax.text( x,Ne-1,'$'+a+'$',rotation=45,verticalalignment='bottom')
		else:
			for i in range(len(fine_states)):
				a=fine_states[i]._latex_()[18:]

				y=Ne- (index_list_fine[i][1]+index_list_fine[i][0])/2 -1.2
				x=    (index_list_fine[i][1]+index_list_fine[i][0])/2 -3.5
				if y_labels: ax.text(-0.75,y,'$'+a+'$',horizontalalignment='right')
				ax.text( x,Ne-1,'$'+a+'$',rotation=45,verticalalignment='bottom')
	else:
		for i in range(1,Ne+1):
			#a=hyperfine_states[i]._latex_()
			y=i -1 #(index_list_hyperfine[i][1]+index_list_hyperfine[i][0])/2 -1.2
			x=i -1 #   (index_list_hyperfine[i][1]+index_list_hyperfine[i][0])/2 -0.5
			ax.text(-0.5,y,'$'+str(Ne-i+1)+'$',horizontalalignment='right')
			ax.text( x,Ne-0.5,'$'+str(i)+'$',verticalalignment='bottom')


	#pyplot.savefig(path+name,bbox_inches='tight')
	#~ if colorbar:
		#~ if take_log:
			#~ pyplot.colorbar(format=ticker.FuncFormatter(fmt_log))#,norm=Normalize(vmin=smallest, vmax=largest))
		#~ else:
			#~ pyplot.colorbar(format=ticker.FuncFormatter(fmt))#,norm=Normalize(vmin=smallest, vmax=largest))
	#pyplot.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
	pyplot.savefig(path+name,bbox_inches='tight')
	#pyplot.close('all')

def plot_Lij(ax,Lij,Nl,states=None,path='',name='default.png',**kwds):
	Ne=len(Lij)

	mats=[]
	colors=[(0,0,1),(1,0,0),(0,1,0),(1,0.5,0)]
	mat=[]
	for i in range(Ne):
		row=[]
		for j in range(Ne):
			band=False
			for l in range(Nl):
				if l+1 in Lij[i][j]:
					row+=[colors[l]]
					band=True
					break
			if not band:
				row+=[(0,0,0)]

		mat+=[row]

	fancy_matrix_plot(ax,mat,states,path,name)

def fancy_r_plot(r,states=None,path='',name='default.png',y_labels=True,
		complex_matrix=False,take_abs=False,take_log=False,hyperfine_labels=False,**kwds):

	nn=3
	size=5
	corr=0.95
	Ne=len(r[0][0])
	f, plot = pyplot.subplots(1, nn, sharey=True,figsize=[3*size,size*corr])
	for p in range(nn):
		mat=r[p]
		rmat=[mat[i] for i in reversed(range(len(mat[0])))]
		if take_abs:
			rmat=[[abs(i) for i in row] for row in rmat]
		if complex_matrix:
			rmat=complex_matrix_plot(rmat,normalize=True,plot=False)

		ax=plot[p]
		ax.imshow(rmat,interpolation='none',**kwds)
		ax.axis('off')

		ax.set_xlim([-0.5,Ne])
		ax.set_ylim([-0.5,Ne-0.5])
		ax.text(Ne/2.0-1,-0.5,'$p='+str(p-1)+'$',verticalalignment='top')

		if states != None:
			fine_states=find_fine_states(states)
			hyperfine_states=make_list_of_states(fine_states,'hyperfine',verbose=0)#split_fine_to_hyperfine(fine_states)
			full_magnetic_states=make_list_of_states(fine_states,'magnetic',verbose=0)
			index_list_fine,index_list_hyperfine=calculate_boundaries(fine_states,full_magnetic_states)

			i_fine=[i[0] for i in index_list_fine[1:]]
			i_hyperfine=[i[0] for i in index_list_hyperfine[1:] if i[0] not in i_fine]

			for i in i_fine:
				ax.plot([-0.5,Ne-0.5],[Ne-(i+0.5),Ne-(i+0.5)],'r-')
				ax.plot([i-0.5,i-0.5],[-0.5,Ne-0.5],'r-')
			for i in i_hyperfine:
				ax.plot([-0.5,Ne-0.5],[Ne-(i+0.5),Ne-(i+0.5)],'b-')
				ax.plot([i-0.5,i-0.5],[-0.5,Ne-0.5],'b-')

			if hyperfine_labels:
				for i in range(len(hyperfine_states)):
					a=hyperfine_states[i]._latex_()
					y=Ne- (index_list_hyperfine[i][1]+index_list_hyperfine[i][0])/2 -1.2
					x=    (index_list_hyperfine[i][1]+index_list_hyperfine[i][0])/2 -0.5
					if p==0: ax.text(-0.75,y,'$'+a+'$',horizontalalignment='right')
					ax.text( x,Ne-1,'$'+a+'$',rotation=45,verticalalignment='bottom')
			else:
				for i in range(len(fine_states)):
					a=fine_states[i]._latex_()[18:]
					y=Ne- (index_list_fine[i][1]+index_list_fine[i][0])/2 -1.2
					x=    (index_list_fine[i][1]+index_list_fine[i][0])/2 -3.5
					if p==0: ax.text(-0.75,y,'$'+a+'$',horizontalalignment='right')
					ax.text( x,Ne-1,'$'+a+'$',rotation=45,verticalalignment='bottom')

	f.subplots_adjust(wspace=0)
	pyplot.savefig(path+name,bbox_inches='tight')
	#pyplot.close()

def make_video(path,name,Ne,states=None,duration=120,fps=6,digs=6,**kwds):
	data=read_result(path,name)
	Nt=len(data)
	for t in range(Nt):
		dati=data[t]

		mati=[[0.0 for j in range(Ne)] for i in range(Ne)]
		for i in range(2,Ne+1):
			for j in range(1,i+1):
				mu=Mu(i,j,1,Ne)
				rho_mu=dati[mu]
				mati[i-1][j-1]=rho_mu
				mati[j-1][i-1]=rho_mu
				if i!=j:
					mu=Mu(i,j,-1,Ne)
					rho_mu=dati[mu]
					mati[i-1][j-1]=mati[i-1][j-1]+1j*rho_mu
					mati[j-1][i-1]=mati[j-1][i-1]-1j*rho_mu
		mati[0][0]=1-sum([mati[i][i] for i in range(1,Ne)])
		n=str(t)
		fancy_matrix_plot(mati,states=states,path=path,name=name+'0'*(digs-len(n))+n,complex_matrix=True,**kwds)


	try:
		os.system('rm '+path+name+'.mp4')
	except:
		pass

	com='avconv -r '+str(fps)+' -i '+path+name+'%0'+str(digs)+'d.png -b:v 1000k '+path+name+'.mp4'
	print com
	os.system(com)

def fit_lorentizan(curve,p0=None,N_points=1000):
	'''Fits a lorentzian curve using p0=[x0,A,gamma] as an initial guess.
	It returns a curve with N_points.'''
	def lorentzian(x,x0,A,gamma): return A*gamma**2/((x-x0)**2+gamma**2)
	N=len(curve)
	x=[curve[i][0] for i in range(N)]
	y=[curve[i][1] for i in range(N)]

	from scipy.optimize import curve_fit
	popt,pcov = curve_fit(lorentzian,x,y,p0=p0)
	x0=popt[0]; A=popt[1]; gamma=popt[2]

	a=x[0]; b=x[-1]
	x_step=(b-a)/(N_points-1)
	x_fit=[a+i*x_step for i in range(N_points)]
	fited_curve=[(xi,lorentzian(xi,x0,A,gamma)) for xi in x_fit]
	return fited_curve,A,x0,gamma

def fit_lorentizan_with_background(curve,p0=None,N_points=1000):
	'''Fits a lorentzian curve using p0=[x0,A,gamma] as an initial guess.
	It returns a curve with N_points.'''
	def lorentzian(x,x0,A,gamma,B): return A*gamma**2/((x-x0)**2+gamma**2)+B
	N=len(curve)
	x=[curve[i][0] for i in range(N)]
	y=[curve[i][1] for i in range(N)]

	from scipy.optimize import curve_fit
	popt,pcov = curve_fit(lorentzian,x,y,p0=p0)
	x0=popt[0]; A=popt[1]; gamma=popt[2]; B=popt[3]

	a=x[0]; b=x[-1]
	x_step=(b-a)/(N_points-1)
	x_fit=[a+i*x_step for i in range(N_points)]
	fited_curve=[(xi,lorentzian(xi,x0,A,gamma,B)) for xi in x_fit]
	return fited_curve,A,x0,gamma,B

########################################################################
#Drawing 3D beam diagrams.
########################################################################

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

def bar_chart_mf(data,path_name):
	N=len(data)

	ind = np.arange(N)  # the x locations for the groups
	width = 0.8       # the width of the bars

	fig, ax = pyplot.subplots()
	rects1 = ax.bar(ind, data, width, color='g')

	# add some text for labels, title and axes ticks
	ax.set_ylabel('Population')
	ax.set_xticks(ind+width/2)
	labs=['m='+str(i) for i in range(-N/2+1,N/2+1)]
	ax.set_xticklabels( labs )

	def autolabel(rects):
		# attach some text labels
		for rect in rects:
			height = rect.get_height()

	autolabel(rects1)
	pyplot.savefig(path_name)
	pyplot.close()

def draw_atom3d(ax):
	#fig = pyplot.figure()
	#ax = fig.gca(projection='3d')

	ax.plot([0],[0],[0],'k.')
	Nt=100
	s_step=2*pi/(Nt-1)
	v1=[[],[],[]]; v2=[[],[],[]]; v3=[[],[],[]]
	for i in range(Nt):
		s=i*s_step
		v1[0]+=[0.250000000000000*cos(s)]
		v1[1]+=[0.250000000000000*sin(s)]
		v1[2]+=[0]
		v2[0]+=[0.176776695296637*cos(s) - 0.125000000000000*sin(s)]
		v2[1]+=[ -0.176776695296637*cos(s) - 0.125000000000000*sin(s)]
		v2[2]+=[-0.176776695296637*sin(s)]
		v3[0]+=[-0.176776695296637*cos(s) + 0.125000000000000*sin(s)]
		v3[1]+=[ 0.176776695296637*cos(s) + 0.125000000000000*sin(s)]
		v3[2]+=[-0.176776695296637*sin(s)]


	ax.plot(v1[0],v1[1],v1[2],'k-')
	ax.plot(v2[0],v2[1],v2[2],'k-')
	ax.plot(v3[0],v3[1],v3[2],'k-')

def draw_plane_wave_3d(ax,beam,dist_to_center=0):
	Ex=[]; Ey=[]; Ez=[]

	k=[cos(beam.phi)*sin(beam.theta),sin(beam.phi)*sin(beam.theta),cos(beam.theta)]
	kx,ky,kz=k

	Nt=1000
	tstep=7*pi/4/(Nt-1)

	alpha=beam.alpha
	beta =beam.beta
	phi  =beam.phi
	theta=beam.theta
	omega=1

	for i in range(Nt):
		t=i*tstep
		#~ Ex+=[(cos(beam.alpha)*cos(beam.phi)*cos(beam.theta) + sin(beam.alpha)*sin(beam.phi))*cos(beam.beta + t) -(cos(beam.phi)*cos(beam.theta)*sin(beam.alpha) - cos(beam.alpha)*sin(beam.phi))*cos(t)-dist_to_center*kx]
		#~ Ey+=[(cos(beam.alpha)*cos(beam.theta)*sin(beam.phi) - cos(beam.phi)*sin(beam.alpha))*cos(beam.beta + t) -(cos(beam.theta)*sin(beam.alpha)*sin(beam.phi) + cos(beam.alpha)*cos(beam.phi))*cos(t)-dist_to_center*ky]
		#~ Ez+=[-cos(beam.alpha)*cos(beam.beta + t)*sin(beam.theta) + cos(t)*sin(beam.alpha)*sin(beam.theta)-dist_to_center*kz]

		Ex+=[(cos(2*alpha)*cos(phi)*cos(theta) - sin(2*alpha)*sin(phi))*cos(omega*t)*cos(2*beta) - (cos(phi)*cos(theta)*sin(2*alpha) + cos(2*alpha)*sin(phi))*sin(omega*t)*sin(2*beta) - dist_to_center*kx]

		Ey+=[(cos(2*alpha)*cos(theta)*sin(phi) + cos(phi)*sin(2*alpha))*cos(omega*t)*cos(2*beta) - (cos(theta)*sin(2*alpha)*sin(phi) - cos(2*alpha)*cos(phi))*sin(omega*t)*sin(2*beta) - dist_to_center*ky]

		Ez+=[-cos(omega*t)*cos(2*alpha)*cos(2*beta)*sin(theta) + sin(omega*t)*sin(2*alpha)*sin(2*beta)*sin(theta) - dist_to_center*kz]


	ax.plot(Ex,Ey,Ez,beam.color+'-')
	ff=dist_to_center-1.0

	arrx=[-kx*dist_to_center,-kx*ff]
	arry=[-ky*dist_to_center,-ky*ff]
	arrz=[-kz*dist_to_center,-kz*ff]
	arrow=Arrow3D(arrx,arry,arrz, mutation_scale=20, lw=1, arrowstyle="-|>", color=beam.color)
	ax.add_artist(arrow)
	ax.plot([Ex[-1]],[Ey[-1]],[Ez[-1]],'.',markersize=8,color=beam.color)

def draw_mot_field_3d(ax,mot_field,dist_to_center=0):
	for beam in mot_field.beams:
		draw_plane_wave_3d(ax,beam,dist_to_center=dist_to_center)

def draw_lasers_3d(ax,lasers,name='default.png',distances=None,lim=None):

	draw_atom3d(ax)

	if distances==None: distances=[1.0 for i in range(len(lasers))]
	for i in range(len(lasers)):
		if type(lasers[i])==PlaneWave:
			draw_plane_wave_3d( ax,lasers[i],distances[i])
		elif type(lasers[i])==MotField:
			draw_mot_field_3d(  ax,lasers[i],distances[i])

	ax.set_xlabel(r"$x$",fontsize=20)
	ax.set_ylabel(r"$y$",fontsize=20)
	ax.set_zlabel(r"$z$",fontsize=20)

	if lim==None: lim=sqrt(2.0)


	ax.set_xlim(-lim, lim)
	ax.set_ylim(-lim, lim)
	ax.set_zlim(-lim, lim)
	ax.set_aspect("equal")

	pyplot.savefig(name,bbox_inches='tight')


def plot_eigenvalues(path,name,Ne,Omega=1,filename='a.png'):
	times=characteristic_times(path,name,Omega=Omega)

	Nd=Ne**2-1
	for i in range(Nd):
		if float('inf') in times[1+i]:
			#print i,times[1+i],times[1+Nd+i],'asd'
			pass
		else:
			#print i,times[1+i],times[1+Nd+i],'ert'
			if len(times[1+Nd+i])>1:
				pyplot.loglog(times[1+Nd+i],times[1+i])
			else:
				pyplot.loglog(times[1+Nd+i],times[1+i],'+')
	pyplot.savefig(filename,bbox_inches='tight')
	pyplot.close('all')

def plot_populations(path,name,Ne,states=None,filename='a.png',fontsize=12,absolute_frequency=True,
					save_path='',use_netcdf=True):
	pyplot.close("all")
	dat=read_result(path,name,N=Ne,use_netcdf=use_netcdf)
	x=dat[0]
	if absolute_frequency:
		x=[xi/2/pi for xi in x]
	pop=dat[1:Ne]
	Nd=len(pop[0]);	Nr=Ne**2-1

	pop1=[1-sum([pop[j][i] for j in range(Ne-1)]) for i in range(Nd)]
	pop =[pop1]+pop

	#We do different things depending on what states we are given.

	if states==None:
		#If we recieve no states all populations are ploted.
		for i in range(Ne):
			pyplot.plot(x,pop[i],label=r"$\mathrm{Poblaci\'on} \ " +str(i+1)+"$")
		pyplot.legend(fontsize=fontsize)
		pyplot.savefig(save_path+filename,bbox_inches='tight')
		pyplot.close('all')

	elif len(states[0].quantum_numbers)>=5:
		#If we recieve magnetic states we make a plot for each hyperfine state.
		magnetic_plots=len(states[0].quantum_numbers)==6

		if not magnetic_plots:
			N_plots=len(states)
			states=split_hyperfine_to_magnetic(states)
			hyperfine_colors=list(reversed([hsv_to_rgb(m*0.8/(N_plots-1),1.0,1.0) for m in range(N_plots)]))
			conta=0

		fine_states=find_fine_states(states)
		boundaries=list(reversed(calculate_boundaries(fine_states,states)[1]))

		for pair in boundaries:
			f=states[pair[0]].f

			if f==0:
				colors=[(0.8,0.0,1.0)]
			else:
				colors=list(reversed([hsv_to_rgb(m*0.8/f,1.0,1.0) for m in range(f+1)]))

			for i in range(pair[0],pair[1]):
				m=states[i].m
				if m<0:
					color=colors[-m]
					style=':'
				else:
					color=colors[m]
					style='-'
				if magnetic_plots:
					pyplot.plot(x,pop[i],style,label=r"$\mathrm{Poblaci\'on} \ M_F=" +str(states[i].m)+"$",color=color)

			if magnetic_plots:
				if f!=0:
					suma=[ sum([ pop[i][j] for i in range(pair[0],pair[1])]) for j in range(len(pop[0]))]
					pyplot.plot(x,suma,'k-',label=r'$\mathrm{suma}$')

				filenamei=str(states[i]).split()[1].replace('_','').replace('/','_').replace('^','F=')
				s=filenamei.find(',')
				filenamei=filenamei[:s]
				filenamei=name+'_'+filenamei+'.png'

				title=find_fine_states([states[i]])[0]._latex_()+'\ F='+str(states[i].f)

				pyplot.title(r"$"+title+"$")
				pyplot.ylim([0,None])
				pyplot.xlim([x[0],x[-1]])
				pyplot.legend(fontsize=fontsize,loc=0)

				pyplot.savefig(save_path+filenamei,bbox_inches='tight')
				pyplot.close('all')
			else:
				suma=[ sum([ pop[i][j] for i in range(pair[0],pair[1])]) for j in range(len(pop[0]))]
				label=states[i]._latex_()
				label=label[label.find(' ')+1:]
				label=label[:label.find('^')]
				label+=r"\ F="+str(states[i].f)
				label ="$"+label+"$"

				pyplot.plot(x,suma,'-',color=hyperfine_colors[conta],label=label)
				#for i in range(len(x)): print i,x[i],suma[i]
				conta+=1
				#pyplot.savefig('a.png')

		if not magnetic_plots:
			title=states[0]._latex_()
			title="$"+title[:title.find(' ')-1]+"$"
			pyplot.title(title,fontsize=20)
			pyplot.xlim([x[0],x[-1]])
			pyplot.legend(fontsize=fontsize,loc=0)
			pyplot.savefig(save_path+name+'_pops.png',bbox_inches='tight')
			pyplot.close('all')

########################################################################
#Drawing of experiment diagrams.
########################################################################
from matplotlib import pyplot
import numpy as np
from math import sin, cos, sqrt, atan2, pi

def rotate_and_traslate(cur,alpha,v0):
	if len(cur)>2 or (type(cur[0][0]) in [list,tuple]):
		cur_list=cur[:]
		for i in range(len(cur_list)):
			curi=cur_list[i]
			curi=rotate_and_traslate(curi,alpha,v0)
			cur_list[i]=curi
		return cur_list
	else:
		x0,y0=cur
		#print v0,111
		rot=np.matrix([[cos(alpha),-sin(alpha)],[sin(alpha),cos(alpha)]])
		xn=[]; yn=[]
		for i in range(len(x0)):
			v=np.matrix([[x0[i]],[y0[i]]])

			vi=np.dot(rot,v)
			xn+=[float(vi[0][0])+v0[0]]; yn+=[float(vi[1][0])+v0[1]]

		return xn,yn

def mirror(ax,p0,alpha=0,size=2.54,width=0.5,format=None):
	if format==None: format='k-'

	x0=[size/2,-size/2, -size/2,  size/2,size/2]
	y0=[0   ,    0,-width,-width,   0]

	x1=[size/2,size/2-width];                             y1=[0,-width]
	x2=[-size/2+width,-size/2];                           y2=[0,-width]
	x3=[(size/2-size/2+width)/2,(size/2-width-size/2)/2]; y3=[0,-width]

	cur_list=[(x0,y0),(x1,y1),(x2,y2),(x3,y3)]
	cur_list=rotate_and_traslate(cur_list,alpha,p0)
	for curi in cur_list: ax.plot(curi[0],curi[1],format)

def vacuum_chamber(ax,p0,size=4.5,alpha=0,format=None):
	if format==None: format='k-'
	a=size
	b=2.5*size/4.5
	r=sqrt(a**2+b**2)

	x0=[a,a]; y0=[-b,b]

	cur_list0=[(x0,y0)]
	for i in range(1,4):
		cur_list0+=[rotate_and_traslate(cur_list0[0],pi/2*i,(0,0))]
	cur_list0=rotate_and_traslate(cur_list0,alpha,p0)


	N=100
	ang=atan2(b,a)
	ang0=ang; angf=pi/2-ang
	angstep=(angf-ang0)/(N-1)
	x1=[ r*cos(i*angstep+ang0) for i in range(N)]
	y1=[ r*sin(i*angstep+ang0) for i in range(N)]

	cur_list=[(x1,y1)]
	for i in range(1,4):
		cur_list+=[rotate_and_traslate(cur_list[0],pi/2*i,(0,0))]
	cur_list=rotate_and_traslate(cur_list,alpha,p0)

	for curi in cur_list0: ax.plot(curi[0],curi[1],'c-')
	for curi in cur_list : ax.plot(curi[0],curi[1],format)

def cloud(ax,p0,size=1.0,alpha=0,format=None,**kwds):
	if format==None: format='k-'
	size=size/4.0
	n=8
	ang=2*pi/n

	N=100
	ang0=-pi/4-0.4; angf= pi/4+0.4
	angstep=(angf-ang0)/(N-1)
	x1=[ size*(cos(angstep*i+ang0)+1) for i in range(N)]
	y1=[ size*sin(angstep*i+ang0) for i in range(N)]

	cur_list=[(x1,y1)]
	for i in range(1,n):
		cur_list+=[rotate_and_traslate(cur_list[0],ang*i,(0,0))]

	cur_list=rotate_and_traslate(cur_list,alpha,p0)
	for curi in cur_list : ax.plot(curi[0],curi[1],format,**kwds)

def lens(ax,p0,size,focus,format=None,join_with_focus=False,**kwds):
	if format==None: format='k-'
	alpha=atan2(p0[1]-focus[1],p0[0]-focus[0])+pi/2

	f=sqrt((p0[0]-focus[0])**2 + (p0[1]-focus[1])**2)
	r=sqrt((size/2)**2 + f**2)
	ang=atan2(size/2,f)

	N=100
	ang0=-ang; angf= ang
	angstep=(angf-ang0)/(N-1)
	x1=[ r*sin(i*angstep+ang0) for i in range(N)]
	y1=[ r*cos(i*angstep+ang0)-f for i in range(N)]
	y2=[ -y1[i] for i in range(N)]


	cur_list=[(x1,y1),(x1,y2)]
	#print cur_list[0]
	cur_list=rotate_and_traslate(cur_list,alpha,p0)
	for curi in cur_list : ax.plot(curi[0],curi[1],format,**kwds)

	if join_with_focus:
		x1,y1=[p0[0]+size/2*cos(alpha),p0[1]+size/2*sin(alpha)]
		x2,y2=[p0[0]-size/2*cos(alpha),p0[1]-size/2*sin(alpha)]
		ax.plot([x1,focus[0]],[y1,focus[1]],':',**kwds)
		ax.plot([x2,focus[0]],[y2,focus[1]],':',**kwds)

def eye(ax,p0,size=1.0,alpha=0,format=None,**kwds):
	if format==None: format='k-'

	N=100
	ang0=pi-3*pi/16; angf= pi+3*pi/16
	angstep=(angf-ang0)/(N-1)
	x1=[ size*(cos(i*angstep+ang0)+1) for i in range(N)]
	y1=[ size*sin(i*angstep+ang0) for i in range(N)]

	ang2=ang0+pi/16
	x2=[size,     size*(1.2*cos(ang2) +1 )]
	y2=[0,        1.2*size*(sin(ang2)  )]
	y3=[0,       -1.2*size*(sin(ang2)  )]

	N=100
	ang0=ang2; angf= ang2+4*pi/16
	angstep=(angf-ang0)/(N-1)
	x4=[ size*(0.85*cos(i*angstep+ang0)+1) for i in range(N)]
	y4=[ size*0.85*sin(i*angstep+ang0) for i in range(N)]


	cur_list=[(x1,y1),(x2,y2),(x2,y3),(x4,y4)]
	cur_list=rotate_and_traslate(cur_list,alpha,p0)

	for curi in cur_list : ax.plot(curi[0],curi[1],format,**kwds)

def beam_splitter(ax,p0,size=2.54,alpha=0,format=None,**kwds):
	if format==None: format='k-'
	a=size/2
	x0=[ a,-a,-a, a, a,-a]
	y0=[ a, a,-a,-a, a,-a]

	cur_list=[(x0,y0)]
	cur_list=rotate_and_traslate(cur_list,alpha,p0)

	for curi in cur_list : ax.plot(curi[0],curi[1],format,**kwds)

def draw_beam(ax,p1,p2,width=0,beta1=None,beta2=None,format=None,**kwds):
	if format==None: format='k-'

	if width==0:
		x0=[p1[0],p2[0]]
		y0=[p1[1],p2[1]]
		ax.plot(x0,y0,format,**kwds)
	else:
		alpha=atan2(p2[1]-p1[1],p2[0]-p1[0])
		a=width/2

		x1,y1=p1
		x2,y2=p2

		x11 = (a*x1**2*cos(beta1) - 2*a*x1*x2*cos(beta1) + a*x2**2*cos(beta1) + a*y1**2*cos(beta1) + a*y2**2*cos(beta1) - (2*a*y1*cos(beta1) - sqrt((x1 - x2)**2 + (y1 - y2)**2)*x1*cos(beta1))*y2 - (x1*y1*cos(beta1) - x1**2*sin(beta1) + x1*x2*sin(beta1))*sqrt((x1 - x2)**2 + (y1 - y2)**2))/(sqrt((x1 - x2)**2 + (y1 - y2)**2)*y2*cos(beta1) - sqrt((x1 - x2)**2 + (y1 - y2)**2)*(y1*cos(beta1) - x1*sin(beta1) + x2*sin(beta1)))
		y11 = (a*x1**2*sin(beta1) - 2*a*x1*x2*sin(beta1) + a*x2**2*sin(beta1) + a*y1**2*sin(beta1) + a*y2**2*sin(beta1) - (2*a*y1*sin(beta1) - sqrt((x1 - x2)**2 + (y1 - y2)**2)*y1*cos(beta1))*y2 - (y1**2*cos(beta1) - (x1*sin(beta1) - x2*sin(beta1))*y1)*sqrt((x1 - x2)**2 + (y1 - y2)**2))/(sqrt((x1 - x2)**2 + (y1 - y2)**2)*y2*cos(beta1) - sqrt((x1 - x2)**2 + (y1 - y2)**2)*(y1*cos(beta1) - x1*sin(beta1) + x2*sin(beta1)))

		x21 = (a*x1**2*cos(beta2) - 2*a*x1*x2*cos(beta2) + a*x2**2*cos(beta2) + a*y1**2*cos(beta2) + a*y2**2*cos(beta2) - (2*a*y1*cos(beta2) - sqrt((x1 - x2)**2 + (y1 - y2)**2)*x2*cos(beta2))*y2 - (x2*y1*cos(beta2) - x1*x2*sin(beta2) + x2**2*sin(beta2))*sqrt((x1 - x2)**2 + (y1 - y2)**2))/(sqrt((x1 - x2)**2 + (y1 - y2)**2)*y2*cos(beta2) - sqrt((x1 - x2)**2 + (y1 - y2)**2)*(y1*cos(beta2) - x1*sin(beta2) + x2*sin(beta2)))
		y21 = (a*x1**2*sin(beta2) - 2*a*x1*x2*sin(beta2) + a*x2**2*sin(beta2) + a*y1**2*sin(beta2) + (a*sin(beta2) + sqrt((x1 - x2)**2 + (y1 - y2)**2)*cos(beta2))*y2**2 - (2*a*y1*sin(beta2) + sqrt((x1 - x2)**2 + (y1 - y2)**2)*(y1*cos(beta2) - x1*sin(beta2) + x2*sin(beta2)))*y2)/(sqrt((x1 - x2)**2 + (y1 - y2)**2)*y2*cos(beta2) - sqrt((x1 - x2)**2 + (y1 - y2)**2)*(y1*cos(beta2) - x1*sin(beta2) + x2*sin(beta2)))

		ax.plot([x11,x21],[y11,y21],format,**kwds)

		x12 = -(a*x1**2*cos(beta2) - 2*a*x1*x2*cos(beta2) + a*x2**2*cos(beta2) + a*y1**2*cos(beta2) + a*y2**2*cos(beta2) - (2*a*y1*cos(beta2) + sqrt((x1 - x2)**2 + (y1 - y2)**2)*x2*cos(beta2))*y2 + (x2*y1*cos(beta2) - x1*x2*sin(beta2) + x2**2*sin(beta2))*sqrt((x1 - x2)**2 + (y1 - y2)**2))/(sqrt((x1 - x2)**2 + (y1 - y2)**2)*y2*cos(beta2) - sqrt((x1 - x2)**2 + (y1 - y2)**2)*(y1*cos(beta2) - x1*sin(beta2) + x2*sin(beta2)))
		y12 = -(a*x1**2*sin(beta2) - 2*a*x1*x2*sin(beta2) + a*x2**2*sin(beta2) + a*y1**2*sin(beta2) + (a*sin(beta2) - sqrt((x1 - x2)**2 + (y1 - y2)**2)*cos(beta2))*y2**2 - (2*a*y1*sin(beta2) - sqrt((x1 - x2)**2 + (y1 - y2)**2)*(y1*cos(beta2) - x1*sin(beta2) + x2*sin(beta2)))*y2)/(sqrt((x1 - x2)**2 + (y1 - y2)**2)*y2*cos(beta2) - sqrt((x1 - x2)**2 + (y1 - y2)**2)*(y1*cos(beta2) - x1*sin(beta2) + x2*sin(beta2)))

		x22 = -(a*x1**2*cos(beta1) - 2*a*x1*x2*cos(beta1) + a*x2**2*cos(beta1) + a*y1**2*cos(beta1) + a*y2**2*cos(beta1) - (2*a*y1*cos(beta1) + sqrt((x1 - x2)**2 + (y1 - y2)**2)*x1*cos(beta1))*y2 + (x1*y1*cos(beta1) - x1**2*sin(beta1) + x1*x2*sin(beta1))*sqrt((x1 - x2)**2 + (y1 - y2)**2))/(sqrt((x1 - x2)**2 + (y1 - y2)**2)*y2*cos(beta1) - sqrt((x1 - x2)**2 + (y1 - y2)**2)*(y1*cos(beta1) - x1*sin(beta1) + x2*sin(beta1)))
		y22 = -(a*x1**2*sin(beta1) - 2*a*x1*x2*sin(beta1) + a*x2**2*sin(beta1) + a*y1**2*sin(beta1) + a*y2**2*sin(beta1) - (2*a*y1*sin(beta1) + sqrt((x1 - x2)**2 + (y1 - y2)**2)*y1*cos(beta1))*y2 + (y1**2*cos(beta1) - (x1*sin(beta1) - x2*sin(beta1))*y1)*sqrt((x1 - x2)**2 + (y1 - y2)**2))/(sqrt((x1 - x2)**2 + (y1 - y2)**2)*y2*cos(beta1) - sqrt((x1 - x2)**2 + (y1 - y2)**2)*(y1*cos(beta1) - x1*sin(beta1) + x2*sin(beta1)))

		ax.plot([x12,x22],[y12,y22],format,**kwds)

def simple_beam_splitter(ax,p0,size=2.54,width=0.1,alpha=0,format=None,**kwds):
	if format==None: format='k-'
	a=size/2
	b=width/2
	x0=[ a,-a,-a, a,a]
	y0=[ b, b,-b,-b,b]

	cur_list=[(x0,y0)]
	cur_list=rotate_and_traslate(cur_list,alpha,p0)

	for curi in cur_list : ax.plot(curi[0],curi[1],format,**kwds)

def draw_laser(ax,p0,width,height,alpha=0,format=None,**kwds):
	if format==None: format='k-'
	simple_beam_splitter(ax,p0,size=width,width=height,alpha=alpha,format=format,**kwds)

def cable(ax,p1,p2,format=None,**kwds):
	if format==None: format='k-'
	pyplot.plot([p1[0],p2[0]],[p1[1],p2[1]],format,**kwds)

########################################################################
#Drawing DAQ.
########################################################################

def draw_box(ax,p0,width,height,alpha=0,format=None,**kwds):
	if format==None: format='k-'
	simple_beam_splitter(ax,p0,size=width,width=height,alpha=alpha,format=format,**kwds)

def draw_for(ax,p0,pf,r=10,format=None,**kwds):
	if format==None: format='k-'

	#We draw lines
	x0,y0=p0
	xf,yf=pf

	xtop=[x0+r,xf-r]
	ytop=[y0  ,y0  ]

	xbot=[x0+r,xf-r]
	ybot=[yf  ,yf  ]

	xl  =[x0  ,x0  ]
	yl  =[y0-r,yf+r]

	xr  =[xf  ,xf  ]
	yr  =[y0-r,yf+r]

	cur_list=[(xtop,ytop),(xbot,ybot),(xl,yl),(xr,yr)]
	cur_list=rotate_and_traslate(cur_list,0.0,(0,0))
	for curi in cur_list : ax.plot(curi[0],curi[1],format,**kwds)


	#We draw a round corner

	p1=[xf-r,y0-r]
	p2=[x0+r,y0-r]
	p3=[x0+r,yf+r]
	p4=[xf-r,yf+r]

	dtheta=pi/2/(100-1)
	theta=[i*dtheta for i in range(100+1)]
	xc=[r*cos(thetai) for thetai in theta]
	yc=[r*sin(thetai) for thetai in theta]

	cur_list1=[(xc,yc)]
	cur_list1=rotate_and_traslate(cur_list1,0.0,p1)
	for curi in cur_list1 : ax.plot(curi[0],curi[1],format,**kwds)

	cur_list2=[(xc,yc)]
	cur_list2=rotate_and_traslate(cur_list2,pi/2,p2)
	for curi in cur_list2 : ax.plot(curi[0],curi[1],format,**kwds)

	cur_list3=[(xc,yc)]
	cur_list3=rotate_and_traslate(cur_list3,pi  ,p3)
	for curi in cur_list3 : ax.plot(curi[0],curi[1],format,**kwds)

	cur_list4=[(xc,yc)]
	cur_list4=rotate_and_traslate(cur_list4,-pi/2,p4)
	for curi in cur_list4 : ax.plot(curi[0],curi[1],format,**kwds)

def draw_arith(ax,p0,size=1,alpha=0,arith=None,format=None,fontsize=10,**kwds):
	if format==None: format='k-'
	a=size/2.0
	x0=[ 0,2.5*a, 0, 0]
	y0=[ a,  0,-a, a]

	cur_list=[(x0,y0)]
	cur_list=rotate_and_traslate(cur_list,alpha,p0)

	for curi in cur_list : ax.plot(curi[0],curi[1],format,**kwds)

	if arith!=None:
		pyplot.text( p0[0]+0.75*a,p0[1],arith,horizontalalignment='center',verticalalignment='center',fontsize=fontsize)



########################################################################
#Drawing energy levels.
########################################################################

def draw_state(ax,p,text='',l=0.5,alignment='left',label_displacement=1.0,fontsize=25,**kwds):

    plo=ax.plot([p[0]-l/2.0,p[0]+l/2.0],[p[1],p[1]],color='black',**kwds)
    if text!='':
		if alignment=='left':
			ax.text(p[0] -l/2.0 - label_displacement,p[1],text,
					horizontalalignment='right',verticalalignment='center',
					color='black',fontsize=fontsize)
		elif alignment=='right':
			ax.text(p[0] +l/2.0 + label_displacement,p[1],text,horizontalalignment='left' ,color='black',fontsize=fontsize)

def draw_multiplet(ax,fine_state,p,hmin,w, fside='right',label_separation=1,label_fontsize=15,fsize=10,deltanu_fontsize=6,
					proportional=False,text='',text_pos='top',magnetic_lines=False,**kwds):
	#We determine the vertical positions, calculated from p[1] up.
	hyperfine_states=make_list_of_states([fine_state],'hyperfine')

	h_list=[ ei.nu -hyperfine_states[0].nu for ei in hyperfine_states]
	h_list=[ i/h_list[-1]      for i in h_list]
	h_min=min([h_list[i+1]-h_list[i] for i in range(len(h_list)-1)])
	h_list=[ hmin*i/h_min + p[1]  for i in h_list]

	if proportional:
		h_list=[p[1]+i*hmin for i in range(len(hyperfine_states))]

	omegaij=[ (hyperfine_states[i+1].nu-hyperfine_states[i].nu)/1e6 for i in range(len(hyperfine_states)-1)]

	for i in range(len(h_list)):
		label='$\mathrm{F}='+str(hyperfine_states[i].f)+'$'
		if magnetic_lines:
			maxf=max([eee.f for eee in hyperfine_states])
			f=hyperfine_states[i].f
			nm=2*maxf+1
			for mf in range(-f,f+1):
				draw_state(ax,[p[0]+mf*w/nm,h_list[i]],"",w/nm*0.5,alignment=fside,fontsize=fsize)

			if fside=='right':
				ax.text(p[0]+w+label_separation,h_list[i],label,fontsize=fsize,horizontalalignment="right",verticalalignment="center")
			elif fside=='left':
				ax.text(p[0]-w-label_separation,h_list[i],label,fontsize=fsize,horizontalalignment="left",verticalalignment="center")
		else:
			draw_state(ax,[p[0],h_list[i]],label,w,alignment=fside,fontsize=fsize)

	for i in range(len(h_list)-1):
		hmid=(h_list[i+1]+h_list[i])/2.0-0.5
		nu=str(omegaij[i])[:5]
		if fside=='left':
			ax.text(p[0]-w/2.0,hmid,r'$'+nu+' \ \mathrm{MHz}$',fontsize=deltanu_fontsize,
					horizontalalignment=fside,verticalalignment='bottom')
		else:
			ax.text(p[0]+w/2.0,hmid,r'$'+nu+' \ \mathrm{MHz}$',fontsize=deltanu_fontsize,
					horizontalalignment=fside,verticalalignment='bottom')

	a=label_separation

	if text !='':
		if text_pos=='top':
			labelx=p[0]
			labely=h_list[-1]+a

			ax.text(labelx,labely,'$'+text+'$',
					verticalalignment='bottom',
					horizontalalignment='center',fontsize=label_fontsize)
		elif text_pos=='right':
			labelx=p[0]+w/2+2.0*a
			if fside=='right': labelx=labelx+a*5.0
			labely=(h_list[-1]+h_list[0])/2.0

			ax.text(labelx,labely,'$'+text+'$',
					verticalalignment='center',
					horizontalalignment='left',fontsize=label_fontsize)

		elif text_pos=='left':
			labelx=p[0]-w/2-2.0*a
			if fside=='left': labelx=labelx-a*5.0
			labely=(h_list[-1]+h_list[0])/2.0

			ax.text(labelx,labely,'$'+text+'$',
					verticalalignment='center',
					horizontalalignment='right',fontsize=label_fontsize)

	return [[p[0],i] for i in h_list]

def excitation(ax,p1,p2,**kwargs):
	x1,y1=p1
	x2,y2=p2
	dx=p2[0]-x1
	dy=p2[1]-y1

	ax.arrow( x1, y1, dx, dy,length_includes_head=True,**kwargs)
	#ax.arrow( x2, y2,-dx,-dy,length_includes_head=True,**kwargs)

def decay(ax,p0,pf,A,n,format=None,**kwds):
	if format==None: format='k-'
	T=sqrt((p0[0]-pf[0])**2+(p0[1]-pf[1])**2)
	alpha=atan2(pf[1]-p0[1],pf[0]-p0[0])

	x=[ i*T/400.0 for i in range(401)]
	y=[A*sin(xi * 2*pi*n/T) for xi in x]

	cur_list=[(x,y)]
	cur_list=rotate_and_traslate(cur_list,alpha,p0)

	for curi in cur_list : ax.plot(curi[0],curi[1],format,**kwds)
