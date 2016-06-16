# -*- coding: utf-8 -*-

#This flag tells whether we are running in Sage or not.
sage_included = 'sage' in globals().keys()
if not sage_included:
	from math import pi,sqrt
else:
	from symbolic import define_density_matrix,define_substitutions

#We make all needed imports.
#~ pyplot_included = 'pyplot' in globals().keys()
#~ if not pyplot_included:
from matplotlib import use; use('Agg')# This allows plots to be made remotely via ssh.

from electric_field import PlaneWave, MotField
from electric_field import electric_field_amplitude_gaussian, electric_field_amplitude_top, electric_field_amplitude_intensity
from misc import Mu, IJ, find_phase_transformation
from misc import formatLij, convolve_with_gaussian, read_result

from graphic import complex_matrix_plot, plot_Lij
from graphic import Arrow3D,bar_chart_mf,draw_atom3d,draw_mot_field_3d,draw_plane_wave_3d,draw_lasers_3d
from graphic import draw_state, excitation, decay, draw_multiplet
from graphic import fancy_matrix_plot, fancy_r_plot, plot_populations

from evolution import write_evolution, run_evolution
from stationary import write_stationary, run_stationary
from misc import compile_code

from atomic_structure_Rb import S,P,D,F,State,Transition
from atomic_structure_Rb import split_fine_to_hyperfine, split_fine_to_magnetic, split_hyperfine_to_magnetic
from atomic_structure_Rb import calculate_matrices, make_list_of_states, calculate_gamma_matrix, calculate_omega_matrix
from atomic_structure_Rb import calculate_boundaries

from sympy.core.numbers import Rational as Integer
