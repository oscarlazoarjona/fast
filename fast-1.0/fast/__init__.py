# -*- coding: utf-8 -*-
#Oscar Gerardo Lazo Arjona
"""
FAST is an acronym for FAST Atomic Spectroscopy from Theory. It calculates
the density matrix dynamics and steady states of generic alkali atoms,
rubidium, and cesium using an arbitrary number of states and radiation fields.

"""
__version__="1.0"

#This flag tells whether we are running in Sage or not.
sage_included = 'sage' in globals().keys()
if not sage_included:
	from math import pi,sqrt
else:
	from symbolic import define_density_matrix,define_substitutions


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

from atomic_structure import S,P,D,F,State,Transition
from atomic_structure import split_fine_to_hyperfine, split_fine_to_magnetic, split_hyperfine_to_magnetic
from atomic_structure import calculate_matrices, make_list_of_states, calculate_gamma_matrix, calculate_omega_matrix
from atomic_structure import calculate_r_matrices, calculate_boundaries

from sympy.core.numbers import Rational as Integer

from sympy import init_printing,pprint
from sympy import Symbol,Matrix,symbols
from sympy import I,conjugate,re,im
from sympy import simplify, KroneckerDelta, Function, Derivative, solve

from symbolic import define_density_matrix, define_laser_variables, polarization_vector
from symbolic import cartesian_to_helicity, helicity_to_cartesian, helicity_dot_product
from symbolic import define_r_components, define_frequencies
from symbolic import delta_greater, delta_lesser
from symbolic import ket,bra,lindblad_operator,lindblad_terms
from symbolic import define_psi_coefficients
from symbolic import define_rho_vector,calculate_A_b
from symbolic import vector_element
