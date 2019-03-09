# -*- coding: utf-8 -*-
# Oscar Gerardo Lazo Arjona
"""FAST is an acronym for FAST Atomic Spectroscopy from Theory.

It calculates the density matrix dynamics and steady states of generic atoms,
rubidium, and cesium using an arbitrary number of states and radiation fields.

>>> from fast import __version__
>>> print __version__
1.1.2

>>> from fast import all_atoms
>>> all_atoms
[85Rb, 87Rb, 133Cs]

"""
# This allows plots to be made remotely via ssh.
from matplotlib import use; use('Agg')

from electric_field import PlaneWave, MotField
from electric_field import electric_field_amplitude_gaussian
from electric_field import electric_field_amplitude_top
from electric_field import electric_field_amplitude_intensity
from misc import Mu, IJ, find_phase_transformation
from misc import formatLij, convolve_with_gaussian, read_result, fprint

from graphic import complex_matrix_plot, plot_Lij
from graphic import Arrow3D, bar_chart_mf, draw_atom3d, draw_mot_field_3d
from graphic import draw_plane_wave_3d, draw_lasers_3d
from graphic import draw_state, excitation, decay, draw_multiplet
from graphic import fancy_matrix_plot, fancy_r_plot, plot_populations

from evolution import write_evolution, run_evolution
from stationary import write_stationary, run_stationary
from misc import compile_code

from atomic_structure import Atom, State, Transition
from atomic_structure import split_fine_to_hyperfine, split_fine_to_magnetic
from atomic_structure import split_hyperfine_to_magnetic
from atomic_structure import calculate_matrices, make_list_of_states
from atomic_structure import calculate_gamma_matrix, calculate_omega_matrix
from atomic_structure import calculate_r_matrices, calculate_boundaries
from atomic_structure import calculate_reduced_matrix_elements
from atomic_structure import vapour_pressure, vapour_number_density
from atomic_structure import vapour_density
from atomic_structure import speed_likely, speed_average, speed_rms
from atomic_structure import collision_rate

from sympy.core.numbers import Rational as Integer

from sympy import init_printing, pprint
from sympy import Symbol, Matrix, symbols
from sympy import conjugate, re, im
from sympy import simplify, KroneckerDelta, Function, Derivative, solve

from symbolic import define_density_matrix, define_laser_variables
from symbolic import polarization_vector
from symbolic import cartesian_to_helicity, helicity_to_cartesian
from symbolic import helicity_dot_product, cartesian_dot_product
from symbolic import define_r_components, define_frequencies
from symbolic import delta_greater, delta_lesser
from symbolic import ket, bra, ketbra, lindblad_operator, lindblad_terms
from symbolic import define_psi_coefficients
from symbolic import define_rho_vector, calculate_A_b
from symbolic import vector_element, phase_transformation

from error_propagation import Measurement

from matplotlib import rcParams

# We set matplotlib to use a nice latex font.
rcParams['mathtext.fontset'] = 'cm'
rcParams['mathtext.rm'] = 'serif'

__version__ = "1.1.2"
all_atoms = [Atom("Rb", 85), Atom("Rb", 87), Atom("Cs", 133)]
