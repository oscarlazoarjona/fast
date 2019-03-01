# -*- coding: utf-8 -*-
# Oscar Gerardo Lazo Arjona
"""FAST is an acronym for FAST Atomic Spectroscopy from Theory.

It calculates the density matrix dynamics and steady states of generic atoms,
rubidium, and cesium using an arbitrary number of states and radiation fields.

>>> from fast import __version__
>>> print __version__
1.3

>>> from fast import all_atoms
>>> all_atoms
[85Rb, 87Rb, 133Cs]

"""

from electric_field import (PlaneWave, MotField,
                            electric_field_amplitude_gaussian,
                            electric_field_amplitude_top,
                            electric_field_amplitude_intensity)

from magnetic_field import (lande_g_factors, zeeman_energies,
                            paschen_back_energies)

from angular_momentum import (coupling_matrix_2j, coupling_matrix_3j,
                              angular_momentum_matrix,
                              orbital_spin_nuclear_matrices,
                              spherical_tensor, wigner_d_small,
                              wigner_d, density_matrix_rotation,
                              perm_j, perm_m)

from misc import (Mu, IJ, find_phase_transformation,
                  formatLij, convolve_with_gaussian, read_result, fprint,
                  compile_code)

from graphic import (complex_matrix_plot, plot_Lij,
                     Arrow3D, bar_chart_mf, draw_atom3d, draw_mot_field_3d,
                     draw_plane_wave_3d, draw_lasers_3d,
                     draw_state, excitation, decay, draw_multiplet,
                     fancy_matrix_plot, fancy_r_plot, plot_populations)

from evolution import write_evolution, run_evolution
from stationary import write_stationary, run_stationary

from atomic_structure import (Atom, State, Transition,
                              split_fine_to_hyperfine, split_fine_to_magnetic,
                              split_hyperfine_to_magnetic,
                              calculate_matrices, make_list_of_states,
                              unperturbed_hamiltonian,
                              calculate_gamma_matrix, calculate_omega_matrix,
                              calculate_r_matrices, calculate_boundaries,
                              calculate_reduced_matrix_elements,
                              vapour_pressure, vapour_number_density,
                              vapour_density,
                              speed_likely, speed_average, speed_rms,
                              collision_rate, matrix_element,
                              thermal_state)

from sympy.core.numbers import Rational as Integer

from sympy import (init_printing, pprint, Symbol, Matrix, symbols,
                   conjugate, re, im, simplify, KroneckerDelta,
                   Function, Derivative, solve)

from symbolic import (define_density_matrix, define_laser_variables,
                      polarization_vector, cartesian_to_helicity,
                      helicity_to_cartesian, helicity_dot_product,
                      cartesian_dot_product, define_r_components,
                      define_frequencies, delta_greater, delta_lesser,
                      ket, bra, ketbra, lindblad_operator, lindblad_terms,
                      define_psi_coefficients, define_rho_vector,
                      calculate_A_b, vector_element, phase_transformation)

from error_propagation import Measurement

from bloch import (Unfolding, fast_hamiltonian, fast_rabi_terms,
                   fast_detuning_terms, fast_lindblad_terms,
                   fast_hamiltonian_terms, fast_bloch_equations,
                   fast_steady_state, fast_time_evolution,
                   fast_sweep_steady_state, fast_sweep_time_evolution,
                   observable, electric_succeptibility, radiated_intensity)

from matplotlib import rcParams

# We set matplotlib to use a nice latex font.
rcParams['mathtext.fontset'] = 'cm'
rcParams['mathtext.rm'] = 'serif'

__version__ = "1.3"
all_atoms = [Atom("Rb", 85), Atom("Rb", 87), Atom("Cs", 133)]
