# -*- coding: utf-8 -*-
r"""This script runs all doctests in FAST."""

########################################################################
# Run all doctests for the core code:
from doctest import testmod
import fast

import doctest_01___Two_level_atom
import doctest_02___Three_level_atom_ladder
import doctest_03___Rb87_one_photon
import doctest_04___Vectors_in_the_helicity_basis_and_the_electric_field
import doctest_05___Two_level_atom_symbolic
import doctest_06___Three_level_atom_ladder_symbolic
import doctest_07___Three_level_atom_Lambda_symbolic
import doctest_08___Three_level_atom_V_symbolic
import doctest_09___Thermal_States
import doctest_10___States_database

verbose = True
# verbose = False
print testmod(fast, verbose=verbose)
print testmod(fast.atomic_structure, verbose=verbose)
print testmod(fast.electric_field, verbose=verbose)
print testmod(fast.magnetic_field, verbose=verbose)
print testmod(fast.symbolic, verbose=verbose)
print testmod(fast.error_propagation, verbose=verbose)
print testmod(fast.evolution, verbose=verbose)
print testmod(fast.graphic, verbose=verbose)
print testmod(fast.misc, verbose=verbose)
print testmod(fast.rk4, verbose=verbose)
print testmod(fast.stationary, verbose=verbose)
print testmod(fast.bloch, verbose=verbose)
print testmod(fast.inhomo, verbose=verbose)

########################################################################
# Run all codes from notebook-generated doctests.

print testmod(doctest_01___Two_level_atom, verbose=verbose)
print testmod(doctest_02___Three_level_atom_ladder, verbose=verbose)
print testmod(doctest_03___Rb87_one_photon, verbose=verbose)
s = doctest_04___Vectors_in_the_helicity_basis_and_the_electric_field
print testmod(s, verbose=verbose)
print testmod(doctest_05___Two_level_atom_symbolic, verbose=verbose)
print testmod(doctest_06___Three_level_atom_ladder_symbolic, verbose=verbose)
print testmod(doctest_07___Three_level_atom_Lambda_symbolic, verbose=verbose)
print testmod(doctest_08___Three_level_atom_V_symbolic, verbose=verbose)
print testmod(doctest_09___Thermal_States, verbose=verbose)
print testmod(doctest_10___States_database, verbose=verbose)
# ########################################################################
# # Toy examples.
# # from toy.two_levels import suite
# # from toy.three_levels_ladder import suite
#
# # Real atom examples.
# # from real.one_photon_780nm import rb87_suite
