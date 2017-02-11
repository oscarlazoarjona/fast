# -*- coding: utf-8 -*-

########################################################################
# Run all doctests for the core code:
from doctest import testmod
import fast

verbose=True; verbose=False
print testmod(fast                   , verbose=verbose)
print testmod(fast.atomic_structure  , verbose=verbose)
print testmod(fast.symbolic          , verbose=verbose)

print testmod(fast.error_propagation , verbose=verbose)
print testmod(fast.evolution         , verbose=verbose)
print testmod(fast.graphic           , verbose=verbose)
print testmod(fast.misc              , verbose=verbose)
print testmod(fast.rk4               , verbose=verbose)
print testmod(fast.stationary        , verbose=verbose)

########################################################################
# Run all codes from notebook-generated examples.
import doctest_08___Three_level_atom_V_symbolic
import doctest_09___Thermal_States
import doctest_10___States_database

print testmod(doctest_08___Three_level_atom_V_symbolic  , verbose=verbose)
print testmod(doctest_09___Thermal_States  , verbose=verbose)
print testmod(doctest_10___States_database , verbose=verbose)
########################################################################
#Toy examples.
from toy.two_levels import suite
from toy.three_levels_ladder import suite

#Real atom examples.
from real.one_photon_780nm import rb87_suite
