# -*- coding: utf-8 -*-

########################################################################
#Run all doctests:
from doctest import testmod
import fast.atomic_structure

res=testmod(fast.atomic_structure,verbose=True)
print res

########################################################################
#Toy examples.
from toy.two_levels import suite
from toy.three_levels_ladder import suite

#Real atom examples.
from real.one_photon_780nm import rb87_suite

