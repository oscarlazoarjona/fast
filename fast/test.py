# -*- coding: utf-8 -*-
# Copyright (C) 2017 Oscar Gerardo Lazo Arjona
# mailto: oscar.lazoarjona@physics.ox.ac.uk
r"""This is a template."""
from symbolic import (define_laser_variables, polarization_vector, hamiltonian,
                      define_r_components, helicity_to_cartesian,
                      define_frequencies, vector_element, cartesian_dot_product)
from sympy import symbols, zeros, Matrix, pi, pprint

from sympy import zeros, Matrix, symbols
from fast.bloch import phase_transformation

Ne = 4
Nl = 3

Ep, omega_laser = define_laser_variables(Nl)
detuning_knob = symbols("delta1 delta2", real=True)

epsilonp = [polarization_vector(0, -pi/2, 0, 0, 1) for l in range(Nl)]
xi = [zeros(Ne, Ne) for i in range(Nl)]
coup = [[(2, 0)], [(3, 2)], [(2, 1)]]

xi = [Matrix([[1 if (i, j) in coup[l] or (j, i) in coup else 0
               for j in range(Ne)]
               for i in range(Ne)])
               for l in range(Nl)]

rm = define_r_components(Ne, xi, explicitly_hermitian=True, p=-1, helicity=True)
rm = helicity_to_cartesian(rm)

rp = define_r_components(Ne, xi, explicitly_hermitian=True, p=1, helicity=True)
rp = helicity_to_cartesian(rp)

r = rm + rp
omega_level, omega, gamma = define_frequencies(Ne, True)

RF = phase_transformation(Ne, Nl, rm, xi, return_equations=False)
t = symbols("t", real=True)
RF = [RF[i]*t for i in range(Ne)]

H = hamiltonian(Ep, epsilonp, detuning_knob, rm, omega_level, omega_laser, xi, RWA=True, RF=RF)
H.simplify()
pprint(H)

# i = 0; j = 2; l = 0
# rmij = vector_element(rm, i, j)
# rpij = vector_element(rm, j, i).conjugate()
#
# epsilonpl = epsilonp[l]
# epsilonml = epsilonpl.conjugate()
#
# print i, j, l, rpij, epsilonml
# print 222, rpij[2]
#
# print cartesian_dot_product(epsilonml, rpij)
