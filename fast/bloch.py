# -*- coding: utf-8 -*-
# ***********************************************************************
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
# ***********************************************************************

# TODO:
#      Make doctests of fast_rabi_terms that show the new features.
r"""This module contains all the fast routines for optical Bloch equations.

Here is an example with rubidum 87.

>>> import fast
>>> element = "Rb"; isotope = 87; N = 5
>>> fine_states = [fast.State(element, isotope, N, 0, 1/fast.Integer(2)),
...                fast.State(element, isotope, N, 1, 3/fast.Integer(2))]

>>> magnetic_states = fast.make_list_of_states(fine_states, "magnetic")
>>> Ne = len(magnetic_states)
>>> Nl = 2
>>> E0 = [1e2, 1e2]
>>> epsilonp = [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]]

>>> omega, gamma, r = fast.calculate_matrices(magnetic_states)
>>> omega_level = [omega[i][0] for i in range(Ne)]
>>> r = fast.helicity_to_cartesian(r, numeric=True)
>>> rm = np.array([[[r[p, i, j]*fast.symbolic.delta_greater(i, j)*a0
...                for j in range(Ne)] for i in range(Ne)]
...                for p in range(3)])

>>> def coupled2(l, i, j):
...     if r[0, i, j] != 0 or \
...        r[1, i, j] != 0 or \
...        r[2, i, j] != 0:
...         if i < j:
...             i, j = j, i
...         if magnetic_states[j].f == 1 and l == 0:
...             return 1.0
...         if magnetic_states[j].f == 2 and l == 1:
...             return 1.0
...         else:
...             return 0.0
...     else:
...         return 0.0

>>> xi = np.array([[[coupled2(l, i, j)
...                for j in range(Ne)] for i in range(Ne)]
...                for l in range(Nl)])

>>> phase = phase_transformation(Ne, Nl, rm, xi, return_equations=False)
>>> detuning_knob = [fast.symbols("delta"+str(i+1)) for i in range(Nl)]
>>> hamiltonian = fast_hamiltonian(E0, epsilonp, detuning_knob, rm,
...                                omega_level, xi, phase)

>>> detuning_knob = [0.0, 0.0]
>>> print hamiltonian(detuning_knob)/hbar_num/2/np.pi*1e-6
[[  0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
   -1.23434184e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j  -9.56117080e-01+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j]
 [  0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   1.10402891e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j  -1.10402891e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j]
 [  0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   1.23434184e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
   -9.56117080e-01+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j]
 [  0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    7.22218019e+01+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
   -1.10402891e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
   -1.10402891e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j]
 [  0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   7.22218019e+01+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    4.27588557e-01+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j  -5.52014453e-01+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j  -1.39649838e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j]
 [  0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   7.22218019e+01+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   4.93736737e-01+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j  -1.48121021e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j]
 [  0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    7.22218019e+01+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   4.27588557e-01+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    5.52014453e-01+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
   -1.39649838e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j]
 [  0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   7.22218019e+01+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   1.10402891e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j  -1.10402891e+00+0.j   0.00000000e+00+0.j]
 [  0.00000000e+00+0.j   1.10402891e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j]
 [ -1.23434184e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   4.27588557e-01+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    7.22218019e+01+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j]
 [  0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   4.93736737e-01+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   7.22218019e+01+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j]
 [  0.00000000e+00+0.j   0.00000000e+00+0.j   1.23434184e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    4.27588557e-01+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   7.22218019e+01+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j]
 [  0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
   -1.10402891e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    2.29161955e+02+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j]
 [ -9.56117080e-01+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j  -5.52014453e-01+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   2.29161955e+02+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j]
 [  0.00000000e+00+0.j  -1.10402891e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   2.29161955e+02+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j]
 [  0.00000000e+00+0.j   0.00000000e+00+0.j  -9.56117080e-01+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    5.52014453e-01+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    2.29161955e+02+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j]
 [  0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   1.10402891e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   2.29161955e+02+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j]
 [  0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   4.95813555e+02+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j]
 [  0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
   -1.10402891e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    4.95813555e+02+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j]
 [  0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j  -1.39649838e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   4.95813555e+02+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j]
 [  0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j  -1.48121021e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   4.95813555e+02+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j]
 [  0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
   -1.39649838e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    4.95813555e+02+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j]
 [  0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j  -1.10402891e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   4.95813555e+02+0.j   0.00000000e+00+0.j]
 [  0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   0.00000000e+00+0.j
    0.00000000e+00+0.j   0.00000000e+00+0.j   4.95813555e+02+0.j]]

"""

from sympy import Symbol, diff, IndexedBase, re, im, symbols
from fast.symbolic import cartesian_dot_product, define_frequencies
from fast.symbolic import define_laser_variables, define_density_matrix

from symbolic import Vector3D, dot, delta_greater, delta_lesser
from misc import part, symbolic_part

import numpy as np
import sympy
from scipy.constants import physical_constants

hbar_num = physical_constants["Planck constant over 2 pi"][0]
e_num = physical_constants["elementary charge"][0]
a0 = physical_constants["Bohr radius"][0]
epsilon_0_num = physical_constants["electric constant"][0]
alpha_num = physical_constants["fine-structure constant"][0]
c_num = physical_constants["speed of light in vacuum"][0]
k_B_num = physical_constants["Boltzmann constant"][0]


def phase_transformation(Ne, Nl, rm, xi, return_equations=False):
    """Returns a phase transformation theta_i.

    The phase transformation is defined in a way such that
        theta1 + omega_level1 = 0.

    >>> xi = np.zeros((1, 2, 2))
    >>> xi[0, 1, 0] = 1.0
    >>> xi[0, 0, 1] = 1.0
    >>> rm = np.zeros((3, 2, 2))
    >>> rm[0, 1, 0] = 1.0
    >>> rm[1, 1, 0] = 1.0
    >>> rm[2, 1, 0] = 1.0
    >>> phase_transformation(2, 1, rm, xi)
    [-omega_1, -omega_1 - varpi_1]

    """
    # We first define the needed variables
    E0, omega_laser = define_laser_variables(Nl)
    theta = [Symbol('theta'+str(i+1)) for i in range(Ne)]

    # We check for the case of xi being a list of matrices.
    if type(xi) == list:
        xi = np.array([[[xi[l][i, j]
                       for j in range(Ne)] for i in range(Ne)]
                       for l in range(Nl)])

    # We find all the equations that the specified problem has to fulfil.
    eqs = []
    for i in range(Ne):
        for j in range(0, i):
            if (rm[0][i, j] != 0) or \
               (rm[1][i, j] != 0) or \
               (rm[2][i, j] != 0):
                for l in range(Nl):
                    if xi[l, i, j] == 1:
                        eqs += [-omega_laser[l] + theta[j] - theta[i]]

    if return_equations:
        return eqs

    # We solve the system of equations.
    # print eqs
    # print theta
    sol = sympy.solve(eqs, theta, dict=True)
    # print sol
    sol = sol[0]
    # We add any missing theta that may be left outside if the system is
    # under determined.

    extra_thetas = []
    for i in range(Ne):
        if theta[i] not in sol.keys():
            sol.update({theta[i]: theta[i]})
            extra_thetas += [theta[i]]

    # We make the solution such that theta1 + omega_level1 = 0.
    omega_level, omega, gamma = define_frequencies(Ne)
    eq_crit = sol[theta[0]] + omega_level[0]
    ss = sympy.solve(eq_crit, extra_thetas[0])[0]
    ss = {extra_thetas[0]: ss}

    sol_simple = [sol[theta[i]].subs(ss)
                  for i in range(Ne)]

    # sol = []
    # for i in range(Ne):
    #     soli = []
    #     for l in range(Nl):
    #         soli += [sympy.diff(sol_simple[theta[i]], omega_laser[l])]
    #     sol += [soli]

    return sol_simple


def define_simplification(omega_level, xi, Nl):
    """Return a simplifying function, its inverse, and simplified frequencies.

    This implements an index iu that labels energies in a non-degenerate
    way.

    >>> Ne = 6
    >>> Nl = 2
    >>> omega_level = [0.0, 100.0, 100.0, 200.0, 200.0, 300.0]
    >>> xi = np.zeros((Nl, Ne, Ne))
    >>> coup = [[(1, 0), (2, 0)], [(3, 0), (4, 0), (5, 0)]]
    >>> for l in range(Nl):
    ...     for pair in coup[l]:
    ...         xi[l, pair[0], pair[1]] = 1.0
    ...         xi[l, pair[1], pair[0]] = 1.0

    >>> aux = define_simplification(omega_level, xi, Nl)
    >>> u, invu, omega_levelu, Neu, xiu = aux
    >>> print omega_levelu
    [0.0, 100.0, 200.0, 300.0]
    >>> print Neu
    4
    >>> print xiu
    [[[ 0.  1.  0.  0.]
      [ 1.  0.  0.  0.]
      [ 0.  0.  0.  0.]
      [ 0.  0.  0.  0.]]
    <BLANKLINE>
     [[ 0.  0.  1.  1.]
      [ 0.  0.  0.  0.]
      [ 1.  0.  0.  0.]
      [ 1.  0.  0.  0.]]]

    """
    try:
        Ne = len(omega_level)
    except:
        Ne = omega_level.shape[0]
    #####################################
    # 1 We calculate the symplifying functions.
    om = omega_level[0]
    iu = 0; Neu = 1
    omega_levelu = [om]
    d = {}; di = {0: 0}
    for i in range(Ne):
        if omega_level[i] != om:
            iu += 1
            om = omega_level[i]
            Neu += 1
            omega_levelu += [om]
            di.update({iu: i})
        d.update({i: iu})

    def u(i):
        return d[i]

    def invu(iu):
        return di[iu]

    #####################################
    # 2 We build the simplified xi.
    Neu = len(omega_levelu)
    xiu = np.array([[[xi[l, invu(i), invu(j)]
                    for j in range(Neu)] for i in range(Neu)]
                    for l in range(Nl)])
    #####################################
    return u, invu, omega_levelu, Neu, xiu


def find_omega_min(omega_levelu, Neu, Nl, xiu):
    r"""Find the smallest transition frequency for each field.

    >>> Ne = 6
    >>> Nl = 2
    >>> omega_level = [0.0, 100.0, 100.0, 200.0, 200.0, 300.0]
    >>> xi = np.zeros((Nl, Ne, Ne))
    >>> coup = [[(1, 0), (2, 0)], [(3, 0), (4, 0), (5, 0)]]
    >>> for l in range(Nl):
    ...     for pair in coup[l]:
    ...         xi[l, pair[0], pair[1]] = 1.0
    ...         xi[l, pair[1], pair[0]] = 1.0

    >>> aux = define_simplification(omega_level, xi, Nl)
    >>> u, invu, omega_levelu, Neu, xiu = aux
    >>> find_omega_min(omega_levelu, Neu, Nl, xiu)
    ([100.0, 200.0], [1, 2], [0, 0])

    """
    omega_min = []; iu0 = []; ju0 = []
    for l in range(Nl):
        omegasl = []
        for iu in range(Neu):
            for ju in range(iu):
                if xiu[l, iu, ju] == 1:
                    omegasl += [(omega_levelu[iu]-omega_levelu[ju], iu, ju)]
        omegasl = list(sorted(omegasl))
        omega_min += [omegasl[0][0]]
        iu0 += [omegasl[0][1]]
        ju0 += [omegasl[0][2]]

    return omega_min, iu0, ju0


def detunings_indices(Neu, Nl, xiu):
    r"""Get the indices of the transitions of all fields.

    They are returned in the form
    [[(i1, j1), (i2, j2)], ...,[(i1, j1)]].
    that is, one list of pairs of indices for each field.

    >>> Ne = 6
    >>> Nl = 2
    >>> omega_level = [0.0, 100.0, 100.0, 200.0, 200.0, 300.0]
    >>> xi = np.zeros((Nl, Ne, Ne))
    >>> coup = [[(1, 0), (2, 0)], [(3, 0), (4, 0), (5, 0)]]
    >>> for l in range(Nl):
    ...     for pair in coup[l]:
    ...         xi[l, pair[0], pair[1]] = 1.0
    ...         xi[l, pair[1], pair[0]] = 1.0

    >>> aux = define_simplification(omega_level, xi, Nl)
    >>> u, invu, omega_levelu, Neu, xiu = aux
    >>> detunings_indices(Neu, Nl, xiu)
    [[(1, 0)], [(2, 0), (3, 0)]]

    """
    pairs = []
    for l in range(Nl):
        ind = []
        for iu in range(Neu):
            for ju in range(iu):
                if xiu[l, iu, ju] == 1:
                    ind += [(iu, ju)]
        pairs += [ind]
    return pairs


def detunings_code(Neu, Nl, pairs, omega_levelu, iu0, ju0):
    r"""Get the code to calculate the simplified detunings.

    >>> Ne = 6
    >>> Nl = 2
    >>> omega_level = [0.0, 100.0, 100.0, 200.0, 200.0, 300.0]
    >>> xi = np.zeros((Nl, Ne, Ne))
    >>> coup = [[(1, 0), (2, 0)], [(3, 0), (4, 0), (5, 0)]]
    >>> for l in range(Nl):
    ...     for pair in coup[l]:
    ...         xi[l, pair[0], pair[1]] = 1.0
    ...         xi[l, pair[1], pair[0]] = 1.0

    >>> aux = define_simplification(omega_level, xi, Nl)
    >>> u, invu, omega_levelu, Neu, xiu = aux
    >>> omega_min, iu0, ju0 = find_omega_min(omega_levelu, Neu, Nl, xiu)
    >>> pairs = detunings_indices(Neu, Nl, xiu)
    >>> print detunings_code(Neu, Nl, pairs, omega_levelu, iu0, ju0)
        delta1_2_1 = detuning_knob[0]
        delta2_3_1 = detuning_knob[1]
        delta2_4_1 = detuning_knob[1] + (-100.0)
    <BLANKLINE>

    """
    code_det = ""
    for l in range(Nl):
        for pair in pairs[l]:
            iu, ju = pair
            code_det += "    delta"+str(l+1)
            code_det += "_"+str(iu+1)
            code_det += "_"+str(ju+1)
            code_det += " = detuning_knob["+str(l)+"]"
            corr = -omega_levelu[iu]+omega_levelu[iu0[l]]
            corr = -omega_levelu[ju0[l]]+omega_levelu[ju] + corr

            if corr != 0:
                code_det += " + ("+str(corr)+")"
            code_det += "\n"
    return code_det


def detunings_combinations(pairs):
    r"""Return all combinations of detunings.

    >>> Ne = 6
    >>> Nl = 2
    >>> omega_level = [0.0, 100.0, 100.0, 200.0, 200.0, 300.0]
    >>> xi = np.zeros((Nl, Ne, Ne))
    >>> coup = [[(1, 0), (2, 0)], [(3, 0), (4, 0), (5, 0)]]
    >>> for l in range(Nl):
    ...     for pair in coup[l]:
    ...         xi[l, pair[0], pair[1]] = 1.0
    ...         xi[l, pair[1], pair[0]] = 1.0

    >>> aux = define_simplification(omega_level, xi, Nl)
    >>> u, invu, omega_levelu, Neu, xiu = aux
    >>> pairs = detunings_indices(Neu, Nl, xiu)
    >>> detunings_combinations(pairs)
    [[(1, 0), (2, 0)], [(1, 0), (3, 0)]]

    """
    def iter(pairs, combs, l):
        combs_n = []
        for i in range(len(combs)):
            for j in range(len(pairs[l])):
                combs_n += [combs[i] + [pairs[l][j]]]
        return combs_n

    Nl = len(pairs)
    combs = [[pairs[0][k]] for k in range(len(pairs[0]))]
    for l in range(1, Nl):
        combs = iter(pairs, combs, 1)

    return combs


def detunings_rewrite(expr, combs, omega_laser, symb_omega_levelu,
                      omega_levelu, iu0, ju0):
    r"""Rewrite a symbolic expression in terms of allowed transition detunings.

    >>> Ne = 6
    >>> Nl = 2
    >>> omega_level = [0.0, 100.0, 100.0, 200.0, 200.0, 300.0]
    >>> xi = np.zeros((Nl, Ne, Ne))
    >>> coup = [[(1, 0), (2, 0)], [(3, 0), (4, 0), (5, 0)]]
    >>> for l in range(Nl):
    ...     for pair in coup[l]:
    ...         xi[l, pair[0], pair[1]] = 1.0
    ...         xi[l, pair[1], pair[0]] = 1.0

    >>> aux = define_simplification(omega_level, xi, Nl)
    >>> u, invu, omega_levelu, Neu, xiu = aux
    >>> omega_min, iu0, ju0 = find_omega_min(omega_levelu, Neu, Nl, xiu)
    >>> pairs = detunings_indices(Neu, Nl, xiu)
    >>> combs = detunings_combinations(pairs)
    >>> symb_omega_levelu, omega, gamma = define_frequencies(Neu)
    >>> E0, omega_laser = define_laser_variables(Nl)

    Most times it is possible to express these combinations of optical
    frequencies in terms of allowed transition detunings.

    >>> expr = +(omega_laser[0]-(symb_omega_levelu[1]-symb_omega_levelu[0]))
    >>> expr += -(omega_laser[1]-(symb_omega_levelu[3]-symb_omega_levelu[0]))
    >>> expr
    -omega_2 + omega_4 + varpi_1 - varpi_2
    >>> detunings_rewrite(expr, combs, omega_laser, symb_omega_levelu,
    ...                   omega_levelu, iu0, ju0)
    '+delta1_2_1-delta2_4_1'

    But some times it is not possible:

    >>> expr = +(omega_laser[1]-(symb_omega_levelu[1]-symb_omega_levelu[0]))
    >>> expr += -(omega_laser[0]-(symb_omega_levelu[3]-symb_omega_levelu[0]))
    >>> expr
    -omega_2 + omega_4 - varpi_1 + varpi_2
    >>> detunings_rewrite(expr, combs, omega_laser, symb_omega_levelu,
    ...                   omega_levelu, iu0, ju0)
    '300.000000000000-detuning_knob[0]+detuning_knob[1]'

    """
    Nl = len(omega_laser)
    Neu = len(symb_omega_levelu)
    # We find the coefficients a_i of the field frequencies.
    a = [diff(expr, omega_laser[l]) for l in range(Nl)]

    # We look for a combination of the detunings obtained with the
    # function detunings_code. For each combination we sum the
    # detunings weighed by a_i.
    success = False
    for comb in combs:
        expr_try = 0
        for l in range(Nl):
            expr_try += a[l]*(omega_laser[l] -
                              symb_omega_levelu[comb[l][0]] +
                              symb_omega_levelu[comb[l][1]])
        if expr-expr_try == 0:
            success = True
            break

    assign = ""
    if success:
        for l in range(Nl):
            if a[l] != 0:
                if a[l] == 1:
                    assign += "+"
                elif a[l] == -1:
                    assign += "-"
                assign += "delta"+str(l+1)
                assign += "_"+str(comb[l][0]+1)
                assign += "_"+str(comb[l][1]+1)
    else:
        # We get the code for Hii using detuning knobs.
        # We find out the remainder terms.
        _remainder = expr - sum([a[l]*omega_laser[l]
                                 for l in range(Nl)])
        # We find the coefficients of the remainder.
        b = [diff(_remainder, symb_omega_levelu[j]) for j in range(Neu)]
        # We calculate the remainder numerically.
        remainder = sum([b[j]*omega_levelu[j] for j in range(Neu)])
        # We add the contributions from the detuning knobs.
        remainder += sum([a[l]*(omega_levelu[iu0[l]] -
                                omega_levelu[ju0[l]])
                          for l in range(Nl)])
        assign = str(remainder)
        # We get the code for Hii using detuning knobs.
        for l in range(Nl):
            if a[l] != 0:
                if a[l] == 1:
                    assign += "+"
                elif a[l] == -1:
                    assign += "-"
                assign += "detuning_knob["+str(l)+"]"

    return assign


def fast_hamiltonian(Ep, epsilonp, detuning_knob, rm, omega_level, xi, theta,
                     file_name=None):
    r"""Return a fast function that returns a Hamiltonian as an array.

    INPUT:

    -  ``Ep`` - A list with the electric field amplitudes (real or complex).
    -  ``epsilonp`` - A list of the polarization vectors of the fields \
    (real or complex).
    -  ``detuning_knob`` - A list of the detunings of each field (relative \
    to the transition of lowest energy).
    -  ``rm`` -     The below-diagonal components
        of the position operator in the cartesian basis:

        .. math::
            \vec{r}^{(-)}_{i j} = [ x_{ij}, y_{ij}, z_{ij} ]
            \hspace{1cm} \forall \hspace{1cm} 0 < j < i
    -  ``omega_level`` - The angular frequencies of each state.
    -  ``xi`` - An array whose ``xi[l, i, j]`` element is 1 if the \
     transition :math:`|i\rangle \rightarrow |j\rangle`\ is driven by field \
     ``l`` and 0 otherwise.
    -  ``theta`` - A list of symbolic expressions representing a phase \
    transformation.
    -  ``file_name`` - A string indicating a file to save the function's \
    code.

    If the arguments Ep, epsilonp, and detuning_knob are symbolic amounts, \
    the returned function will accept numeric values of Ep, epsilonp, and \
    detuning_knob as arguments.

    All quantities should be in SI units.

    EXAMPLES:

    We build an example using states coupled like this:

     --- |4>        --- |5>          --- |6>
      ^              ^                ^
      |              |                |
      |    --- |2>   |     --- |3>    |
    2 |     ^      2 |      ^         | 2
      |   1 |        |    1 |         |
      |     |        |      |         |
    ------------------------------------- |1>

    With the numbers on kets labeling states and the plain numbers labeling
    fields.

    The number of states and fields:
    >>> Ne = 6
    >>> Nl = 2

    We invent some energy levels:
    >>> omega_level = np.array([0.0, 100.0, 100.0, 200.0, 200.0, 300.0])
    >>> omega_level = omega_level*1e6*2*np.pi

    We build the symbol xi, that chooses which laser couples which
    transition.
    >>> xi = np.zeros((Nl, Ne, Ne))
    >>> coup = [[(1, 0), (2, 0)], [(3, 0), (4, 0), (5, 0)]]
    >>> for l in range(Nl):
    ...     for pair in coup[l]:
    ...         xi[l, pair[0], pair[1]] = 1.0
    ...         xi[l, pair[1], pair[0]] = 1.0

    We invent some electric dipole matrix elements:
    >>> from scipy.constants import physical_constants
    >>> a0 = physical_constants["Bohr radius"][0]
    >>> rm = np.zeros((3, Ne, Ne))
    >>> for l in range(Nl):
    ...     for i in range(Ne):
    ...         for j in range(i):
    ...             if xi[l, i, j] != 0:
    ...                 rm[2, i, j] = float(i)*a0

    The phase transformation:
    >>> theta = phase_transformation(Ne, Nl, rm, xi)

    We define the possible arguments:
    >>> from sympy import symbols, pi
    >>> from fast.symbolic import polarization_vector
    >>> detuning_knob = symbols("delta1 delta2")
    >>> detuning_knob_vals = np.array([-1.0, 3.0])*1e6*2*np.pi

    >>> Ep, omega_laser = define_laser_variables(Nl)
    >>> Ep_vals = [1e2, 1e2]

    >>> alpha = symbols("alpha")
    >>> epsilon = polarization_vector(0, pi/2, alpha, 0, 1)
    >>> epsilonp = [epsilon, epsilon]
    >>> epsilonp_vals = [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]]

    There are 8 ways to call fast_hamiltonian:

    1 .- Get a function of detunings, field amplitudes, polarizations:
    >>> H1 = fast_hamiltonian(Ep, epsilonp, detuning_knob, rm,
    ...                       omega_level, xi, theta)

    2 .- Get a function of field amplitudes, polarizations:
    >>> H2 = fast_hamiltonian(Ep, epsilonp, detuning_knob_vals, rm,
    ...                       omega_level, xi, theta)

    3 .- Get a function of detunings, polarizations:
    >>> H3 = fast_hamiltonian(Ep_vals, epsilonp, detuning_knob, rm,
    ...                       omega_level, xi, theta)

    4 .- Get a function of detunings, field amplitudes:
    >>> H4 = fast_hamiltonian(Ep, epsilonp_vals, detuning_knob, rm,
    ...                       omega_level, xi, theta)

    5 .- Get a function of detunings:
    >>> H5 = fast_hamiltonian(Ep_vals, epsilonp_vals, detuning_knob, rm,
    ...                       omega_level, xi, theta)

    6 .- Get a function of field amplitudes:
    >>> H6 = fast_hamiltonian(Ep, epsilonp_vals, detuning_knob_vals, rm,
    ...                       omega_level, xi, theta)

    7 .- Get a function of polarizations:
    >>> H7 = fast_hamiltonian(Ep_vals, epsilonp, detuning_knob_vals, rm,
    ...                       omega_level, xi, theta)

    8 .- Get a function of nothing:
    >>> H8 = fast_hamiltonian(Ep_vals, epsilonp_vals, detuning_knob_vals, rm,
    ...                       omega_level, xi, theta)


    We test all of these combinations.
    >>> print H1(Ep_vals, epsilonp_vals, detuning_knob_vals) \
    ...     /hbar_num/2/np.pi*1e-6
    [[  0.00000000+0.j   0.63977241+0.j   1.27954481+0.j   1.91931722+0.j
        2.55908963+0.j   3.19886203+0.j]
     [  0.63977241+0.j   1.00000000+0.j   0.00000000+0.j   0.00000000+0.j
        0.00000000+0.j   0.00000000+0.j]
     [  1.27954481+0.j   0.00000000+0.j   1.00000000+0.j   0.00000000+0.j
        0.00000000+0.j   0.00000000+0.j]
     [  1.91931722+0.j   0.00000000+0.j   0.00000000+0.j  -3.00000000+0.j
        0.00000000+0.j   0.00000000+0.j]
     [  2.55908963+0.j   0.00000000+0.j   0.00000000+0.j   0.00000000+0.j
       -3.00000000+0.j   0.00000000+0.j]
     [  3.19886203+0.j   0.00000000+0.j   0.00000000+0.j   0.00000000+0.j
        0.00000000+0.j  97.00000000+0.j]]

    >>> print H2(Ep_vals, epsilonp_vals)/hbar_num/2/np.pi*1e-6
    [[  0.00000000+0.j   0.63977241+0.j   1.27954481+0.j   1.91931722+0.j
        2.55908963+0.j   3.19886203+0.j]
     [  0.63977241+0.j   1.00000000+0.j   0.00000000+0.j   0.00000000+0.j
        0.00000000+0.j   0.00000000+0.j]
     [  1.27954481+0.j   0.00000000+0.j   1.00000000+0.j   0.00000000+0.j
        0.00000000+0.j   0.00000000+0.j]
     [  1.91931722+0.j   0.00000000+0.j   0.00000000+0.j  -3.00000000+0.j
        0.00000000+0.j   0.00000000+0.j]
     [  2.55908963+0.j   0.00000000+0.j   0.00000000+0.j   0.00000000+0.j
       -3.00000000+0.j   0.00000000+0.j]
     [  3.19886203+0.j   0.00000000+0.j   0.00000000+0.j   0.00000000+0.j
        0.00000000+0.j  97.00000000+0.j]]

    >>> print H3(epsilonp_vals, detuning_knob_vals)/hbar_num/2/np.pi*1e-6
    [[  0.00000000+0.j   0.63977241+0.j   1.27954481+0.j   1.91931722+0.j
        2.55908963+0.j   3.19886203+0.j]
     [  0.63977241+0.j   1.00000000+0.j   0.00000000+0.j   0.00000000+0.j
        0.00000000+0.j   0.00000000+0.j]
     [  1.27954481+0.j   0.00000000+0.j   1.00000000+0.j   0.00000000+0.j
        0.00000000+0.j   0.00000000+0.j]
     [  1.91931722+0.j   0.00000000+0.j   0.00000000+0.j  -3.00000000+0.j
        0.00000000+0.j   0.00000000+0.j]
     [  2.55908963+0.j   0.00000000+0.j   0.00000000+0.j   0.00000000+0.j
       -3.00000000+0.j   0.00000000+0.j]
     [  3.19886203+0.j   0.00000000+0.j   0.00000000+0.j   0.00000000+0.j
        0.00000000+0.j  97.00000000+0.j]]

    >>> print H4(Ep_vals, detuning_knob_vals)/hbar_num/2/np.pi*1e-6
    [[  0.00000000+0.j   0.63977241+0.j   1.27954481+0.j   1.91931722+0.j
        2.55908963+0.j   3.19886203+0.j]
     [  0.63977241+0.j   1.00000000+0.j   0.00000000+0.j   0.00000000+0.j
        0.00000000+0.j   0.00000000+0.j]
     [  1.27954481+0.j   0.00000000+0.j   1.00000000+0.j   0.00000000+0.j
        0.00000000+0.j   0.00000000+0.j]
     [  1.91931722+0.j   0.00000000+0.j   0.00000000+0.j  -3.00000000+0.j
        0.00000000+0.j   0.00000000+0.j]
     [  2.55908963+0.j   0.00000000+0.j   0.00000000+0.j   0.00000000+0.j
       -3.00000000+0.j   0.00000000+0.j]
     [  3.19886203+0.j   0.00000000+0.j   0.00000000+0.j   0.00000000+0.j
        0.00000000+0.j  97.00000000+0.j]]

    >>> print H5(detuning_knob_vals)/hbar_num/2/np.pi*1e-6
    [[  0.00000000+0.j   0.63977241+0.j   1.27954481+0.j   1.91931722+0.j
        2.55908963+0.j   3.19886203+0.j]
     [  0.63977241+0.j   1.00000000+0.j   0.00000000+0.j   0.00000000+0.j
        0.00000000+0.j   0.00000000+0.j]
     [  1.27954481+0.j   0.00000000+0.j   1.00000000+0.j   0.00000000+0.j
        0.00000000+0.j   0.00000000+0.j]
     [  1.91931722+0.j   0.00000000+0.j   0.00000000+0.j  -3.00000000+0.j
        0.00000000+0.j   0.00000000+0.j]
     [  2.55908963+0.j   0.00000000+0.j   0.00000000+0.j   0.00000000+0.j
       -3.00000000+0.j   0.00000000+0.j]
     [  3.19886203+0.j   0.00000000+0.j   0.00000000+0.j   0.00000000+0.j
        0.00000000+0.j  97.00000000+0.j]]

    >>> print H6(Ep_vals)/hbar_num/2/np.pi*1e-6
    [[  0.00000000+0.j   0.63977241+0.j   1.27954481+0.j   1.91931722+0.j
        2.55908963+0.j   3.19886203+0.j]
     [  0.63977241+0.j   1.00000000+0.j   0.00000000+0.j   0.00000000+0.j
        0.00000000+0.j   0.00000000+0.j]
     [  1.27954481+0.j   0.00000000+0.j   1.00000000+0.j   0.00000000+0.j
        0.00000000+0.j   0.00000000+0.j]
     [  1.91931722+0.j   0.00000000+0.j   0.00000000+0.j  -3.00000000+0.j
        0.00000000+0.j   0.00000000+0.j]
     [  2.55908963+0.j   0.00000000+0.j   0.00000000+0.j   0.00000000+0.j
       -3.00000000+0.j   0.00000000+0.j]
     [  3.19886203+0.j   0.00000000+0.j   0.00000000+0.j   0.00000000+0.j
        0.00000000+0.j  97.00000000+0.j]]

    >>> print H7(epsilonp_vals)/hbar_num/2/np.pi*1e-6
    [[  0.00000000+0.j   0.63977241+0.j   1.27954481+0.j   1.91931722+0.j
        2.55908963+0.j   3.19886203+0.j]
     [  0.63977241+0.j   1.00000000+0.j   0.00000000+0.j   0.00000000+0.j
        0.00000000+0.j   0.00000000+0.j]
     [  1.27954481+0.j   0.00000000+0.j   1.00000000+0.j   0.00000000+0.j
        0.00000000+0.j   0.00000000+0.j]
     [  1.91931722+0.j   0.00000000+0.j   0.00000000+0.j  -3.00000000+0.j
        0.00000000+0.j   0.00000000+0.j]
     [  2.55908963+0.j   0.00000000+0.j   0.00000000+0.j   0.00000000+0.j
       -3.00000000+0.j   0.00000000+0.j]
     [  3.19886203+0.j   0.00000000+0.j   0.00000000+0.j   0.00000000+0.j
        0.00000000+0.j  97.00000000+0.j]]

    >>> print H8()/hbar_num/2/np.pi*1e-6
    [[  0.00000000+0.j   0.63977241+0.j   1.27954481+0.j   1.91931722+0.j
        2.55908963+0.j   3.19886203+0.j]
     [  0.63977241+0.j   1.00000000+0.j   0.00000000+0.j   0.00000000+0.j
        0.00000000+0.j   0.00000000+0.j]
     [  1.27954481+0.j   0.00000000+0.j   1.00000000+0.j   0.00000000+0.j
        0.00000000+0.j   0.00000000+0.j]
     [  1.91931722+0.j   0.00000000+0.j   0.00000000+0.j  -3.00000000+0.j
        0.00000000+0.j   0.00000000+0.j]
     [  2.55908963+0.j   0.00000000+0.j   0.00000000+0.j   0.00000000+0.j
       -3.00000000+0.j   0.00000000+0.j]
     [  3.19886203+0.j   0.00000000+0.j   0.00000000+0.j   0.00000000+0.j
        0.00000000+0.j  97.00000000+0.j]]

    """
    # We determine which arguments are constants.
    if True:
        Nl = len(Ep)
        Ne = np.array(rm[0]).shape[0]
        try:
            Ep = np.array([complex(Ep[l]) for l in range(Nl)])
            variable_Ep = False
        except:
            variable_Ep = True

        try:
            epsilonp = [np.array([complex(epsilonp[l][i]) for i in range(3)])
                        for l in range(Nl)]
            variable_epsilonp = False
        except:
            variable_epsilonp = True

        try:
            detuning_knob = np.array([float(detuning_knob[l])
                                      for l in range(Nl)])
            variable_detuning_knob = False
        except:
            variable_detuning_knob = True

        # We convert rm to a numpy array
        rm = np.array([[[complex(rm[k][i, j])
                       for j in range(Ne)] for i in range(Ne)]
                       for k in range(3)])
    # We establish the arguments of the output function.
    if True:
        code = ""
        code += "def hamiltonian("
        if variable_Ep: code += "Ep, "
        if variable_epsilonp: code += "epsilonp, "
        if variable_detuning_knob: code += "detuning_knob, "
        if code[-2:] == ", ":
            code = code[:-2]
        code += "):\n"

        code += '    r"""A fast calculation of the hamiltonian."""\n'
        code += "    H = np.zeros(("+str(Ne)+", "+str(Ne)+"), complex)\n\n"
    # We get the code for the below-diagonal elements
    # (Rabi frequencies).
    if True:
        code += "    # We calculate the below-diagonal elements.\n"
        for i in range(Ne):
            for j in range(i):
                for l in range(Nl):
                    if xi[l, i, j] == 1.0:
                        # We get the below-diagonal terms.
                        code += "    H["+str(i)+", "+str(j)+"] = "
                        # We get the code for Ep.
                        if variable_Ep:
                            code += "0.5*Ep["+str(l)+"]"

                        else:
                            code += str(0.5*Ep[l])
                        # We get the code for epsilonp dot rm
                        rmij = rm[:, i, j]
                        if variable_epsilonp:
                            code += "*cartesian_dot_product("
                            code += "epsilonp["+str(l)+"],"
                            code += str(list(rmij*e_num))+" )"
                        else:
                            dp = cartesian_dot_product(epsilonp[l], rmij)
                            dp = dp*e_num
                            code += "*("+str(dp)+")"

                        code += "\n"
    # We get the code for the above-diagonal elements
    # (Conjugate Rabi frequencies).
    if True:
        code += "\n"
        code += """    # We calculate the above-diagonal elements.\n"""
        code += """    for i in range("""+str(Ne)+"""):\n"""
        code += """        for j in range(i+1, """+str(Ne)+"""):\n"""
        code += """            H[i, j] = H[j, i].conjugate()\n\n"""
    # We get the code for the diagonal elements (detunings).
    if True:
        code += "    # We calculate the diagonal elements.\n"
        # We build the degeneration simplification and is inverse (to avoid
        # large combinatorics).
        aux = define_simplification(omega_level, xi, Nl)
        u, invu, omega_levelu, Neu, xiu = aux
        # For each field we find the smallest transition frequency, and its
        # simplified indices.
        omega_min, iu0, ju0 = find_omega_min(omega_levelu, Neu, Nl, xiu)
        #####################################
        # We get the code to calculate the non degenerate detunings.
        pairs = detunings_indices(Neu, Nl, xiu)
        if not variable_detuning_knob:
            code += "    detuning_knob = np.zeros("+str(Nl)+")\n"
            for l in range(Nl):
                code += "    detuning_knob["+str(l)+"] = " +\
                    str(detuning_knob[l])+"\n"

        code_det = detunings_code(Neu, Nl, pairs, omega_levelu, iu0, ju0)
        code += code_det
        code += "\n"
        #####################################
        # We find the coefficients a_l that multiply omega_laser_l in
        # H_ii = omega_level_iu + theta_iu = \sum_i a_i varpi_i + remainder
        _omega_level, omega, gamma = define_frequencies(Ne)
        _omega_levelu, omega, gamma = define_frequencies(Neu)
        E0, omega_laser = define_laser_variables(Nl)
        # So we build all combinations.
        combs = detunings_combinations(pairs)
        for i in range(Ne):
            _Hii = theta[i] + _omega_levelu[u(i)]
            aux = (_Hii, combs, omega_laser,
                   _omega_levelu, omega_levelu, iu0, ju0)
            assign = detunings_rewrite(*aux)

            if assign != "":
                code += "    H["+str(i)+", "+str(i)+"] = "+assign+"\n"

    code += "\n"
    code += """    for i in range("""+str(Ne)+"""):\n"""
    code += """        H[i, i] = H[i, i]*"""+str(hbar_num)+"\n"

    code += "    return H\n"

    if file_name is not None:
        f = file(file_name, "w")
        f.write(code)
        f.close()

    hamiltonian = code
    # print code
    exec hamiltonian
    return hamiltonian


class Unfolding(object):
    r"""A class defining a matrix unfolding."""

    def __init__(self, Ne, real=False, lower_triangular=True,
                 normalized=False):
        r"""Return functions to map matrix element indices to unfolded \
        indices.

        This class has as atributes


        -  ``Nrho`` - the size of the unfolded density matrix.
        -  ``real`` - whether the elements of the unfolded density matrix \
        are real.
        -  ``lower_triangular`` - whether the unfolding includes only the \
        lower triangular elements.
        -  ``normalized`` - whether the populations are normalized to 1 (the \
        first population is not included).

        and two functions ``Mu``, ``IJ``. If ``real=True``, ``Mu`` is a \
        function that takes arguments

        -  ``s`` 1 for the real part, -1 for the imaginary part.
        -  ``i`` an index spanning the ``Ne`` states.
        -  ``j`` an index spanning the ``Ne`` states.
        -  ``k`` an index spanning the ``Nv`` velocity classes (0 by default).

        and returns an index ``mu`` spanning the density matrix, ``IJ`` is \
        the inverse of ``Mu``.

        If ``real=False``, ``Mu`` is a function that takes arguments

        -  ``i`` an index spanning the ``Ne`` states.
        -  ``j`` an index spanning the ``Ne`` states.
        -  ``k`` an index spanning the ``Nv`` velocity classes (0 by default).

        and returns an index ``mu`` spanning the density matrix, ``IJ`` is \
        the inverse of ``Mu``.

        For a two-level system:
        >>> Ne = 2
        >>> Nv = 1
        >>> unf = Unfolding(Ne, real=True, lower_triangular=True,
        ...                       normalized=True)

        The second population and the real and imaginary parts of the \
        coherence,
        >>> unf.Mu(1, 1, 1), unf.Mu(1, 1, 0), unf.Mu(-1, 1, 0)
        (0, 1, 2)

        and back again
        >>> unf.IJ(0), unf.IJ(1), unf.IJ(2)
        ((1, 1, 1), (1, 1, 0), (-1, 1, 0))


        >>> Ne = 3
        >>> Nv = 3

        >>> def test_unfolding(Ne, real=False,
        ...                    lower_triangular=True, normalized=False):
        ...
        ...     vect = Unfolding(Ne, real, lower_triangular,
        ...                          normalized)
        ...     if normalized:
        ...         j0 = 1
        ...     else:
        ...         j0 = 0
        ...     for j in range(j0, Ne):
        ...         if real:
        ...             muu = vect.Mu(1, j, j)
        ...             ss, ii, jj = vect.IJ(muu)
        ...             print j, j, muu, j-ii, j-jj, 1-ss
        ...         else:
        ...             muu = vect.Mu(0, j, j)
        ...             ss, ii, jj = vect.IJ(muu)
        ...             print j, j, muu, j-ii, j-jj
        ...     for j in range(Ne):
        ...         for i in range(j+1, Ne):
        ...             if real:
        ...                 muu = vect.Mu(1, i, j)
        ...                 ss, ii, jj = vect.IJ(muu)
        ...                 print i, j, muu, i-ii, j-jj, 1-ss
        ...                 muu = vect.Mu(-1, i, j)
        ...                 ss, ii, jj = vect.IJ(muu)
        ...                 print i, j, muu, i-ii, j-jj, -1-ss
        ...             else:
        ...                 muu = vect.Mu(0, i, j)
        ...                 ss, ii, jj = vect.IJ(muu)
        ...                 print i, j, muu, i-ii, j-jj
        ...             if not lower_triangular:
        ...                 if real:
        ...                     muu = vect.Mu(1, j, i)
        ...                     ss, ii, jj = vect.IJ(muu)
        ...                     print j, i, muu, j-ii, i-jj, 1-ss
        ...                     muu = vect.Mu(-1, j, i)
        ...                     ss, ii, jj = vect.IJ(muu)
        ...                     print j, i, muu, j-ii, i-jj, -1-ss
        ...                 else:
        ...                     muu = vect.Mu(0, j, i)
        ...                     ss, ii, jj = vect.IJ(muu)
        ...                     print i, j, muu, j-ii, i-jj
        ...


        >>> test_unfolding(Ne, False, False, False)
        0 0 0 0 0
        1 1 1 0 0
        2 2 2 0 0
        1 0 3 0 0
        1 0 4 0 0
        2 0 5 0 0
        2 0 6 0 0
        2 1 7 0 0
        2 1 8 0 0


        >>> test_unfolding(Ne, False, False, True)
        1 1 0 0 0
        2 2 1 0 0
        1 0 2 0 0
        1 0 3 0 0
        2 0 4 0 0
        2 0 5 0 0
        2 1 6 0 0
        2 1 7 0 0

        >>> test_unfolding(Ne, False, True, False)
        0 0 0 0 0
        1 1 1 0 0
        2 2 2 0 0
        1 0 3 0 0
        2 0 4 0 0
        2 1 5 0 0

        >>> test_unfolding(Ne, False, True, True)
        1 1 0 0 0
        2 2 1 0 0
        1 0 2 0 0
        2 0 3 0 0
        2 1 4 0 0

        >>> test_unfolding(Ne, True, False, False)
        0 0 0 0 0 0
        1 1 1 0 0 0
        2 2 2 0 0 0
        1 0 3 0 0 0
        1 0 4 0 0 0
        0 1 5 0 0 0
        0 1 6 0 0 0
        2 0 7 0 0 0
        2 0 8 0 0 0
        0 2 9 0 0 0
        0 2 10 0 0 0
        2 1 11 0 0 0
        2 1 12 0 0 0
        1 2 13 0 0 0
        1 2 14 0 0 0

        >>> test_unfolding(Ne, True, False, True)
        1 1 0 0 0 0
        2 2 1 0 0 0
        1 0 2 0 0 0
        1 0 3 0 0 0
        0 1 4 0 0 0
        0 1 5 0 0 0
        2 0 6 0 0 0
        2 0 7 0 0 0
        0 2 8 0 0 0
        0 2 9 0 0 0
        2 1 10 0 0 0
        2 1 11 0 0 0
        1 2 12 0 0 0
        1 2 13 0 0 0

        >>> test_unfolding(Ne, True, True, False)
        0 0 0 0 0 0
        1 1 1 0 0 0
        2 2 2 0 0 0
        1 0 3 0 0 0
        1 0 4 0 0 0
        2 0 5 0 0 0
        2 0 6 0 0 0
        2 1 7 0 0 0
        2 1 8 0 0 0

        >>> test_unfolding(Ne, True, True, True)
        1 1 0 0 0 0
        2 2 1 0 0 0
        1 0 2 0 0 0
        1 0 3 0 0 0
        2 0 4 0 0 0
        2 0 5 0 0 0
        2 1 6 0 0 0
        2 1 7 0 0 0

        """
        real_map = {}; real_map_inv = {}
        comp_map = {}; comp_map_inv = {}
        mu_real = 0
        mu_comp = 0

        # We get the mappings of populations.
        if normalized:
            start = 1
        else:
            start = 0
        for i in range(start, Ne):
            comp_map.update({(0, i, i): mu_comp})
            comp_map_inv.update({mu_comp: (0, i, i)})
            mu_comp += 1
            real_map.update({(1, i, i): mu_real})
            real_map_inv.update({mu_real: (1, i, i)})
            mu_real += 1

            real_map.update({(-1, i, i): None})
            real_map_inv.update({None: (-1, i, i)})

        # We get the mappings for coherences.
        for j in range(Ne):
            for i in range(j+1, Ne):
                comp_map.update({(0, i, j): mu_comp})
                comp_map_inv.update({mu_comp: (0, i, j)})
                mu_comp += 1
                for s in [1, -1]:
                    real_map.update({(s, i, j): mu_real})
                    real_map_inv.update({mu_real: (s, i, j)})
                    mu_real += 1

                if not lower_triangular:
                    comp_map.update({(0, j, i): mu_comp})
                    comp_map_inv.update({mu_comp: (0, j, i)})
                    mu_comp += 1
                    for s in [1, -1]:
                        real_map.update({(s, j, i): mu_real})
                        real_map_inv.update({mu_real: (s, j, i)})
                        mu_real += 1

            # if normalized:
            #     three.pop((0, 0))
            #     four.pop((0, 0, 1))
            #
            #     three_inv.pop(0)
            #     four_inv.pop(0)

        Nrho_real = mu_real
        Nrho_comp = mu_comp

        def Mu_comp(s, i, j):
            return comp_map[(s, i, j)]

        def Mu_real(s, i, j):
            return real_map[(s, i, j)]

        def IJ_comp(mu):
            return comp_map_inv[mu]

        def IJ_real(mu):
            return real_map_inv[mu]

        if real:
            self.Mu = Mu_real
            self.IJ = IJ_real
            self.Nrho = Nrho_real
            self.map = real_map
        else:
            self.Mu = Mu_comp
            self.IJ = IJ_comp
            self.Nrho = Nrho_comp
            self.map = comp_map
        self.lower_triangular = lower_triangular
        self.real = real
        self.normalized = normalized
        self.Ne = Ne

    def __call__(self, rho):
        r"""Unfold a matrix into a vector.

        The input of this function can be a numpy array or a sympy Matrix.

        >>> unfolding = Unfolding(2, real=True, lower_triangular=True,
        ...                       normalized=True)
        >>> rhos = np.array([[0.6, 1+2j], [1-2j, 0.4]])
        >>> print unfolding(rhos)
        [ 0.4  1.  -2. ]

        >>> from fast import define_density_matrix
        >>> from sympy import pprint
        >>> rho = define_density_matrix(2)
        >>> pprint(unfolding(rho), use_unicode=False)
        [  rho22  ]
        [         ]
        [re(rho21)]
        [         ]
        [im(rho21)]

        """
        Nrho = self.Nrho
        IJ = self.IJ
        if isinstance(rho, np.ndarray):
            if self.real:
                rhov = np.zeros(Nrho)
            else:
                rhov = np.zeros(Nrho, complex)
            numeric = True
        elif isinstance(rho, sympy.Matrix):
            rhov = sympy.zeros(Nrho, 1)
            numeric = False
        else:
            raise ValueError

        for mu in range(Nrho):
            s, i, j = IJ(mu)
            if numeric:
                rhomu = part(rho[i, j], s)
            else:
                rhomu = symbolic_part(rho[i, j], s)
            rhov[mu] = rhomu

        return rhov

    def inverse(self, rhov, time_derivative=False):
        r"""Fold a vector into a matrix.

        The input of this function can be a numpy array or a sympy Matrix.

        If the input is understood to represent the time derivative of a
        density matrix, then the flag time_derivative must be set to True.

        >>> unfolding = Unfolding(2, real=True, lower_triangular=True,
        ...                       normalized=True)
        >>> rhos = np.array([[0.6, 1+2j], [1-2j, 0.4]])
        >>> print rhos == unfolding.inverse(unfolding(rhos))
        [[ True  True]
         [ True  True]]

        >>> from fast import define_density_matrix
        >>> from sympy import pprint
        >>> rho = define_density_matrix(2)
        >>> pprint(unfolding.inverse(unfolding(rho)), use_unicode=False)
        [      -rho22 + 1         re(rho21) - I*im(rho21)]
        [                                                ]
        [re(rho21) + I*im(rho21)           rho22         ]


        >>> rhops = np.array([[0.0, 0.0],
        ...                   [0.0, 0.0]])

        >>> print unfolding.inverse(unfolding(rhops), True)
        [[-0.-0.j  0.-0.j]
         [ 0.+0.j  0.+0.j]]

        """
        Ne = self.Ne
        Nrho = self.Nrho
        IJ = self.IJ

        if isinstance(rhov, np.ndarray):
            rho = np.zeros((Ne, Ne), complex)
            numeric = True
        elif isinstance(rhov, sympy.Matrix):
            rho = sympy.zeros(Ne, Ne)
            numeric = False

        for mu in range(Nrho):
            s, i, j = IJ(mu)
            if numeric:
                if s == 1:
                    rho[i, j] += rhov[mu]
                elif s == -1:
                    rho[i, j] += 1j*rhov[mu]
                elif s == 0:
                    rho[i, j] += rhov[mu]
            else:
                if s == 1:
                    rho[i, j] += rhov[mu]
                elif s == -1:
                    rho[i, j] += sympy.I*rhov[mu]
                elif s == 0:
                    rho[i, j] += rhov[mu]

        if self.lower_triangular:
            for i in range(Ne):
                for j in range(i):
                    rho[j, i] = rho[i, j].conjugate()

        if self.normalized:
            if time_derivative:
                rho[0, 0] = -sum([rho[i, i] for i in range(1, Ne)])
            else:
                rho[0, 0] = 1-sum([rho[i, i] for i in range(1, Ne)])

        return rho


def independent_get_coefficients(coef, rhouv, s, i, j, k, u, v,
                                 unfolding, matrix_form):
    r"""Get the indices mu, nu, and term coefficients for linear terms.

    >>> from fast.symbolic import define_density_matrix
    >>> Ne = 2
    >>> coef = 1+2j
    >>> rhouv = define_density_matrix(Ne)[1, 1]
    >>> s, i, j, k, u, v = (1, 1, 0, 1, 1, 1)
    >>> unfolding = Unfolding(Ne, real=True, normalized=True)

    >>> independent_get_coefficients(coef, rhouv, s, i, j, k, u, v,
    ...                              unfolding, False)
    [[1, None, -2.00000000000000, False, False]]

    """
    if matrix_form:
        coef = -coef
    Mu = unfolding.Mu
    mu = Mu(s, i, j)
    rhouv_isconjugated = False
    if s == 1:
        coef_list = [[mu, None, -im(coef), matrix_form, rhouv_isconjugated]]
    elif s == -1:
        coef_list = [[mu, None, re(coef), matrix_form, rhouv_isconjugated]]
    else:
        coef_list = [[mu, None, coef, matrix_form, rhouv_isconjugated]]
    return coef_list


def linear_get_coefficients(coef, rhouv, s, i, j, k, u, v,
                            unfolding, matrix_form):
    r"""Get the indices mu, nu, and term coefficients for linear terms.

    We determine mu and nu, the indices labeling the density matrix components
          d rho[mu] /dt = sum_nu A[mu, nu]*rho[nu]
    for this complex and rho_u,v.

    >>> from fast.symbolic import define_density_matrix
    >>> Ne = 2
    >>> coef = 1+2j
    >>> rhouv = define_density_matrix(Ne)[1, 1]
    >>> s, i, j, k, u, v = (1, 1, 0, 1, 1, 1)
    >>> unfolding = Unfolding(Ne, real=True, normalized=True)

    >>> linear_get_coefficients(coef, rhouv, s, i, j, k, u, v,
    ...                              unfolding, False)
    [[1, 0, -2.00000000000000, False, False]]

    """
    Ne = unfolding.Ne
    Mu = unfolding.Mu
    # We determine mu, the index labeling the equation.
    mu = Mu(s, i, j)

    if unfolding.normalized and u == 0 and v == 0:
        # We find the nu and coefficients for a term of the form.
        # coef*rho_{00} = coef*(1-sum_{i=1}^{Ne-1} rho_{ii})
        if unfolding.real:
            ss = 1
        else:
            ss = 0

        mu11 = Mu(ss, 1, 1)
        muNeNe = Mu(ss, Ne-1, Ne-1)
        rhouv_isconjugated = False
        if s == 1:
            coef_list = [[mu, nu, im(coef), matrix_form, rhouv_isconjugated]
                         for nu in range(mu11, muNeNe+1)]
        elif s == -1:
            coef_list = [[mu, nu, -re(coef), matrix_form, rhouv_isconjugated]
                         for nu in range(mu11, muNeNe+1)]
        elif s == 0:
            coef_list = [[mu, nu, -coef, matrix_form, rhouv_isconjugated]
                         for nu in range(mu11, muNeNe+1)]
        return coef_list

    #####################################################################

    if (unfolding.lower_triangular and
       isinstance(rhouv, sympy.conjugate)):
        u, v = (v, u)
        rhouv_isconjugated = True
    else:
        rhouv_isconjugated = False
    # If the unfolding is real, there are two terms for this
    # component rhouv of equation mu.
    if unfolding.real:
        nur = Mu(1, u, v)
        nui = Mu(-1, u, v)
    else:
        nu = Mu(0, u, v)
    #####################################################################
    # We determine the coefficients for each term.
    if unfolding.real:
        # There are two sets of forumas for the coefficients depending
        # on whether rhouv_isconjugated.
        # re(I*x*conjugate(y)) = -im(x)*re(y) + re(x)*im(y)
        # re(I*x*y)            = -im(x)*re(y) - re(x)*im(y)
        # im(I*x*conjugate(y)) = +re(x)*re(y) + im(x)*im(y)
        # im(I*x*y)            = +re(x)*re(y) - im(x)*im(y)
        if s == 1:
            # The real part
            if rhouv_isconjugated:
                coef_rerhouv = -im(coef)
                coef_imrhouv = re(coef)
            else:
                coef_rerhouv = -im(coef)
                coef_imrhouv = -re(coef)
        elif s == -1:
            if rhouv_isconjugated:
                coef_rerhouv = re(coef)
                coef_imrhouv = im(coef)
            else:
                coef_rerhouv = re(coef)
                coef_imrhouv = -im(coef)

        coef_list = [[mu, nur, coef_rerhouv, matrix_form, rhouv_isconjugated]]
        if nui is not None:
            coef_list += [[mu, nui, coef_imrhouv,
                           matrix_form, rhouv_isconjugated]]
    else:
        coef_list = [[mu, nu, coef, matrix_form, rhouv_isconjugated]]

    return coef_list


def term_code(mu, nu, coef, matrix_form, rhouv_isconjugated, linear=True):
    r"""Get code to calculate a linear term.

    >>> term_code(1, 0, 33, False, False, True)
    '    rhs[1] += (33)*rho[0]\n'

    """
    if coef == 0:
        return ""
    coef = str(coef)

    # We change E_{0i} -> E0[i-1]
    ini = coef.find("E_{0")
    fin = coef.find("}")
    if ini != -1:
        l = int(coef[ini+4: fin])
        coef = coef[:ini]+"Ep["+str(l-1)+"]"+coef[fin+1:]

    # We change r[i, j] -> r[:, i, j]
    coef = coef.replace("rp[", "rp[:, ")
    coef = coef.replace("rm[", "rm[:, ")

    # We change symbolic complex-operations into fast numpy functions.
    coef = coef.replace("conjugate(", "np.conjugate(")
    coef = coef.replace("re(", "np.real(")
    coef = coef.replace("im(", "np.imag(")

    coef = coef.replace("*I", "j")

    if not linear:
        if matrix_form:
            s = "    b["+str(mu)+"] += "+coef+"\n"
        else:
            s = "    rhs["+str(mu)+"] += "+coef+"\n"
        return s

    # We add the syntax to calculate the term and store it in memory.
    s = "    "
    if matrix_form:
        s += "A["+str(mu)+", "+str(nu)+"] += "+coef+"\n"
    else:
        s += "rhs["+str(mu)+"] += ("+coef+")"
        if rhouv_isconjugated:
            s += "*np.conjugate(rho["+str(nu)+'])\n'
        else:
            s += "*rho["+str(nu)+']\n'

    return s


def fast_rabi_terms(Ep, epsilonp, rm, xi, theta, unfolding,
                    matrix_form=False, file_name=None, return_code=False):
    r"""Return a fast function that returns the Rabi frequency terms.

    We test a basic two-level system.

    >>> import numpy as np
    >>> from scipy.constants import physical_constants
    >>> from sympy import Matrix
    >>> from fast.electric_field import electric_field_amplitude_top
    >>> from fast.symbolic import define_laser_variables, polarization_vector

    >>> Ne = 2
    >>> Nl = 1
    >>> a0 = physical_constants["Bohr radius"][0]
    >>> rm = [np.array([[0, 0], [a0, 0]]),
    ...       np.array([[0, 0], [0, 0]]),
    ...       np.array([[0, 0], [0, 0]])]
    >>> xi = np.array([[[0, 1], [1, 0]]])
    >>> theta = phase_transformation(Ne, Nl, rm, xi)

    We define symbolic variables to be used as token arguments.
    >>> Eps = [electric_field_amplitude_top(1e-3, 1e-3, 1, "SI")]
    >>> Ep, omega_laser = define_laser_variables(Nl)
    >>> epsilonps = [polarization_vector(0, 0, 0, 0, 1)]

    An map to unfold the density matrix.
    >>> unfolding = Unfolding(Ne, True, True, True)

    We obtain a function to calculate Rabi frequency terms.
    >>> rabi_terms = fast_rabi_terms(Ep, epsilonps, rm, xi, theta, unfolding)

    Apply this to a density matrix.
    >>> rhos = np.array([[0.6, 3+2j],
    ...                  [3-2j, 0.4]])
    >>> rhosv = unfolding(rhos)
    >>> rhs_rabi = rabi_terms(rhosv, Eps)
    >>> print rhs_rabi
    [-55680831.47404964         0.           2784041.57370248]

    """
    if not unfolding.lower_triangular:
        mes = "It is very inefficient to solve using all components of the "
        mes += "density matrix. Better set lower_triangular=True in Unfolding."
        raise NotImplementedError(mes)
    if matrix_form and (not unfolding.real) and (unfolding.lower_triangular):
        mes = "It is not possible to express the equations in matrix form "
        mes += "for complex lower triangular components only."
        raise ValueError(mes)
    Nl = len(Ep)
    Ne = unfolding.Ne
    # We determine which arguments are constants.
    if True:
        try:
            Ep = np.array([complex(Ep[l]) for l in range(Nl)])
            variable_Ep = False
        except:
            variable_Ep = True

        try:
            epsilonp = [np.array([complex(epsilonp[l][i]) for i in range(3)])
                        for l in range(Nl)]
            variable_epsilonp = False
        except:
            variable_epsilonp = True
    # We unpack variables.
    if True:
        Nrho = unfolding.Nrho
        # Mu = unfolding.Mu
        IJ = unfolding.IJ
        normalized = unfolding.normalized
        lower_triangular = unfolding.lower_triangular

        # The conjugate stuff.
        Em = [ii.conjugate() for ii in Ep]
        if variable_epsilonp:
            epsilonm = Vector3D(epsilonp.args[0].conjugate())
        else:
            epsilonm = [ii.conjugate() for ii in epsilonp]
        # We convert rm to a numpy array
        rm = np.array([[[complex(rm[k][i, j])
                       for j in range(Ne)] for i in range(Ne)]
                       for k in range(3)])
        rp = np.array([rm[ii].T.conjugate() for ii in range(3)])
        rm_aux = Vector3D(IndexedBase("rm", shape=(Ne, Ne)))
        rp_aux = Vector3D(IndexedBase("rp", shape=(Ne, Ne)))

        # We define needed matrices.
        rho = define_density_matrix(Ne, explicitly_hermitian=lower_triangular,
                                    normalized=normalized)

        if variable_epsilonp:
            Omega = [[sum([xi[l, i, j] *
                          (Ep[l]*dot(epsilonp[l], rm_aux[i, j]) *
                          delta_greater(i, j) +
                           Em[l]*dot(epsilonm[l], rp_aux[i, j]) *
                          delta_lesser(i, j))
                           for l in range(Nl)])
                     for j in range(Ne)] for i in range(Ne)]
        else:
            Omega = [[sum([xi[l, i, j] *
                          (Ep[l]*dot(epsilonp[l], rm[:, i, j]) +
                           Em[l]*dot(epsilonm[l], rp[:, i, j]))
                           for l in range(Nl)])
                     for j in range(Ne)] for i in range(Ne)]
    # We establish the arguments of the output function.
    if True:
        code = ""
        code += "def rabi_terms("
        if not matrix_form: code += "rho, "
        if variable_Ep: code += "Ep, "
        if variable_epsilonp: code += "epsilonp, "
        if code[-2:] == ", ":
            code = code[:-2]
        code += "):\n"
        code += '    r"""A fast calculation of the Rabi terms."""\n'
    # We initialize the output and auxiliaries.
    if True:
        # We introduce the factor that multiplies all terms.
        if unfolding.real:
            code += "    fact = "+str(e_num/hbar_num)+"\n\n"
        else:
            code += "    fact = "+str(1j*e_num/hbar_num)+"\n\n"
        if variable_epsilonp:
            # We put rm and rp into the code
            np.set_printoptions(threshold=np.nan)
            code += "    rm = np."+rm.__repr__()+"\n\n"
            code += "    rp = np."+rp.__repr__()+"\n\n"
            code += "    def dot(epsilon, rij):\n"
            code += "        return epsilon[0]*rij[0]"
            code += " + epsilon[1]*rij[1]"
            code += " + epsilon[2]*rij[2]\n\n"
            np.set_printoptions(threshold=1000)

        if matrix_form:
            code += "    A = np.zeros(("+str(Nrho)+", "+str(Nrho)
            if not unfolding.real:
                code += "), complex)\n\n"
            else:
                code += "))\n\n"
            if unfolding.normalized:
                code += "    b = np.zeros(("+str(Nrho)
                if not unfolding.real:
                    code += "), complex)\n\n"
                else:
                    code += "))\n\n"
        else:
            code += "    rhs = np.zeros(("+str(Nrho)
            if not unfolding.real:
                code += "), complex)\n\n"
            else:
                code += "))\n\n"

    # We write code for the linear terms.
    rho00_terms = []
    for mu in range(Nrho):
        s, i, j = IJ(mu)
        for l in range(Nl):
            for k in range(Ne):
                if xi[l, i, k] == 1:
                    # There is a I* Omega_l,i,k * rho_k,j term.
                    u = k; v = j
                    args = (0.5*Omega[i][k], rho[k, j], s, i, j, k, u, v,
                            unfolding, matrix_form)
                    term_list = linear_get_coefficients(*args)
                    for term in term_list:
                        code += term_code(*term)

                    # We keep note that there was a term with rho00.
                    if k == 0 and j == 0:
                        rho00_terms += [args]

                if xi[l, k, j] == 1:
                    # There is a -I * Omega_l,k,j * rho_i,k term.
                    u = i; v = k
                    args = (-0.5*Omega[k][j], rho[i, k], s, i, j, k, u, v,
                            unfolding, matrix_form)
                    term_list = linear_get_coefficients(*args)
                    for term in term_list:
                        code += term_code(*term)
                    # We keep note that there was a term with rho00.
                    if i == 0 and k == 0:
                        rho00_terms += [args]

    # We write code for the independent terms.
    if unfolding.normalized:
        code += "\n    # Independent terms:\n"
        for term in rho00_terms:
            coef_list = independent_get_coefficients(*term)
            for coef in coef_list:
                code += term_code(*coef, linear=False)
    # We finish the code.
    if True:
        if matrix_form:
            if unfolding.normalized:
                code += "    A *= fact\n"
                code += "    b *= fact\n"
                code += "    return A, b\n"
            else:
                code += "    A *= fact\n"
                code += "    return A\n"
        else:
            code += "    rhs *= fact\n"
            code += "    return rhs\n"
    # We write the code to file if provided, and execute it.
    if True:
        if file_name is not None:
            f = file(file_name+".py", "w")
            f.write(code)
            f.close()

        rabi_terms = code
        if not return_code:
            exec rabi_terms
    return rabi_terms


def fast_detuning_terms(detuning_knob, omega_level, xi, theta, unfolding,
                        matrix_form=False, file_name=None, return_code=False):
    r"""Return a fast function that returns the detuning terms.

    >>> from sympy import Matrix, symbols
    >>> from scipy.constants import physical_constants
    >>> Ne = 2
    >>> Nl = 1
    >>> unfolding = Unfolding(Ne, True, True, True)

    >>> a0 = physical_constants["Bohr radius"][0]
    >>> rm = [Matrix([[0, 0], [a0, 0]]),
    ...       Matrix([[0, 0], [0, 0]]),
    ...       Matrix([[0, 0], [0, 0]])]
    >>> xi = np.array([[[0, 1], [1, 0]]])
    >>> omega_level = [1, 100]
    >>> rhos = np.array([[0.6, 3+2j],
    ...                  [3-2j, 0.4]])
    >>> rhosv = unfolding(rhos)
    >>> detuning_knobs = [1.0]
    >>> theta = phase_transformation(Ne, Nl, rm, xi)
    >>> detuning_knob = [symbols("delta1")]

    >>> detuning_terms = fast_detuning_terms(detuning_knob, omega_level,
    ...                                      xi, theta,
    ...                                      unfolding)

    >>> print detuning_terms(rhosv, detuning_knobs)
    [ 0.  2.  3.]

    """
    # We unpack variables.
    if True:
        Ne = unfolding.Ne
        Nrho = unfolding.Nrho
        Nl = xi.shape[0]
        IJ = unfolding.IJ
        Mu = unfolding.Mu
    # We determine which arguments are constants.
    if True:
        try:
            detuning_knob = np.array([float(detuning_knob[l])
                                      for l in range(Nl)])
            variable_detuning_knob = False
        except:
            variable_detuning_knob = True
    # We establish the arguments of the output function.
    if True:
        code = ""
        code += "def detuning_terms("
        if not matrix_form: code += "rho, "
        if variable_detuning_knob: code += "detuning_knob, "
        if code[-2:] == ", ":
            code = code[:-2]
        code += "):\n"

        code += '    r"""A fast calculation of the detuning terms."""\n'
    # We initialize the output and auxiliaries.
    if True:
        # We introduce the factor that multiplies all terms.
        if unfolding.real:
            code += "    fact = 1.0\n\n"
        else:
            code += "    fact = 1.0j\n\n"

        if matrix_form:
            code += "    A = np.zeros(("+str(Nrho)+", "+str(Nrho)
            if not unfolding.real:
                code += "), complex)\n\n"
            else:
                code += "))\n\n"
            if unfolding.normalized:
                code += "    b = np.zeros(("+str(Nrho)
                if not unfolding.real:
                    code += "), complex)\n\n"
                else:
                    code += "))\n\n"
        else:
            code += "    rhs = np.zeros(("+str(Nrho)
            if not unfolding.real:
                code += "), complex)\n\n"
            else:
                code += "))\n\n"

    # We build the degeneration simplification and is inverse (to avoid
    # large combinatorics).
    aux = define_simplification(omega_level, xi, Nl)
    u, invu, omega_levelu, Neu, xiu = aux
    # For each field we find the smallest transition frequency, and its
    # simplified indices.
    omega_min, iu0, ju0 = find_omega_min(omega_levelu, Neu, Nl, xiu)
    #####################################
    # We get the code to calculate the non degenerate detunings.
    pairs = detunings_indices(Neu, Nl, xiu)
    if not variable_detuning_knob:
        code += "    detuning_knob = np.zeros("+str(Nl)+")\n"
        for l in range(Nl):
            code += "    detuning_knob["+str(l)+"] = " +\
                str(detuning_knob[l])+"\n"

    code_det = detunings_code(Neu, Nl, pairs, omega_levelu, iu0, ju0)
    code += code_det
    code += "\n"
    #####################################
    # There is a term
    # I * Theta_ij * rho_ij = I * (omega_level_j - omega_level_i
    #                              theta_j - theta_i)
    # for all i != j.
    # This term can be re expressed as
    # re(Theta_ij*rho_ij) = - Theta_ij * im(rho_ij)
    # im(Theta_ij*rho_ij) = + Theta_ij * re(rho_ij)
    _omega_level, omega, gamma = define_frequencies(Ne)
    _omega_levelu, omega, gamma = define_frequencies(Neu)
    E0, omega_laser = define_laser_variables(Nl)
    # We build all combinations.
    combs = detunings_combinations(pairs)
    # We add all terms.
    for mu in range(Nrho):
        s, i, j = IJ(mu)
        if i != j:
            _Thetaij = _omega_levelu[u(j)] - _omega_levelu[u(i)]
            _Thetaij += theta[j] - theta[i]
            aux = (_Thetaij, combs, omega_laser,
                   _omega_levelu, omega_levelu, iu0, ju0)
            assign = detunings_rewrite(*aux)

            if assign != "":
                if s == 0:
                    nu = mu
                elif s == 1:
                    assign = "-(%s)" % assign
                    nu = Mu(-s, i, j)
                elif s == -1:
                    assign = "+(%s)" % assign
                    nu = Mu(-s, i, j)

                if matrix_form:
                    term_code = "    A[%s, %s] = %s\n" % (mu, nu, assign)
                else:
                    term_code = "    rhs[%s] = (%s)*rho[%s]\n" % (mu,
                                                                  assign, nu)
            else:
                term_code = ""
            code += term_code
    #####################################

    # We finish the code.
    if True:
        if matrix_form:
            if unfolding.normalized:
                code += "    A *= fact\n"
                code += "    b *= fact\n"
                code += "    return A, b\n"
            else:
                code += "    A *= fact\n"
                code += "    return A\n"
        else:
            code += "    rhs *= fact\n"
            code += "    return rhs\n"
    # We write the code to file if provided, and execute it.
    if True:
        if file_name is not None:
            f = file(file_name+".py", "w")
            f.write(code)
            f.close()

        detuning_terms = code
        if not return_code:
            exec detuning_terms
    return detuning_terms


def fast_lindblad_terms(gamma, unfolding, matrix_form=False, file_name=None,
                        return_code=False):
    r"""Return a fast function that returns the Lindblad terms.

    We test a basic two-level system.

    >>> import numpy as np
    >>> Ne = 2
    >>> gamma21 = 2*np.pi*6e6
    >>> gamma = np.array([[0.0, -gamma21],
    ...                   [gamma21, 0.0]])
    >>> rhos = np.array([[0.6, 3+2j],
    ...                  [3-2j, 0.4]])

    An map to unfold the density matrix.
    >>> unfolding = Unfolding(Ne, True, True, True)

    We obtain a function to calculate Lindblad terms.
    >>> lindblad_terms = fast_lindblad_terms(gamma, unfolding)

    Apply this to a density matrix.
    >>> rhos = np.array([[0.6, 3+2j],
    ...                  [3-2j, 0.4]])
    >>> rhosv = unfolding(rhos)
    >>> rhs_lindblad = lindblad_terms(rhosv)
    >>> print rhs_lindblad
    [-15079644.73724    -56548667.76450001  37699111.843     ]

    """
    Ne = unfolding.Ne
    Nrho = unfolding.Nrho
    Mu = unfolding.Mu

    # We establish the arguments of the output function.
    if True:
        code = ""
        code += "def lindblad_terms("
        if not matrix_form: code += "rho, "
        if code[-2:] == ", ": code = code[:-2]
        code += "):\n"
    # We initialize the output and auxiliaries.
    if True:
        # We introduce the factor that multiplies all terms.
        if matrix_form:
            code += "    A = np.zeros(("+str(Nrho)+", "+str(Nrho)
            if not unfolding.real:
                code += "), complex)\n\n"
            else:
                code += "))\n\n"
            if unfolding.normalized:
                code += "    b = np.zeros(("+str(Nrho)
                if not unfolding.real:
                    code += "), complex)\n\n"
                else:
                    code += "))\n\n"
        else:
            code += "    rhs = np.zeros(("+str(Nrho)
            if not unfolding.real:
                code += "), complex)\n\n"
            else:
                code += "))\n\n"

    for a in range(Ne):
        for b in range(a):
            # The first term is of the from
            # gamma_ab * rho_aa |b><b|
            if not (unfolding.normalized and b == 0):
                coef = gamma[a, b]
                if unfolding.real:
                    mu = Mu(1, b, b)
                    nu = Mu(1, a, a)
                else:
                    mu = Mu(0, b, b)
                    nu = Mu(0, a, a)
                code += term_code(mu, nu, coef, matrix_form, False)

            # The second term is of the form
            #  sum_j -gamma_ab/2 rho_aj |a><j|
            # for a lower triangular unfolding, this j runs from 1 to a.
            for j in range(a):
                coef = -gamma[a, b]*0.5
                if unfolding.real:
                    mur = Mu(1, a, j)
                    code += term_code(mur, mur, coef, matrix_form, False)
                    mui = Mu(-1, a, j)
                    code += term_code(mui, mui, coef, matrix_form, False)
                else:
                    mu = Mu(0, a, j)
                    code += term_code(mu, mu, coef, matrix_form, False)

            # The third term is of the form
            #  - sum_i 1/2 rho_ia |i><a|
            # for a lower triangular unfolding, this i runs from a to Ne.
            for i in range(a+1, Ne):
                coef = -gamma[a, b]*0.5
                if unfolding.real:
                    mur = Mu(1, i, a)
                    code += term_code(mur, mur, coef, matrix_form, False)
                    mui = Mu(-1, i, a)
                    code += term_code(mui, mui, coef, matrix_form, False)
                else:
                    mu = Mu(0, i, a)
                    code += term_code(mu, mu, coef, matrix_form, False)

            # We missed one term in each of the previous fors, that together
            # correspond to
            # -gamma_ab * rho_aa |a><a|
            coef = -gamma[a, b]
            if unfolding.real:
                mu = Mu(1, a, a)
            else:
                mu = Mu(0, a, a)
            code += term_code(mu, mu, coef, matrix_form, False)

    # We finish the code.
    if True:
        if matrix_form:
            if unfolding.normalized:
                code += "    return A, b\n"
            else:
                code += "    return A\n"
        else:
            code += "    return rhs\n"
    # We write the code to file if provided, and execute it.
    if True:
        if file_name is not None:
            f = file(file_name+".py", "w")
            f.write(code)
            f.close()

        lindblad_terms = code
        if not return_code:
            exec lindblad_terms
    return lindblad_terms


def fast_hamiltonian_terms(Ep, epsilonp, detuning_knob,
                           omega_level, rm, xi, theta,
                           unfolding, matrix_form=False, file_name=None,
                           return_code=False):
    r"""Return a fast function that returns the Hamiltonian terms.

    We test a basic two-level system.

    >>> import numpy as np
    >>> from scipy.constants import physical_constants
    >>> from sympy import Matrix, symbols
    >>> from fast.electric_field import electric_field_amplitude_top
    >>> from fast.symbolic import define_laser_variables, polarization_vector

    >>> Ne = 2
    >>> Nl = 1
    >>> a0 = physical_constants["Bohr radius"][0]
    >>> rm = [np.array([[0, 0], [a0, 0]]),
    ...       np.array([[0, 0], [0, 0]]),
    ...       np.array([[0, 0], [0, 0]])]
    >>> xi = np.array([[[0, 1], [1, 0]]])
    >>> omega_level = [0, 1.0e9]
    >>> theta = phase_transformation(Ne, Nl, rm, xi)

    We define symbolic variables to be used as token arguments.
    >>> Ep, omega_laser = define_laser_variables(Nl)
    >>> epsilonps = [polarization_vector(0, 0, 0, 0, 1)]
    >>> detuning_knob = [symbols("delta1", real=True)]

    An map to unfold the density matrix.
    >>> unfolding = Unfolding(Ne, True, True, True)

    We obtain a function to calculate Hamiltonian terms.
    >>> aux = (Ep, epsilonps, detuning_knob, omega_level, rm, xi, theta,
    ...        unfolding, False, None)
    >>> hamiltonian_terms = fast_hamiltonian_terms(*aux)

    Apply this to a density matrix.
    >>> rhos = np.array([[0.6, 3+2j],
    ...                  [3-2j, 0.4]])
    >>> rhosv = unfolding(rhos)

    We specify values for the variables
    >>> detuning_knobs = [100e6]
    >>> Eps = electric_field_amplitude_top(1e-3, 1e-3, 1, "SI")
    >>> Eps *= np.exp(1j*np.pi)
    >>> Eps = [Eps]
    >>> print hamiltonian_terms(rhosv, Eps, detuning_knobs)
    [  5.56808315e+07   2.00000000e+08   2.97215958e+08]


    """
    if not unfolding.lower_triangular:
        mes = "It is very inefficient to solve using all components of the "
        mes += "density matrix. Better set lower_triangular=True in Unfolding."
        raise NotImplementedError(mes)
    if matrix_form and (not unfolding.real) and (unfolding.lower_triangular):
        mes = "It is not possible to express the equations in matrix form "
        mes += "for complex lower triangular components only."
        raise ValueError(mes)
    Nl = len(Ep)
    Nrho = unfolding.Nrho
    # We determine which arguments are constants.
    if True:
        try:
            Ep = np.array([complex(Ep[l]) for l in range(Nl)])
            variable_Ep = False
        except:
            variable_Ep = True

        try:
            epsilonp = [np.array([complex(epsilonp[l][i]) for i in range(3)])
                        for l in range(Nl)]
            variable_epsilonp = False
        except:
            variable_epsilonp = True
        try:
            detuning_knob = np.array([float(detuning_knob[l])
                                      for l in range(Nl)])
            variable_detuning_knob = False
        except:
            variable_detuning_knob = True
    # We obtain code for the two parts.
    if True:
        if file_name is not None:
            file_name_rabi = file_name+"_rabi"
            file_name_detuning = file_name+"_detuning"
        else:
            file_name_rabi = file_name
            file_name_detuning = file_name

        rabi_terms = fast_rabi_terms(Ep, epsilonp, rm, xi, theta, unfolding,
                                     matrix_form=matrix_form,
                                     file_name=file_name_rabi,
                                     return_code=True)

        detuning_terms = fast_detuning_terms(detuning_knob, omega_level, xi,
                                             theta, unfolding,
                                             matrix_form=matrix_form,
                                             file_name=file_name_detuning,
                                             return_code=True)
        code = rabi_terms + "\n\n" + detuning_terms + "\n\n"

        # If these functions have 0 arguments, we call them only once!
        if not variable_Ep and not variable_epsilonp and matrix_form:
            code += "rabi_terms = rabi_terms()\n\n"
        if not variable_detuning_knob and matrix_form:
            code += "detuning_terms = detuning_terms()\n\n"
    # We establish the arguments of the output function.
    if True:
        code += "def hamiltonian_terms("
        if not matrix_form: code += "rho, "
        if variable_Ep: code += "Ep, "
        if variable_epsilonp: code += "epsilonp, "
        if variable_detuning_knob: code += "detuning_knob, "
        code += "rabi_terms=rabi_terms, detuning_terms=detuning_terms"
        # if code[-2:] == ", ": code = code[:-2]
        code += "):\n"

        code += '    r"""A fast calculation of the hamiltonian terms."""\n'
        # if not variable_Ep and not varia
    # We initialize the output and auxiliaries.
    if True:
        # We introduce the factor that multiplies all terms.
        if matrix_form:
            code += "    A = np.zeros(("+str(Nrho)+", "+str(Nrho)
            if not unfolding.real:
                code += "), complex)\n\n"
            else:
                code += "))\n\n"
            if unfolding.normalized:
                code += "    b = np.zeros(("+str(Nrho)
                if not unfolding.real:
                    code += "), complex)\n\n"
                else:
                    code += "))\n\n"
        else:
            code += "    rhs = np.zeros(("+str(Nrho)
            if not unfolding.real:
                code += "), complex)\n\n"
            else:
                code += "))\n\n"
    # We call the Rabi terms.
    if True:
        if not variable_Ep and not variable_epsilonp and matrix_form:
            aux_code = "rabi_terms\n"
        else:
            aux_code = "rabi_terms("
            if not matrix_form: aux_code += "rho, "
            if variable_Ep: aux_code += "Ep, "
            if variable_epsilonp: aux_code += "epsilonp, "
            if aux_code[-2:] == ", ": aux_code = aux_code[:-2]
            aux_code += ")\n"

        if matrix_form:
            if unfolding.normalized:
                code += "    aux = " + aux_code
                code += "    A += aux[0]\n"
                code += "    b += aux[1]\n"
            else:
                code += "    A = " + aux_code
        else:
            code += "    rhs = " + aux_code
    # We call the detuning terms.
    if True:
        if not variable_detuning_knob and matrix_form:
            aux_code = "detuning_terms\n"
        else:
            aux_code = "detuning_terms("
            if not matrix_form: aux_code += "rho, "
            if variable_detuning_knob: aux_code += "detuning_knob, "
            if aux_code[-2:] == ", ": aux_code = aux_code[:-2]
            aux_code += ")\n"

        if matrix_form:
            if unfolding.normalized:
                code += "    aux = " + aux_code
                code += "    A += aux[0]\n"
                code += "    b += aux[1]\n"
            else:
                code += "    A += " + aux_code
        else:
            code += "    rhs += " + aux_code
    # We finish the code.
    if True:
        # code = rabi_code + "\n\n" + code
        if matrix_form:
            if unfolding.normalized:
                code += "    return A, b\n"
            else:
                code += "    return A\n"
        else:
            code += "    return rhs\n"
    # We write the code to file if provided, and execute it.
    if True:
        if file_name is not None:
            f = file(file_name+".py", "w")
            f.write(code)
            f.close()

        hamiltonian_terms = code
        if not return_code:
            exec hamiltonian_terms
    return hamiltonian_terms


def fast_bloch_equations(Ep, epsilonp, detuning_knob, gamma,
                         omega_level, rm, xi, theta,
                         unfolding, matrix_form=False, file_name=None,
                         return_code=False):
    r"""Return a fast function that returns the numeric right-hand sides of \
    Bloch equations.

    We test a basic two-level system.

    >>> import numpy as np
    >>> from scipy.constants import physical_constants
    >>> from sympy import Matrix, symbols
    >>> from fast.electric_field import electric_field_amplitude_top
    >>> from fast.symbolic import (define_laser_variables,
    ...                            polarization_vector)

    >>> Ne = 2
    >>> Nl = 1
    >>> a0 = physical_constants["Bohr radius"][0]
    >>> rm = [np.array([[0, 0], [a0, 0]]),
    ...       np.array([[0, 0], [0, 0]]),
    ...       np.array([[0, 0], [0, 0]])]
    >>> xi = np.array([[[0, 1], [1, 0]]])
    >>> omega_level = [0, 1.0e9]
    >>> gamma21 = 2*np.pi*6e6
    >>> gamma = np.array([[0, -gamma21], [gamma21, 0]])
    >>> theta = phase_transformation(Ne, Nl, rm, xi)

    We define symbolic variables to be used as token arguments.
    >>> Ep, omega_laser = define_laser_variables(Nl)
    >>> epsilonps = [polarization_vector(0, 0, 0, 0, 1)]
    >>> detuning_knob = [symbols("delta1", real=True)]

    An map to unfold the density matrix.
    >>> unfolding = Unfolding(Ne, True, True, True)

    We obtain a function to calculate the Bloch equations.
    >>> aux = (Ep, epsilonps, detuning_knob, gamma,
    ...        omega_level, rm, xi, theta,
    ...        unfolding, False, None)
    >>> bloch_equations = fast_bloch_equations(*aux)

    Apply this to a density matrix.
    >>> rhos = np.array([[0.6, 3+2j],
    ...                  [3-2j, 0.4]])
    >>> rhosv = unfolding(rhos)

    We specify values for the variables
    >>> detuning_knobs = [100e6]
    >>> Eps = electric_field_amplitude_top(1e-3, 1e-3, 1, "SI")
    >>> Eps *= np.exp(1j*np.pi)
    >>> Eps = [Eps]
    >>> print bloch_equations(rhosv, Eps, detuning_knobs)
    [  4.06011867e+07   1.43451332e+08   3.34915070e+08]

    """
    if not unfolding.lower_triangular:
        mes = "It is very inefficient to solve using all components of the "
        mes += "density matrix. Better set lower_triangular=True in Unfolding."
        raise NotImplementedError(mes)
    if matrix_form and (not unfolding.real) and (unfolding.lower_triangular):
        mes = "It is not possible to express the equations in matrix form "
        mes += "for complex lower triangular components only."
        raise ValueError(mes)
    Nl = len(Ep)
    Nrho = unfolding.Nrho
    # We determine which arguments are constants.
    if True:
        try:
            Ep = np.array([complex(Ep[l]) for l in range(Nl)])
            variable_Ep = False
        except:
            variable_Ep = True

        try:
            epsilonp = [np.array([complex(epsilonp[l][i]) for i in range(3)])
                        for l in range(Nl)]
            variable_epsilonp = False
        except:
            variable_epsilonp = True
        try:
            detuning_knob = np.array([float(detuning_knob[l])
                                      for l in range(Nl)])
            variable_detuning_knob = False
        except:
            variable_detuning_knob = True
    # We obtain code for the three parts.
    if True:
        if file_name is not None:
            file_name_rabi = file_name+"_rabi"
            file_name_detuning = file_name+"_detuning"
            file_name_lindblad = file_name+"_lindblad"
        else:
            file_name_rabi = file_name
            file_name_detuning = file_name
            file_name_lindblad = file_name

        rabi_terms = fast_rabi_terms(Ep, epsilonp, rm, xi, theta, unfolding,
                                     matrix_form=matrix_form,
                                     file_name=file_name_rabi,
                                     return_code=True)

        detuning_terms = fast_detuning_terms(detuning_knob, omega_level, xi,
                                             theta, unfolding,
                                             matrix_form=matrix_form,
                                             file_name=file_name_detuning,
                                             return_code=True)

        lindblad_terms = fast_lindblad_terms(gamma, unfolding,
                                             matrix_form=matrix_form,
                                             file_name=file_name_lindblad,
                                             return_code=True)

        code = rabi_terms+"\n\n"
        code += detuning_terms+"\n\n"
        code += lindblad_terms+"\n\n"

        # If these functions have 0 arguments, we call them only once!
        if not variable_Ep and not variable_epsilonp and matrix_form:
            code += "rabi_terms = rabi_terms()\n\n"
        if not variable_detuning_knob and matrix_form:
            code += "detuning_terms = detuning_terms()\n\n"
        if matrix_form:
            code += "lindblad_terms = lindblad_terms()\n\n"
    # We establish the arguments of the output function.
    if True:
        code += "def bloch_equations("
        if not matrix_form: code += "rho, "
        if variable_Ep: code += "Ep, "
        if variable_epsilonp: code += "epsilonp, "
        if variable_detuning_knob: code += "detuning_knob, "
        code += "rabi_terms=rabi_terms, detuning_terms=detuning_terms, "
        code += "lindblad_terms=lindblad_terms):\n"
        code += '    r"""A fast calculation of Bloch equations."""\n'
    # We initialize the output and auxiliaries.
    if True:
        # We introduce the factor that multiplies all terms.
        if matrix_form:
            code += "    A = np.zeros(("+str(Nrho)+", "+str(Nrho)
            if not unfolding.real:
                code += "), complex)\n\n"
            else:
                code += "))\n\n"
            if unfolding.normalized:
                code += "    b = np.zeros(("+str(Nrho)
                if not unfolding.real:
                    code += "), complex)\n\n"
                else:
                    code += "))\n\n"
        else:
            code += "    rhs = np.zeros(("+str(Nrho)
            if not unfolding.real:
                code += "), complex)\n\n"
            else:
                code += "))\n\n"
    # We call the Rabi terms.
    if True:
        if not variable_Ep and not variable_epsilonp and matrix_form:
            aux_code = "rabi_terms\n"
        else:
            aux_code = "rabi_terms("
            if not matrix_form: aux_code += "rho, "
            if variable_Ep: aux_code += "Ep, "
            if variable_epsilonp: aux_code += "epsilonp, "
            if aux_code[-2:] == ", ": aux_code = aux_code[:-2]
            aux_code += ")\n"

        if matrix_form:
            if unfolding.normalized:
                code += "    aux = " + aux_code
                code += "    A += aux[0]\n"
                code += "    b += aux[1]\n"
            else:
                code += "    A = " + aux_code
        else:
            code += "    rhs = " + aux_code
    # We call the detuning terms.
    if True:
        if not variable_detuning_knob and matrix_form:
            aux_code = "detuning_terms\n"
        else:
            aux_code = "detuning_terms("
            if not matrix_form: aux_code += "rho, "
            if variable_detuning_knob: aux_code += "detuning_knob, "
            if aux_code[-2:] == ", ": aux_code = aux_code[:-2]
            aux_code += ")\n"

        if matrix_form:
            if unfolding.normalized:
                code += "    aux = " + aux_code
                code += "    A += aux[0]\n"
                code += "    b += aux[1]\n"
            else:
                code += "    A += " + aux_code
        else:
            code += "    rhs += " + aux_code
    # We call the Lindblad terms.
    if True:
        if matrix_form:
            aux_code = "lindblad_terms\n"
        else:
            aux_code = "lindblad_terms(rho)\n"
        if matrix_form:
            if unfolding.normalized:
                code += "    aux = " + aux_code
                code += "    A += aux[0]\n"
                code += "    b += aux[1]\n"
            else:
                code += "    A += " + aux_code
        else:
            code += "    rhs += " + aux_code
    # We finish the code.
    if True:
        # code = rabi_code + "\n\n" + code
        if matrix_form:
            if unfolding.normalized:
                code += "    return A, b\n"
            else:
                code += "    return A\n"
        else:
            code += "    return rhs\n"
    # We write the code to file if provided, and execute it.
    if True:
        if file_name is not None:
            f = file(file_name+".py", "w")
            f.write(code)
            f.close()

        bloch_equations = code
        if not return_code:
            exec bloch_equations
    return bloch_equations


def fast_steady_state(Ep, epsilonp, detuning_knob, gamma,
                      omega_level, rm, xi, theta,
                      file_name=None, return_code=False):
    r"""Return a fast function that returns a steady state.

    We test a basic two-level system.

    >>> import numpy as np
    >>> from scipy.constants import physical_constants
    >>> from sympy import Matrix, symbols
    >>> from fast.electric_field import electric_field_amplitude_top
    >>> from fast.symbolic import (define_laser_variables,
    ...                            polarization_vector)

    >>> Ne = 2
    >>> Nl = 1
    >>> a0 = physical_constants["Bohr radius"][0]
    >>> rm = [np.array([[0, 0], [a0, 0]]),
    ...       np.array([[0, 0], [0, 0]]),
    ...       np.array([[0, 0], [0, 0]])]
    >>> xi = np.array([[[0, 1], [1, 0]]])
    >>> omega_level = [0, 1.0e9]
    >>> gamma21 = 2*np.pi*6e6
    >>> gamma = np.array([[0, -gamma21], [gamma21, 0]])
    >>> theta = phase_transformation(Ne, Nl, rm, xi)

    We define symbolic variables to be used as token arguments.
    >>> Ep, omega_laser = define_laser_variables(Nl)
    >>> epsilonps = [polarization_vector(0, 0, 0, 0, 1)]
    >>> detuning_knob = [symbols("delta1", real=True)]

    An map to unfold the density matrix.
    >>> unfolding = Unfolding(Ne, True, True, True)

    We obtain a function to calculate Hamiltonian terms.
    >>> aux = (Ep, epsilonps, detuning_knob, gamma,
    ...        omega_level, rm, xi, theta)
    >>> steady_state = fast_steady_state(*aux)

    We specify values for the variables
    >>> detuning_knobs = [100e6]
    >>> Eps = electric_field_amplitude_top(1e-3, 1e-3, 1, "SI")
    >>> Eps *= np.exp(1j*np.pi)
    >>> Eps = [Eps]
    >>> print steady_state(Eps, detuning_knobs)
    [ 0.01803732  0.12957649 -0.02442459]

    """
    # We unpack variables.
    if True:
        Ne = len(omega_level)
        Nl = xi.shape[0]
        unfolding = Unfolding(Ne, True, True, True)
    # We determine which arguments are constants.
    if True:
        try:
            Ep = np.array([complex(Ep[l]) for l in range(Nl)])
            variable_Ep = False
        except:
            variable_Ep = True

        try:
            epsilonp = [np.array([complex(epsilonp[l][i]) for i in range(3)])
                        for l in range(Nl)]
            variable_epsilonp = False
        except:
            variable_epsilonp = True
        try:
            detuning_knob = np.array([float(detuning_knob[l])
                                      for l in range(Nl)])
            variable_detuning_knob = False
        except:
            variable_detuning_knob = True
    # We obtain code for the three parts.
    if True:
        args = (Ep, epsilonp, detuning_knob, gamma,
                omega_level, rm, xi, theta,
                unfolding, True, None, True)

        bloch_equations = fast_bloch_equations(*args)

        code = bloch_equations+"\n\n"

        if ((not variable_Ep) and
           (not variable_epsilonp) and
           (not variable_detuning_knob)):
            # We can call bloch_equations here!
            code += "bloch_equations = bloch_equations()\n"
    # We establish the arguments of the output function.
    if True:
        code += "def steady_state("
        if variable_Ep: code += "Ep, "
        if variable_epsilonp: code += "epsilonp, "
        if variable_detuning_knob: code += "detuning_knob, "
        code += "bloch_equations=bloch_equations):\n"
        code += '    r"""A fast calculation of the steady state."""\n'
    # We call the Bloch equations.
    if True:
        code += r"""    A, b = bloch_equations"""
        if ((not variable_Ep) and
           (not variable_epsilonp) and
           (not variable_detuning_knob)):

            code += "\n"
        else:
            code += "("
            if variable_Ep: code += "Ep, "
            if variable_epsilonp: code += "epsilonp, "
            if variable_detuning_knob: code += "detuning_knob, "
            if code[-2:] == ", ": code = code[:-2]
            code += ")\n"

        code += """    rhox = np.linalg.solve(A, b)\n"""
        code += """    return rhox\n"""
    # We write the code to file if provided, and execute it.
    if True:
        if file_name is not None:
            f = file(file_name+".py", "w")
            f.write(code)
            f.close()

        steady_state = code
        if not return_code:
            exec steady_state
    return steady_state


def fast_time_evolution(Ep, epsilonp, detuning_knob, gamma,
                        omega_level, rm, xi, theta,
                        semi_analytic=True,
                        file_name=None, return_code=False):
    r"""Return a fast function to calculate the time evolution of a state.

    We test a basic two-level system.

    >>> import numpy as np
    >>> from scipy.constants import physical_constants
    >>> from sympy import Matrix, symbols
    >>> from fast.electric_field import electric_field_amplitude_top
    >>> from fast.symbolic import (define_laser_variables,
    ...                            polarization_vector)

    >>> Ne = 2
    >>> Nl = 1
    >>> a0 = physical_constants["Bohr radius"][0]
    >>> rm = [np.array([[0, 0], [a0, 0]]),
    ...       np.array([[0, 0], [0, 0]]),
    ...       np.array([[0, 0], [0, 0]])]
    >>> xi = np.array([[[0, 1], [1, 0]]])
    >>> omega_level = [0, 1.0e9]
    >>> gamma21 = 2*np.pi*6e6
    >>> gamma = np.array([[0, -gamma21], [gamma21, 0]])
    >>> theta = phase_transformation(Ne, Nl, rm, xi)

    We define symbolic variables to be used as token arguments.
    >>> Ep, omega_laser = define_laser_variables(Nl)
    >>> epsilonp = [polarization_vector(0, 0, 0, 0, 1)]
    >>> detuning_knob = [symbols("delta1", real=True)]

    An map to unfold the density matrix.
    >>> unfolding = Unfolding(Ne, True, True, True)

    We obtain a function to calculate time evolution.
    >>> aux = (Ep, epsilonp, detuning_knob, gamma,
    ...        omega_level, rm, xi, theta)
    >>> time_evolution = fast_time_evolution(*aux)

    We specify values for the variables
    >>> detuning_knobs = [100e6]
    >>> Eps = electric_field_amplitude_top(1e-3, 1e-3, 1, "SI")
    >>> Eps *= np.exp(1j*np.pi)
    >>> Eps = [Eps]

    >>> t = np.linspace(0, 1e-6, 11)
    >>> rho0 = np.array([[1, 0], [0, 0]])
    >>> rho0 = unfolding(rho0)

    >>> print time_evolution(t, rho0, Eps, detuning_knobs)
    [[ 0.          0.          0.        ]
     [ 0.02147019  0.14276917 -0.01160453]
     [ 0.01825732  0.12987983 -0.02707568]
     [ 0.01794253  0.1292563  -0.0242422 ]
     [ 0.01804897  0.12962366 -0.02440144]
     [ 0.0180373   0.12957509 -0.02443214]
     [ 0.01803708  0.12957577 -0.0244238 ]
     [ 0.01803736  0.12957664 -0.02442457]
     [ 0.01803731  0.12957648 -0.02442461]
     [ 0.01803732  0.12957649 -0.02442459]
     [ 0.01803732  0.12957649 -0.02442459]]

    >>> print time_evolution(t, rho0, Eps, detuning_knobs, average=True)
    [ 0.0174924   0.12441977 -0.02216659]

    """
    if not semi_analytic:
        s = "The numeric integrator has not been implemented."
        raise NotImplementedError(s)
    # We unpack variables.
    if True:
        Ne = len(omega_level)
        Nl = xi.shape[0]
        unfolding = Unfolding(Ne, True, True, True)
        Nrho = unfolding.Nrho
    # We determine which arguments are constants.
    if True:
        try:
            Ep = np.array([complex(Ep[l]) for l in range(Nl)])
            variable_Ep = False
        except:
            variable_Ep = True

        try:
            epsilonp = [np.array([complex(epsilonp[l][i]) for i in range(3)])
                        for l in range(Nl)]
            variable_epsilonp = False
        except:
            variable_epsilonp = True
        try:
            detuning_knob = np.array([float(detuning_knob[l])
                                      for l in range(Nl)])
            variable_detuning_knob = False
        except:
            variable_detuning_knob = True
    # We obtain code for the three parts.
    if True:
        args = (Ep, epsilonp, detuning_knob, gamma,
                omega_level, rm, xi, theta,
                unfolding, True, None, True)

        bloch_equations = fast_bloch_equations(*args)

        code = bloch_equations+"\n\n"

        if ((not variable_Ep) and
           (not variable_epsilonp) and
           (not variable_detuning_knob)):
            # We can call bloch_equations here!
            code += "bloch_equations = bloch_equations()\n"
    # We establish the arguments of the output function.
    if True:
        code += "def time_evolution(t, rho0, "
        if variable_Ep: code += "Ep, "
        if variable_epsilonp: code += "epsilonp, "
        if variable_detuning_knob: code += "detuning_knob, "
        code += "average=False, "
        code += "bloch_equations=bloch_equations):\n"
        code += '    r"""A fast calculation of time evolution."""\n'
    # We call the Bloch equations.
    if True:
        code += r"""    A, b = bloch_equations"""
        if ((not variable_Ep) and
           (not variable_epsilonp) and
           (not variable_detuning_knob)):

            code += "\n"
        else:
            code += "("
            if variable_Ep: code += "Ep, "
            if variable_epsilonp: code += "epsilonp, "
            if variable_detuning_knob: code += "detuning_knob, "
            if code[-2:] == ", ": code = code[:-2]
            code += ")\n"

        code += "    Nt = t.shape[0]\n"
        code += "    lam, S = np.linalg.eig(A)\n"
        code += "    Sinv = np.linalg.inv(S)\n"
        code += "    d = np.dot(Sinv, b)/lam\n"
        code += "    r = np.dot(Sinv, rho0) - d\n"
        code += "    rho_steady = np.dot(S, d)\n"
        code += "    rho = np.zeros((Nt, %s))\n" % Nrho
        code += "    for i, ti in enumerate(t):\n"
        code += "        rho_prime = r*np.exp(lam*ti)\n"
        code += "        rho[i, :] = np.real(np.dot(S, rho_prime) "
        code += "+ rho_steady)\n"
        code += "    if average:\n"
        code += "        rho = time_average(rho, t)\n"
        code += """    return rho\n"""
    # We write the code to file if provided, and execute it.
    if True:
        if file_name is not None:
            f = file(file_name+".py", "w")
            f.write(code)
            f.close()

        time_evolution = code
        if not return_code:
            exec time_evolution
    return time_evolution


def time_average(rho, t):
    r"""Return a time-averaged density matrix (using trapezium rule).

    We test a basic two-level system.

    >>> import numpy as np
    >>> from scipy.constants import physical_constants
    >>> from sympy import Matrix, symbols
    >>> from fast.electric_field import electric_field_amplitude_top
    >>> from fast.symbolic import (define_laser_variables,
    ...                            polarization_vector)

    >>> Ne = 2
    >>> Nl = 1
    >>> a0 = physical_constants["Bohr radius"][0]
    >>> rm = [np.array([[0, 0], [a0, 0]]),
    ...       np.array([[0, 0], [0, 0]]),
    ...       np.array([[0, 0], [0, 0]])]
    >>> xi = np.array([[[0, 1], [1, 0]]])
    >>> omega_level = [0, 1.0e9]
    >>> gamma21 = 2*np.pi*6e6
    >>> gamma = np.array([[0, -gamma21], [gamma21, 0]])
    >>> theta = phase_transformation(Ne, Nl, rm, xi)

    We define symbolic variables to be used as token arguments.
    >>> Ep, omega_laser = define_laser_variables(Nl)
    >>> epsilonps = [polarization_vector(0, 0, 0, 0, 1)]
    >>> detuning_knob = [symbols("delta1", real=True)]

    An map to unfold the density matrix.
    >>> unfolding = Unfolding(Ne, True, True, True)

    We obtain a function to calculate Hamiltonian terms.
    >>> aux = (Ep, epsilonps, detuning_knob, gamma,
    ...        omega_level, rm, xi, theta)
    >>> time_evolution = fast_time_evolution(*aux)

    We specify values for the variables
    >>> detuning_knobs = [100e6]
    >>> Eps = electric_field_amplitude_top(1e-3, 1e-3, 1, "SI")
    >>> Eps *= np.exp(1j*np.pi)
    >>> Eps = [Eps]

    >>> t = np.linspace(0, 1e-6, 11)
    >>> rho0 = np.array([[1, 0], [0, 0]])
    >>> rho0 = unfolding(rho0)

    >>> rho = time_evolution(t, rho0, Eps, detuning_knobs)
    >>> print time_average(rho, t)
    [ 0.0174924   0.12441977 -0.02216659]

    """
    T = t[-1]-t[0]
    dt = t[1]-t[0]
    rhoav = np.sum(rho[1:-1], axis=0) + 0.5*(rho[0]+rho[-1])
    rhoav = dt/T*rhoav
    return rhoav


def fast_sweep_steady_state(Ep, epsilonp, gamma,
                            omega_level, rm, xi, theta,
                            file_name=None, return_code=False):
    r"""Return an spectrum of density matrices in the steady state.

    We test a basic two-level system.

    >>> import numpy as np
    >>> from sympy import symbols
    >>> from scipy.constants import physical_constants

    >>> e_num = physical_constants["elementary charge"][0]
    >>> hbar_num = physical_constants["Planck constant over 2 pi"][0]

    >>> Ne = 2
    >>> Nl = 1
    >>> Ep = [-1.0]
    >>> epsilonp = [np.array([0, 0, 1.0])]
    >>> delta = symbols("delta")

    >>> detuning_knob = [delta]
    >>> gamma = np.array([[0.0, -1.0], [1.0, 0.0]])
    >>> omega_level = np.array([0.0, 100.0])
    >>> rm = [np.array([[0.0, 0.0], [1.0, 0.0]])*hbar_num/e_num
    ...       for p in range(3)]
    >>> xi = np.array([[[0, 1], [1, 0]]])
    >>> theta = phase_transformation(Ne, Nl, rm, xi)
    >>> sweep_steady_state = fast_sweep_steady_state(Ep, epsilonp, gamma,
    ...                                              omega_level, rm, xi,
    ...                                              theta)
    >>> deltas, rho = sweep_steady_state([[-20, 20, 11]])
    >>> print rho
    [[ 0.00062383 -0.02495321 -0.00062383]
     [ 0.00097371 -0.03115871 -0.00097371]
     [ 0.00172712 -0.04145078 -0.00172712]
     [ 0.003861   -0.06177606 -0.003861  ]
     [ 0.01492537 -0.11940299 -0.01492537]
     [ 0.33333333  0.         -0.33333333]
     [ 0.01492537  0.11940299 -0.01492537]
     [ 0.003861    0.06177606 -0.003861  ]
     [ 0.00172712  0.04145078 -0.00172712]
     [ 0.00097371  0.03115871 -0.00097371]
     [ 0.00062383  0.02495321 -0.00062383]]

    """
    # We unpack variables.
    if True:
        Nl = xi.shape[0]
    # We determine which arguments are constants.
    if True:
        try:
            Ep = np.array([complex(Ep[l]) for l in range(Nl)])
            variable_Ep = False
        except:
            variable_Ep = True

        try:
            epsilonp = [np.array([complex(epsilonp[l][i]) for i in range(3)])
                        for l in range(Nl)]
            variable_epsilonp = False
        except:
            variable_epsilonp = True
    # We obtain code for the steady state.
    if True:
        detuning_knob = symbols("delta1:"+str(Nl))
        args = (Ep, epsilonp, detuning_knob, gamma, omega_level, rm, xi, theta,
                file_name, True)

        steady_state = fast_steady_state(*args)
        code = steady_state+"\n\n"
    # We establish the arguments of the output function.
    if True:
        code += "def sweep_steady_state("
        if variable_Ep: code += "Ep, "
        if variable_epsilonp: code += "epsilonp, "
        code += "detuning_knob, "
        code += "steady_state=steady_state):\n"
        code += '    r"""A fast frequency sweep of the steady state."""\n'
    # Code to determine the sweep range.
    if True:
        code += """    sweepN = -1\n"""
        code += """    for i, delta in enumerate(detuning_knob):\n"""
        code += """        if hasattr(delta, "__getitem__"):\n"""
        code += """            sweepN = i\n"""
        code += """            delta0 = delta[0]\n"""
        code += """            deltaf = delta[1]\n"""
        code += """            Ndelta = delta[2]\n"""
        code += """            break\n\n"""
        code += """    if sweepN == -1:\n"""
        code += """        s = 'One of the detuning knobs '\n"""
        code += """        s += 'must be of the form '\n"""
        code += """        s += '(start, stop, Nsteps)'\n"""
        code += """        raise ValueError(s)\n\n"""
        code += """    deltas = np.linspace(delta0, deltaf, Ndelta)\n\n"""
    # We call steady_state.
    if True:
        code += "    args = [["
        if variable_Ep: code += "Ep, "
        if variable_epsilonp: code += "epsilonp, "
        code += """list(detuning_knob[:sweepN]) +\n"""
        code += """            [deltas[i]] +\n"""
        code += """            list(detuning_knob[sweepN+1:])]\n"""
        code += """          for i in range(Ndelta)]\n\n"""
        code += "    rho = np.array([steady_state(*argsi)\n"
        code += "                   for argsi in args])\n\n"
    # We finish the code.
    if True:
        code += "    return deltas, rho\n"
    # We write the code to file if provided, and execute it.
    if True:
        if file_name is not None:
            f = file(file_name+".py", "w")
            f.write(code)
            f.close()

        sweep_steady_state = code
        if not return_code:
            exec sweep_steady_state
    return sweep_steady_state


def fast_sweep_time_evolution(Ep, epsilonp, gamma,
                              omega_level, rm, xi, theta,
                              semi_analytic=True,
                              file_name=None, return_code=False):
    r"""Return a spectrum of time evolutions of the density matrix.

    We test a basic two-level system.

    >>> import numpy as np
    >>> from sympy import symbols
    >>> from scipy.constants import physical_constants

    >>> e_num = physical_constants["elementary charge"][0]
    >>> hbar_num = physical_constants["Planck constant over 2 pi"][0]

    >>> Ne = 2
    >>> Nl = 1
    >>> Ep = [-1.0]
    >>> epsilonp = [np.array([0, 0, 1.0])]
    >>> delta = symbols("delta")

    >>> detuning_knob = [delta]
    >>> gamma = np.array([[0.0, -1.0], [1.0, 0.0]])
    >>> omega_level = np.array([0.0, 100.0])
    >>> rm = [np.array([[0.0, 0.0], [1.0, 0.0]])*hbar_num/e_num
    ...       for p in range(3)]
    >>> xi = np.array([[[0, 1], [1, 0]]])
    >>> theta = phase_transformation(Ne, Nl, rm, xi)
    >>> sweep_time_evolution = fast_sweep_time_evolution(Ep, epsilonp, gamma,
    ...                                                  omega_level, rm, xi,
    ...                                                  theta)

    >>> t = np.linspace(0, 1e1, 11)
    >>> unfolding = Unfolding(Ne, True, True, True)
    >>> rho0 = np.array([[1, 0], [0, 0]])
    >>> rho0 = unfolding(rho0)

    >>> deltas, rho = sweep_time_evolution(t, rho0, [[-20, 20, 5]])
    >>> print rho
    [[[  0.00000000e+00   0.00000000e+00   0.00000000e+00]
      [  5.62048975e-04  -1.87738708e-02  -1.44368294e-02]
      [  1.03023067e-03  -3.12259136e-02  -7.30305635e-03]
      [  9.12184547e-04  -3.01487805e-02   1.33250669e-03]
      [  6.37106612e-04  -2.50730799e-02   2.74370534e-03]
      [  5.34378513e-04  -2.30997060e-02   2.29774290e-04]
      [  5.80976828e-04  -2.40435159e-02  -1.46256227e-03]
      [  6.38084021e-04  -2.52093923e-02  -1.32909716e-03]
      [  6.46745698e-04  -2.54070480e-02  -6.44976127e-04]
      [  6.29482707e-04  -2.50709886e-02  -3.74574189e-04]
      [  6.18119824e-04  -2.48414541e-02  -4.99674790e-04]]
    <BLANKLINE>
     [[  0.00000000e+00   0.00000000e+00   0.00000000e+00]
      [  5.81424860e-03  -7.46503637e-02   1.38589101e-02]
      [  2.24583358e-03  -4.30271127e-02  -1.94356220e-02]
      [  2.27876199e-03  -4.68669200e-02   8.17093102e-03]
      [  3.05714698e-03  -5.47244372e-02  -6.72997591e-03]
      [  2.09801690e-03  -4.56261983e-02  -2.21208118e-03]
      [  2.68657227e-03  -5.16850335e-02  -1.19060077e-03]
      [  2.43506413e-03  -4.90723263e-02  -3.84674854e-03]
      [  2.45717014e-03  -4.94192610e-02  -1.61405170e-03]
      [  2.52412284e-03  -5.00355681e-02  -2.83273897e-03]
      [  2.44906189e-03  -4.93038867e-02  -2.45410426e-03]]
    <BLANKLINE>
     [[  0.00000000e+00   0.00000000e+00   0.00000000e+00]
      [  1.43610413e-01   0.00000000e+00  -3.44581796e-01]
      [  3.06127967e-01   0.00000000e+00  -4.13732743e-01]
      [  3.61099850e-01   0.00000000e+00  -3.73871185e-01]
      [  3.54270513e-01   0.00000000e+00  -3.37098714e-01]
      [  3.38348041e-01   0.00000000e+00  -3.26304189e-01]
      [  3.31346609e-01   0.00000000e+00  -3.28729595e-01]
      [  3.31150662e-01   0.00000000e+00  -3.32436798e-01]
      [  3.32607978e-01   0.00000000e+00  -3.33880540e-01]
      [  3.33431981e-01   0.00000000e+00  -3.33826365e-01]
      [  3.33548000e-01   0.00000000e+00  -3.33475215e-01]]
    <BLANKLINE>
     [[  0.00000000e+00   0.00000000e+00   0.00000000e+00]
      [  5.81424860e-03   7.46503637e-02   1.38589101e-02]
      [  2.24583358e-03   4.30271127e-02  -1.94356220e-02]
      [  2.27876199e-03   4.68669200e-02   8.17093102e-03]
      [  3.05714698e-03   5.47244372e-02  -6.72997591e-03]
      [  2.09801690e-03   4.56261983e-02  -2.21208118e-03]
      [  2.68657227e-03   5.16850335e-02  -1.19060077e-03]
      [  2.43506413e-03   4.90723263e-02  -3.84674854e-03]
      [  2.45717014e-03   4.94192610e-02  -1.61405170e-03]
      [  2.52412284e-03   5.00355681e-02  -2.83273897e-03]
      [  2.44906189e-03   4.93038867e-02  -2.45410426e-03]]
    <BLANKLINE>
     [[  0.00000000e+00   0.00000000e+00   0.00000000e+00]
      [  5.62048975e-04   1.87738708e-02  -1.44368294e-02]
      [  1.03023067e-03   3.12259136e-02  -7.30305635e-03]
      [  9.12184547e-04   3.01487805e-02   1.33250669e-03]
      [  6.37106612e-04   2.50730799e-02   2.74370534e-03]
      [  5.34378513e-04   2.30997060e-02   2.29774290e-04]
      [  5.80976828e-04   2.40435159e-02  -1.46256227e-03]
      [  6.38084021e-04   2.52093923e-02  -1.32909716e-03]
      [  6.46745698e-04   2.54070480e-02  -6.44976127e-04]
      [  6.29482707e-04   2.50709886e-02  -3.74574189e-04]
      [  6.18119824e-04   2.48414541e-02  -4.99674790e-04]]]

    >>> deltas, rho = sweep_time_evolution(t, rho0, [[-20, 20, 11]],
    ...                                    average=True)
    >>> print rho
    [[ 0.00064803 -0.0240473  -0.00214949]
     [ 0.00105421 -0.03076586 -0.00072189]
     [ 0.00157427 -0.0375314   0.00239996]
     [ 0.00412612 -0.06039333 -0.00605042]
     [ 0.01603272 -0.11753686 -0.01178172]
     [ 0.2998768   0.         -0.32911995]
     [ 0.01603272  0.11753686 -0.01178172]
     [ 0.00412612  0.06039333 -0.00605042]
     [ 0.00157427  0.0375314   0.00239996]
     [ 0.00105421  0.03076586 -0.00072189]
     [ 0.00064803  0.0240473  -0.00214949]]

    """
    # We unpack variables.
    if True:
        Nl = xi.shape[0]
    # We determine which arguments are constants.
    if True:
        try:
            Ep = np.array([complex(Ep[l]) for l in range(Nl)])
            variable_Ep = False
        except:
            variable_Ep = True

        try:
            epsilonp = [np.array([complex(epsilonp[l][i]) for i in range(3)])
                        for l in range(Nl)]
            variable_epsilonp = False
        except:
            variable_epsilonp = True
    # We obtain code for the steady state.
    if True:
        detuning_knob = symbols("delta1:"+str(Nl))
        args = (Ep, epsilonp, detuning_knob, gamma, omega_level, rm, xi, theta,
                file_name, True)

        args = (Ep, epsilonp, detuning_knob, gamma, omega_level, rm, xi,
                theta, True, file_name, True)
        time_evolution = fast_time_evolution(*args)
        code = time_evolution+"\n\n"
    # We establish the arguments of the output function.
    if True:
        code += "def sweep_time_evolution(t, rho0, "
        if variable_Ep: code += "Ep, "
        if variable_epsilonp: code += "epsilonp, "
        code += "detuning_knob, average=False, "
        code += "time_evolution=time_evolution):\n"
        code += '    r"""A fast frequency sweep of the steady state."""\n'
    # Code to determine the sweep range.
    if True:
        code += """    sweepN = -1\n"""
        code += """    for i, delta in enumerate(detuning_knob):\n"""
        code += """        if hasattr(delta, "__getitem__"):\n"""
        code += """            sweepN = i\n"""
        code += """            delta0 = delta[0]\n"""
        code += """            deltaf = delta[1]\n"""
        code += """            Ndelta = delta[2]\n"""
        code += """            break\n\n"""
        code += """    if sweepN == -1:\n"""
        code += """        s = 'One of the detuning knobs '\n"""
        code += """        s += 'must be of the form '\n"""
        code += """        s += '(start, stop, Nsteps)'\n"""
        code += """        raise ValueError(s)\n\n"""
        code += """    deltas = np.linspace(delta0, deltaf, Ndelta)\n\n"""
    # We call steady_state.
    if True:
        code += "    args = [[t, rho0, "
        if variable_Ep: code += "Ep, "
        if variable_epsilonp: code += "epsilonp, "
        code += """list(detuning_knob[:sweepN]) +\n"""
        code += """            [deltas[i]] +\n"""
        code += """            list(detuning_knob[sweepN+1:]), average]\n"""
        code += """          for i in range(Ndelta)]\n\n"""
        code += "    rho = np.array([time_evolution(*argsi)\n"
        code += "                   for argsi in args])\n\n"
    # We finish the code.
    if True:
        code += "    return deltas, rho\n"
    # We write the code to file if provided, and execute it.
    if True:
        if file_name is not None:
            f = file(file_name+".py", "w")
            f.write(code)
            f.close()

        sweep_time_evolution = code
        if not return_code:
            exec sweep_time_evolution
    return sweep_time_evolution


def observable(operator, rho, unfolding, complex=False):
    r"""Return an observable ammount.

    INPUT:

    -  ``operator`` - An square matrix representing a hermitian operator \
    in thesame basis as the density matrix.
    -  ``rho`` - A density matrix in unfolded format, or a list of such \
    density matrices.
    -  ``unfolding`` - A mapping from matrix element indices to unfolded \
    indices.

    >>> Ne = 2
    >>> unfolding = Unfolding(Ne, True, True, True)
    >>> rho = np.array([[0.6, 1+2j], [1-2j, 0.4]])
    >>> rho = unfolding(rho)
    >>> sx = np.array([[0, 1], [1, 0]])
    >>> print observable(sx, rho, unfolding)
    2.0

    """
    if len(rho.shape) == 2:
        return np.array([observable(operator, i, unfolding) for i in rho])

    Ne = unfolding.Ne
    Mu = unfolding.Mu
    obs = 0

    if unfolding.normalized:
        rho11 = 1 - sum([rho[Mu(1, i, i)] for i in range(1, Ne)])

    for i in range(Ne):
        for k in range(Ne):
            if unfolding.real:
                if k == 0 and i == 0:
                    obs += operator[i, k]*rho11
                else:
                    if k < i:
                        u, v = (i, k)
                    else:
                        u, v = (k, i)
                    obs += operator[i, k]*rho[Mu(1, u, v)]
                    if k != i:
                        if k < i:
                            obs += 1j*operator[i, k]*rho[Mu(-1, u, v)]
                        else:
                            obs += -1j*operator[i, k]*rho[Mu(-1, u, v)]

            else:
                if k == 0 and i == 0:
                    obs += operator[i, k]*rho11
                else:
                    obs += operator[i, k]*rho[Mu(0, k, i)]

    if not complex:
        obs = np.real(obs)
    return obs


def electric_succeptibility(l, Ep, epsilonp, rm, n, rho, unfolding, part=0):
    r"""Return the electric succeptibility for a given field.

    INPUT:

    -  ``l`` - The index labeling the probe field.
    -  ``Ep`` - A list of the amplitudes of all pump fields.
    -  ``epsilonp`` - The polarization vector of the probe field.
    -  ``rm`` -     The below-diagonal components of the position operator \
    in the cartesian basis:
    -  ``n`` - The number density of atoms.
    -  ``rho`` - A density matrix in unfolded format, or a list of such \
    density matrices.
    -  ``unfolding`` - A mapping from matrix element indices to unfolded \
    indices.


    >>> import numpy as np
    >>> from sympy import symbols
    >>> from scipy.constants import physical_constants
    >>> from fast import vapour_number_density

    >>> e_num = physical_constants["elementary charge"][0]
    >>> hbar_num = physical_constants["Planck constant over 2 pi"][0]

    >>> Ne = 2
    >>> Nl = 1
    >>> Ep = [-1.0]
    >>> epsilonp = np.array([[0, 0, 1.0]])
    >>> delta = symbols("delta")

    >>> detuning_knob = [delta]
    >>> gamma = np.array([[0.0, -1.0], [1.0, 0.0]])
    >>> omega_level = np.array([0.0, 100.0])
    >>> rm = [np.array([[0.0, 0.0], [1.0, 0.0]])*hbar_num/e_num
    ...       for p in range(3)]
    >>> xi = np.array([[[0, 1], [1, 0]]])
    >>> theta = phase_transformation(Ne, Nl, rm, xi)
    >>> sweep_steady_state = fast_sweep_steady_state(Ep, epsilonp, gamma,
    ...                                              omega_level, rm, xi,
    ...                                              theta)
    >>> deltas, rho = sweep_steady_state([[-20, 20, 11]])
    >>> n = vapour_number_density(273.15+20, "Rb")
    >>> unfolding = Unfolding(Ne, True, True, True)
    >>> chire = electric_succeptibility(0, Ep, epsilonp, rm, n,
    ...                                 rho, unfolding)
    >>> print chire
    [  4.48238603e-09 -1.12059651e-10j   5.59709041e-09 -1.74909075e-10j
       7.44587027e-09 -3.10244594e-10j   1.10969341e-08 -6.93558379e-10j
       2.14485517e-08 -2.68106896e-09j   0.00000000e+00 -5.98772067e-08j
      -2.14485517e-08 -2.68106896e-09j  -1.10969341e-08 -6.93558379e-10j
      -7.44587027e-09 -3.10244594e-10j  -5.59709041e-09 -1.74909075e-10j
      -4.48238603e-09 -1.12059651e-10j]

    """
    epsilonm = epsilonp.conjugate()
    rp = np.array([rm[i].transpose().conjugate() for i in range(3)])
    if part == 1:
        op = cartesian_dot_product(rp, epsilonm[0])
        op += cartesian_dot_product(rm, epsilonp[0])
        op = -e_num*n/epsilon_0_num/np.abs(Ep[0])*op
    elif part == -1:
        op = cartesian_dot_product(rm, epsilonp[0])
        op += - cartesian_dot_product(rp, epsilonm[0])
        op = -1j*e_num*n/epsilon_0_num/np.abs(Ep[0])*op
    elif part == 0:
        chire = electric_succeptibility(l, Ep, epsilonp, rm,
                                        n, rho, unfolding, +1)
        chiim = electric_succeptibility(l, Ep, epsilonp, rm,
                                        n, rho, unfolding, -1)
        return chire + 1j*chiim

    return np.real(observable(op, rho, unfolding))


def radiated_intensity(rho, i, j, epsilonp, rm, omega_level, xi,
                       N, D, unfolding):
    r"""Return the radiated intensity in a given direction.

    >>> from fast import State, Integer, split_hyperfine_to_magnetic
    >>> g = State("Rb", 87, 5, 1, 3/Integer(2), 0)
    >>> e = State("Rb", 87, 4, 2, 5/Integer(2), 1)
    >>> magnetic_states = split_hyperfine_to_magnetic([g, e])
    >>> omega0 = magnetic_states[0].omega
    >>> omega_level = [ei.omega - omega0 for ei in magnetic_states]
    >>> Ne = len(magnetic_states)

    >>> N = 4e6
    >>> D = 0.1
    >>> unfolding = Unfolding(Ne, True, True, True)

    >>> rho = np.zeros((Ne, Ne))
    >>> rho[0, 0] = 0.8
    >>> rho[3, 3] = 0.2
    >>> rho[3, 0] = 0.3
    >>> rho[0, 3] = 0.3
    >>> rho = unfolding(rho)

    >>> ep = np.array([1, 1j, 0])/np.sqrt(2.0)
    >>> ex = np.array([1, 0, 0])
    >>> r0 = 4.75278521538619e-11
    >>> rm = np.zeros((3, Ne, Ne), complex)
    >>> rm[0, 1, 0] = -r0
    >>> rm[0, 3, 0] = r0
    >>> rm[1, 1, 0] = -1j*r0
    >>> rm[1, 3, 0] = -1j*r0
    >>> rm[1, 2, 0] = -np.sqrt(2)*r0

    >>> xi = np.zeros((1, Ne, Ne))
    >>> xi[0, 1, 0] = 1
    >>> xi[0, 2, 0] = 1
    >>> xi[0, 3, 0] = 1
    >>> xi[0, :, :] += xi[0, :, :].transpose()

    >>> print radiated_intensity(rho, 1, 0, ex, rm,
    ...                          omega_level, xi, N, D, unfolding)
    4.60125990174e-22

    """
    def inij(i, j, ilist, jlist):
        if (i in ilist) and (j in jlist):
            return 1
        else:
            return 0

    rm = np.array(rm)
    Nl = xi.shape[0]
    Ne = xi.shape[1]
    aux = define_simplification(omega_level, xi, Nl)
    u = aux[0]
    omega_levelu = aux[2]

    ui = u(i)
    uj = u(j)
    omegaij = omega_levelu[ui] - omega_levelu[uj]

    ilist = [ii for ii in range(Ne) if u(ii) == ui]
    jlist = [jj for jj in range(Ne) if u(jj) == uj]

    rp = np.array([rm[ii].conjugate().transpose() for ii in range(3)])

    rm = np.array([[[rm[p, ii, jj]*inij(ii, jj, ilist, jlist)
                   for jj in range(Ne)]
                   for ii in range(Ne)]
                   for p in range(3)])
    rp = np.array([[[rp[p, ii, jj]*inij(jj, ii, ilist, jlist)
                   for jj in range(Ne)]
                   for ii in range(Ne)]
                   for p in range(3)])

    epsilonm = epsilonp.conjugate()

    Adag = cartesian_dot_product(rm, epsilonp)
    A = cartesian_dot_product(rp, epsilonm)

    fact = alpha_num*N*hbar_num*omegaij**3/2/np.pi/c_num**2/D**2

    Iop = fact * np.dot(Adag, A)
    intensity = observable(Iop, rho, unfolding)
    intensity = float(np.real(intensity))
    return intensity


if __name__ == "__main__":
    import doctest
    print doctest.testmod(verbose=False)
