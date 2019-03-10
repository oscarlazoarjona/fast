# -*- coding: utf-8 -*-
# ***********************************************************************
#       Copyright (C) 2014 - 2019 Oscar Gerardo Lazo Arjona             *
#               <oscar.lazoarjona@physics.ox.ac.uk>                     *
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

"""This module contains all routines related to magnetic field interaction.

>>> import numpy as np
>>> np.set_printoptions(precision=4)

"""
from atomic_structure import Atom, State
from sympy import Integer
from numpy import array, pi
from scipy.constants import physical_constants, hbar

muB = physical_constants["Bohr magneton"][0]


def lande_g_factors(element, isotope, L=None, J=None, F=None):
    r"""Return the Lande g-factors for a given atom or level.

    >>> element = "Rb"
    >>> isotope = 87
    >>> print(lande_g_factors(element, isotope))
    [ 9.9999e-01  2.0023e+00 -9.9514e-04]

    The spin-orbit g-factor for a certain J
    >>> print(lande_g_factors(element, isotope, L=0, J=1/Integer(2)))
    [0.9999936864200584 2.0023193043622 -0.0009951414 2.00231930436220]

    The nuclear-coupled g-factor for a certain F
    >>> print(lande_g_factors(element, isotope, L=0, J=1/Integer(2), F=1))
    [0.9999936864200584 2.0023193043622 -0.0009951414 2.00231930436220
     -0.501823752840550]

    """
    atom = Atom(element, isotope)
    gL = atom.gL
    gS = atom.gS
    gI = atom.gI

    res = [gL, gS, gI]
    if J is not None:
        if L is None:
            raise ValueError("A value of L must be specified.")
        S = 1/Integer(2)
        gJ = gL*(J*(J+1)-S*(S+1)+L*(L+1))/(2*J*(J+1))
        gJ += gS*(J*(J+1)+S*(S+1)-L*(L+1))/(2*J*(J+1))
        res += [gJ]
    if F is not None:
        II = atom.nuclear_spin
        if F == 0:
            gF = gJ
        else:
            gF = gJ*(F*(F+1)-II*(II+1)+J*(J+1))/(2*F*(F+1))
            gF += gI*(F*(F+1)+II*(II+1)-J*(J+1))/(2*F*(F+1))
        res += [gF]

    return array(res)


def zeeman_energies(fine_state, Bz):
    r"""Return Zeeman effect energies for a given fine state and\
    magnetic field.

    >>> ground_state = State("Rb", 87, 5, 0, 1/Integer(2))
    >>> Bz = 200.0
    >>> Bz = Bz/10000
    >>> for f_group in zeeman_energies(ground_state, Bz):
    ...     print(f_group)
    [-2.73736448508248e-24 -2.83044285506388e-24 -2.92352122504527e-24]
    [1.51284728917866e-24 1.60555650110849e-24 1.69826571303833e-24
     1.79097492496816e-24 1.88368413689800e-24]

    """
    element = fine_state.element
    isotope = fine_state.isotope
    N = fine_state.n
    L = fine_state.l
    J = fine_state.j

    energiesZeeman = []
    for i, F in enumerate(fine_state.fperm):
        gL, gS, gI, gJ, gF = lande_g_factors(element, isotope, L, J, F)
        energiesF = []
        hyperfine_level = State(element, isotope, N, L, J, F)
        for MF in range(-F, F+1):
            unperturbed_energy = hbar*hyperfine_level.omega
            energyMF = unperturbed_energy + muB*gF*MF*Bz
            energiesF += [energyMF]
        energiesZeeman += [array(energiesF)]

    return energiesZeeman


def paschen_back_energies(fine_state, Bz):
    r"""Return Paschen-Back regime energies for a given fine state and\
    magnetic field.

    >>> ground_state = State("Rb", 87, 5, 0, 1/Integer(2))
    >>> Bz = 200.0
    >>> Bz = Bz/10000
    >>> for f_group in paschen_back_energies(ground_state, Bz):
    ...     print(f_group)
    [1.51284728917866e-24 3.80485568127324e-25 -7.51876152924007e-25
     -1.88423787397534e-24]
    [-1.51229355210131e-24 -3.80300989101543e-25 7.51691573898227e-25
     1.88368413689800e-24]

    """
    element = fine_state.element
    isotope = fine_state.isotope

    N = fine_state.n
    L = fine_state.l
    J = fine_state.j
    II = Atom(element, isotope).nuclear_spin
    MJ = [-J+i for i in range(2*J+1)]
    MI = [-II+i for i in range(2*II+1)]

    Ahfs = fine_state.Ahfs
    Bhfs = fine_state.Bhfs
    gL, gS, gI, gJ = lande_g_factors(element, isotope, L, J)

    energiesPBack = []
    for mj in MJ:
        energiesMJ = []
        unperturbed_energy = hbar*State(element, isotope, N, L, J).omega
        for mi in MI:
            energyMI = unperturbed_energy
            energyMI += 2*pi*hbar*Ahfs*mi*mj
            if J != 1/Integer(2) and II != 1/Integer(2):
                num = 9*(mi*mj)**2 - 3*J*(J+1)*mi**2
                num += -3*II*(II+1)*mj**2 + II*(II+1)*J*(J+1)
                den = 4*J*(2*J-1)*II*(2*II-1)
                energyMI += 2*pi*hbar*Bhfs*num/den

            energyMI += muB*(gJ*mj+gI*mi)*Bz

            energiesMJ += [energyMI]
        energiesPBack += [energiesMJ]
    return array(energiesPBack)


if __name__ == "__main__":
    import doctest
    print(doctest.testmod(verbose=False))
