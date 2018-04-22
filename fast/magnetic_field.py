# -*- coding: utf-8 -*-
# ***********************************************************************
#       Copyright (C) 2014 - 2018 Oscar Gerardo Lazo Arjona             *
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
"""This module contains all routines related to magnetic field interaction."""
from atomic_structure import Atom
from sympy import Integer
from numpy import array


def Lande_g_factors(element, isotope, L=None, J=None, F=None):
    r"""Return the Lande g-factors for a given atom or level.

    >>> element = "Rb"
    >>> isotope = 87
    >>> print Lande_g_factors(element, isotope)
    [  9.99993686e-01   2.00231930e+00  -9.95141400e-04]

    The spin-orbit g-factor for a certain J
    >>> print Lande_g_factors(element, isotope, L=0, J=1/Integer(2))
    [0.9999936864200584 2.0023193043622 -0.0009951414 2.00231930436220]

    The nuclear-coupled g-factor for a certain F
    >>> print Lande_g_factors(element, isotope, L=0, J=1/Integer(2), F=1)
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
        gF = gJ*(F*(F+1)-II*(II+1)+J*(J+1))/(2*F*(F+1))
        gF += gI*(F*(F+1)+II*(II+1)-J*(J+1))/(2*F*(F+1))
        res += [gF]

    return array(res)


if __name__ == "__main__":
    import doctest
    print doctest.testmod(verbose=False)
