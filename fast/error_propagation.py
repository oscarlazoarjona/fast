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

r"""
The class Measurement is an extention of floats
that have associated standard deviations or errors (x,sigmax).
Basic operations on these numbers are defined such that one
can comfortably make operations on these numbers to obtain the
corresponding errors.
"""

from math import sqrt, log


class Measurement(object):
    r"""A class for error propagation arithmetic."""

    def __init__(self, value, sigma):
        r"""A class for error propagation arithmetic."""
        self.value = float(value)
        self.sigma = sigma

    def __str__(self):
        r"""The string method for Measurement."""
        return '('+str(self.value)+', '+str(self.sigma)+')'

    def __mul__(self, other, cov=0.0):
        r"""Multiplication."""
        # Scalar multiplication
        if isinstance(other, float) or isinstance(other, int):
            return Measurement(other*self.value, abs(other)*self.sigma)
        # Measurement multiplication
        elif isinstance(other, Measurement):
            sigmaf = self.value**2 * other.sigma**2
            sigmaf += other.value**2 * self.sigma**2
            sigmaf += 2*self.value*other.value*cov
            sigmaf = sqrt(sigmaf)
            return Measurement(self.value*other.value, sigmaf)

    def __rmul__(self, other):
        r"""Reverse multiplication."""
        return self.__mul__(other)

    def __add__(self, other, cov=0.0):
        r"""Addition."""
        # Scalar addition
        if isinstance(other, float) or isinstance(other, int):
            return Measurement(other+self.value, self.sigma)
        # Measurement addition
        elif isinstance(other, Measurement):
            sigmaf = self.sigma**2 + other.sigma**2 + 2*cov
            sigmaf = sqrt(sigmaf)
            return Measurement(self.value + other.value, sigmaf)

    def __radd__(self, other):
        r"""Reverse addition."""
        return self.__add__(other)

    def __sub__(self, other, cov=0.0):
        r"""Substraction."""
        # Scalar substraction
        if isinstance(other, float) or isinstance(other, int):
            return Measurement(-other+self.value, self.sigma)
        # Measurement substraction
        elif isinstance(other, Measurement):
            sigmaf = self.sigma**2 + other.sigma**2 - 2*cov
            sigmaf = sqrt(sigmaf)
            return Measurement(self.value - other.value, sigmaf)

    def __rsub__(self, other):
        r"""Reverse substraction."""
        if isinstance(other, float) or isinstance(other, int):
            other = Measurement(other, 0.0)

        return other.__sub__(self)

    def __div__(self, other, cov=0.0):
        r"""Division."""
        # Scalar division.
        if isinstance(other, float) or isinstance(other, int):
            other = Measurement(other, 0.0)
        # Measurement division.
        sigmaf = (self.sigma/self.value)**2
        sigmaf += (other.sigma/other.value)**2 - 2*cov/(self.value*other.value)
        sigmaf = sqrt(sigmaf)
        sigmaf = sqrt((self.value/other.value)**2)*sigmaf

        return Measurement(self.value / other.value, sigmaf)

    def __rdiv__(self, other):
        r"""Reverse division."""
        if isinstance(other, float) or isinstance(other, int):
            other = Measurement(other, 0.0)

        return other.__div__(self)

    def __neg__(self):
        r"""Negative."""
        return Measurement(-self.value, self.sigma)

    def __pow__(self, other, cov=0.0):
        r"""Power."""
        # Scalar power.
        if isinstance(other, float) or isinstance(other, int):
            other = Measurement(other, 0.0)
        # Measurement power.
        sigmaf = (other.value*self.sigma/self.value)**2
        sigmaf += (log(self.value)*other.sigma)**2
        sigmaf += 2*other.value*log(self.value)*cov/self.value
        sigmaf = sqrt(sigmaf)

        return Measurement(self.value ** other.value, sigmaf)

    def __rpow__(self, other):
        r"""Reverse power."""
        if isinstance(other, float) or isinstance(other, int):
            other = Measurement(other, 0.0)

        return other.__pow__(self)


def rel(m):
    r"""Relative error."""
    return m.sigma/m.value


P = 100  # uW
err = 0.2

P = Measurement(P, P*err)
a = Measurement(4.7, 0.0)
k = 5

E02 = k*P/a**2
E0 = (k*P/a**2)**0.5
