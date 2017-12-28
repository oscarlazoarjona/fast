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
"""This module contains all routines related to the optical fields."""


from scipy.constants import physical_constants
from math import pi, sqrt, cos, sin
from fast.symbolic import polarization_vector
import numpy as np

Pi = pi
# Physical constants (SI units):
#                                                           units
c = physical_constants["speed of light in vacuum"][0]      # m/s
a0 = physical_constants["Bohr radius"][0]                  # m
hbar = physical_constants["Planck constant over 2 pi"][0]  # J s
e = physical_constants["elementary charge"][0]             # C
mu0 = physical_constants["mag. constant"][0]               # N / A^2


class PlaneWave(object):
    """This class implements plane waves propagating in an arbitrary direction
    and with an arbitrary (and well-defined) polarization. It takes as input:

    INPUT:

    - ``phi`` - The azimutal angle in spherical coordinates of the wave vector.
    - ``theta`` - The polar angle in spherical coordinates of the wave vector.
    - ``alpha`` - The angle between the fast axis of a half-wave plate and an\
 incident linearly polarized beam.
    - ``beta`` - The angle between the fast axis of a quarter-wave plate and\
 an incident linearly polarized beam.

    The easiest way to understand what this means is to use the function
    draw_lasers to get a 3d view of the direction and polarization of the
    wave.

    An explicit way to understand what alpha and beta are is to see how the
    overall form of the electric field is contructed. We begin with a
    polarization vector for a plane wave propagating towards z:

    >>> from sympy import symbols, Matrix, I, exp, sin, cos, pprint
    >>> E0=symbols("E0",real=True)
    >>> t,x,y,z,omega,kx,ky,kz=symbols("t,x,y,z,omega,kx,ky,kz",real=True)
    >>> phi, theta, alpha, beta = symbols("phi, theta, alpha, beta",real=True)
    >>> ep=Matrix([cos(2*beta),I*sin(2*beta),0])

    >>> pprint(ep,use_unicode=False)
    [ cos(2*beta) ]
    [             ]
    [I*sin(2*beta)]
    [             ]
    [      0      ]

    Where beta specifies the circularity of the polarization. We define
    the following rotation matrices

    >>> R1=Matrix([[cos(2*alpha), -sin(2*alpha), 0],
    ...            [sin(2*alpha), cos(2*alpha), 0], [0, 0, 1]])
    >>> pprint(R1, use_unicode=False)
    [cos(2*alpha)  -sin(2*alpha)  0]
    [                              ]
    [sin(2*alpha)  cos(2*alpha)   0]
    [                              ]
    [     0              0        1]

    >>> R2=Matrix([[cos(theta), 0, sin(theta)],
    ...           [0, 1, 0], [-sin(theta),0,cos(theta)]])
    >>> pprint(R2, use_unicode=False)
    [cos(theta)   0  sin(theta)]
    [                          ]
    [     0       1      0     ]
    [                          ]
    [-sin(theta)  0  cos(theta)]

    >>> R3=Matrix([[cos(phi), -sin(theta), 0],
    ...           [sin(theta), cos(theta), 0], [0, 0, 1]])
    >>> pprint(R3, use_unicode=False)
    [ cos(phi)   -sin(theta)  0]
    [                          ]
    [sin(theta)  cos(theta)   0]
    [                          ]
    [    0            0       1]

    Then we create a general polarization vector applying this rotation to ep

    >>> epsilonp=R3*R2*R1*ep

    And we build the plane wave in the following way

    >>> plane_wave=E0/2 * exp(I*(kx*x+ky*y+kz*z-omega*t)) * epsilonp
    >>> plane_wave=plane_wave + plane_wave.conjugate()

    """

    def __init__(self, phi, theta, alpha, beta,
                 omega=1, E0=1, color='blue', numeric=True):
        r"""A plane wave.

        >>> PlaneWave(0, 0, 0, 0, symbolical=False)
        PlaneWave with phi=0, theta=0, alpha=0, beta=0, E0=1

        """
        self.phi = phi
        self.theta = theta
        self.alpha = alpha
        self.beta = beta
        self.E0 = E0
        self.omega = omega

        self.epsilonp = polarization_vector(phi, theta, alpha, beta, 1,
                                            numeric=numeric)
        self.epsilonm = polarization_vector(phi, theta, alpha, beta, -1,
                                            numeric=numeric)
        self.k = [cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)]
        if numeric:
            self.k = np.array([float(ki) for ki in self.k])

        if color == 'blue':
            self.color = 'b'
        elif color == 'red':
            self.color = 'r'
        elif color == 'green':
            self.color = 'g'
        Ypm1 = 0.5*sqrt(2)*((cos(2*alpha)*cos(phi)*cos(theta) -
                            sin(2*alpha)*sin(phi))*cos(2*beta) -
                            1j*(cos(phi)*cos(theta)*sin(2*alpha) +
                            cos(2*alpha)*sin(phi))*sin(2*beta)) -\
            0.5*1j*sqrt(2)*((cos(2*alpha)*cos(theta)*sin(phi) +
                            cos(phi)*sin(2*alpha))*cos(2*beta) -
                            1j*(cos(theta)*sin(2*alpha)*sin(phi) -
                            cos(2*alpha)*cos(phi))*sin(2*beta))

        Yp0 = -cos(2*alpha)*cos(2*beta)*sin(theta) +\
            1j*sin(2*alpha)*sin(2*beta)*sin(theta)

        Yp1 = -0.5*sqrt(2)*((cos(2*alpha)*cos(phi)*cos(theta) -
                            sin(2*alpha)*sin(phi))*cos(2*beta) -
                            1j*(cos(phi)*cos(theta)*sin(2*alpha) +
                            cos(2*alpha)*sin(phi))*sin(2*beta)) -\
            0.5*1j*sqrt(2)*((cos(2*alpha)*cos(theta)*sin(phi) +
                            cos(phi)*sin(2*alpha))*cos(2*beta) -
                            1j*(cos(theta)*sin(2*alpha)*sin(phi) -
                            cos(2*alpha)*cos(phi))*sin(2*beta))

        Ymm1 = 0.5*sqrt(2)*((cos(2*alpha)*cos(phi)*cos(theta) -
                            sin(2*alpha)*sin(phi))*cos(2*beta) +
                            1j*(cos(phi)*cos(theta)*sin(2*alpha) +
                            cos(2*alpha)*sin(phi))*sin(2*beta)) -\
            0.5*1j*sqrt(2)*((cos(2*alpha)*cos(theta)*sin(phi) +
                            cos(phi)*sin(2*alpha))*cos(2*beta) +
                            1j*(cos(theta)*sin(2*alpha)*sin(phi) -
                            cos(2*alpha)*cos(phi))*sin(2*beta))

        Ym0 = -cos(2*alpha)*cos(2*beta)*sin(theta) -\
            1j*sin(2*alpha)*sin(2*beta)*sin(theta)

        Ym1 = -0.5*sqrt(2)*((cos(2*alpha)*cos(phi)*cos(theta) -
                            sin(2*alpha)*sin(phi))*cos(2*beta) +
                            1j*(cos(phi)*cos(theta)*sin(2*alpha) +
                            cos(2*alpha)*sin(phi))*sin(2*beta)) -\
            0.5*1j*sqrt(2)*((cos(2*alpha)*cos(theta)*sin(phi) +
                            cos(phi)*sin(2*alpha))*cos(2*beta) +
                            1j*(cos(theta)*sin(2*alpha)*sin(phi) -
                            cos(2*alpha)*cos(phi))*sin(2*beta))

        # We set the components of the polarization vectors epsilon^(+-)
        # in the helicity basis.
        self.Yp = [Ypm1, Yp0, Yp1]
        self.Ym = [Ymm1, Ym0, Ym1]

        for i in range(3):

            a = self.Yp[i]; rea = a.real; ima = a.imag
            b = self.Ym[i]; reb = b.real; imb = b.imag

            if abs(rea) < 1e-15: rea = 0
            if abs(ima) < 1e-15: ima = 0
            if abs(reb) < 1e-15: reb = 0
            if abs(imb) < 1e-15: imb = 0

            self.Yp[i] = rea+1j*ima
            self.Ym[i] = reb+1j*imb

    def __str__(self):
        r"""The str routine for plane waves.

        >>> w = PlaneWave(0, 0, 0, 0, symbolical=False)
        >>> w.__str__()
        'PlaneWave with phi=0, theta=0, alpha=0, beta=0, E0=1'

        """
        s = 'PlaneWave with phi='+str(self.phi)+', theta='+str(self.theta)
        s += ', alpha='+str(self.alpha)+', beta='+str(self.beta)
        s += ', E0='+str(self.E0)
        return s

    def __repr__(self):
        r"""The repr routine for plane waves.

        >>> w = PlaneWave(0, 0, 0, 0, symbolical=False)
        >>> w.__repr__()
        'PlaneWave with phi=0, theta=0, alpha=0, beta=0, E0=1'

        """
        return self.__str__()


class MotField(object):
    r"""The optical field of a MOT scheme."""

    def __init__(self, relative_intensities, parity=1, color='blue'):
        r"""A MOT field.

        >>> MotField([1, 1, 1, 1, 1, 1])
        MotField with relative intensities [1, 1, 1, 1, 1, 1]

        """
        lx = PlaneWave(Pi, Pi/2, 0, parity*Pi/8, symbolical=False, color=color)
        lx_r = PlaneWave(0, Pi/2, 0, parity*Pi/8, symbolical=False,
                         color=color)

        ly = PlaneWave(-Pi/2, Pi/2, 0, parity*Pi/8, symbolical=False,
                       color=color)
        ly_r = PlaneWave(Pi/2, Pi/2, 0, parity*Pi/8, symbolical=False,
                         color=color)

        lz = PlaneWave(0, Pi, 0, -parity*Pi/8, symbolical=False, color=color)
        lz_r = PlaneWave(0, 0, 0, -parity*Pi/8, symbolical=False, color=color)

        self.lx = lx; self.ly = ly; self.lz = lz
        self.lx_r = lx_r; self.ly_r = ly_r; self.lz_r = lz_r
        self.beams = [lx, ly, lz, lx_r, ly_r, lz_r]

        Yp = [0 for i in range(3)]; Ym = [0 for i in range(3)]
        for i in range(3):
            Yp[i] += relative_intensities[0]*lx.Yp[i]
            Yp[i] += relative_intensities[1]*ly.Yp[i]
            Yp[i] += relative_intensities[2]*lz.Yp[i]
            Yp[i] += relative_intensities[3]*lx_r.Yp[i]
            Yp[i] += relative_intensities[4]*ly_r.Yp[i]
            Yp[i] += relative_intensities[5]*lz_r.Yp[i]

            Ym[i] += relative_intensities[0]*lx.Ym[i]
            Ym[i] += relative_intensities[1]*ly.Ym[i]
            Ym[i] += relative_intensities[2]*lz.Ym[i]
            Ym[i] += relative_intensities[3]*lx_r.Ym[i]
            Ym[i] += relative_intensities[4]*ly_r.Ym[i]
            Ym[i] += relative_intensities[5]*lz_r.Ym[i]

        self.relative_intensities = relative_intensities
        self.Yp = Yp
        self.Ym = Ym

    def __str__(self):
        r"""The string routine for MotField

        >>> lasers = MotField([1, 1, 1, 1, 1, 1])
        >>> lasers.__str__()
        'MotField with relative intensities [1, 1, 1, 1, 1, 1]'

        """
        s = "MotField with relative intensities "
        s += str(self.relative_intensities)
        return s

    def __repr__(self):
        r"""The repr routine for MotField

        >>> lasers = MotField([1, 1, 1, 1, 1, 1])
        >>> lasers.__repr__()
        'MotField with relative intensities [1, 1, 1, 1, 1, 1]'

        """
        return self.__str__()

    def plot(self, **kwds):
        r"""The plotting function for MOT fields."""
        plo = self.lx.plot(dist_to_center=3, **kwds)
        plo += self.ly.plot(dist_to_center=3, **kwds)
        plo += self.lz.plot(dist_to_center=3, **kwds)
        plo += self.lx_r.plot(dist_to_center=3, **kwds)
        plo += self.ly_r.plot(dist_to_center=3, **kwds)
        plo += self.lz_r.plot(dist_to_center=3, **kwds)
        return plo


def electric_field_amplitude_gaussian(P, sigmax, sigmay=None, Omega=1.0e6,
                                      units="ad-hoc"):
    """Return the amplitude of the electric field for a Gaussian beam.

    This the amplitude at the center of a laser beam of power P (in Watts) and\
 a Gaussian intensity distribution of standard deviations sigmax, sigmay\
 (in meters). The value of E0 is given in rescaled units according to the\
 frequency scale  Omega (in Hertz) understood as absolute frequency\
 (as opposed to angular frequency).

    >>> print electric_field_amplitude_gaussian(0.001, 0.001)
    19.6861467587

    """
    e0 = hbar*Omega/(e*a0)  # This is the electric field scale.

    if sigmay is None: sigmay = sigmax
    E0 = sqrt((c*mu0*P)/(2*Pi))/sqrt(sigmax*sigmay)
    if units == "ad-hoc":
        E0 = E0/e0
    return E0


def electric_field_amplitude_top(P, a, Omega=1e6, units="ad-hoc"):
    """Return the amplitude of the electric field for a top hat beam.

    This is the amplitude of a laser beam of power P (in Watts) and a top-hat\
 intensity distribution of radius a (in meters). The value of E0 is given in\
 rescaled units according to the frequency scale  Omega (in Hertz)\
 understood as absolute frequency (as opposed to angular frequency).

    >>> print electric_field_amplitude_top(0.001, 0.001)
    27.8404157371

    """
    e0 = hbar*Omega/(e*a0)  # This is the electric field scale.
    E0 = sqrt((c*mu0*P)/(Pi*a**2))
    if units == "ad-hoc":
        E0 = E0/e0
    return E0


def electric_field_amplitude_intensity(s0, Isat=16.6889462814,
                                       Omega=1e6, units="ad-hoc"):
    """Return the amplitude of the electric field for saturation parameter.

    This is at a given saturation parameter s0=I/Isat, where I0 is by default \
Isat=16.6889462814 m/m^2 is the saturation intensity of the D2 line of \
rubidium for circularly polarized light. Optionally, a frequency scale \
`Omega` can be provided.

    >>> print electric_field_amplitude_intensity(1.0, units="ad-hoc")
    9.0152984553
    >>> print electric_field_amplitude_intensity(1.0, Omega=1.0, units="SI")
    112.135917207
    >>> print electric_field_amplitude_intensity(1.0, units="SI")
    0.000112135917207

    """
    E0_sat = sqrt(2*mu0*c*Isat)/Omega

    if units == "ad-hoc":
        e0 = hbar/(e*a0)  # This is the electric field scale.
        E0_sat = E0_sat/e0

    return E0_sat*sqrt(s0)


if __name__ == "__main__":
    import doctest
    print doctest.testmod(verbose=False)
