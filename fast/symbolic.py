# -*- coding: utf-8 -*-

# ************************************************************************
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
r"""This module contains all the routines to produce symbolic equations.

A simple example:

>>> define_density_matrix(2)
Matrix([
[rho11, rho12],
[rho21, rho22]])

>>> import numpy as np
>>> np.set_printoptions(precision=4)

References
~~~~~~~~~~
.. [Edmonds74] A. R. Edmonds. Angular momentum in quantum mechanics.
  Investigations in physics, 4.; Investigations in physics, no. 4.
  Princeton, N.J., Princeton University Press, 1957.

"""

from sympy import Symbol, Matrix, symbols
from sympy import I, conjugate, diff
from sympy import sin, cos, sqrt, exp
from sympy import KroneckerDelta
from sympy import Function, Derivative
from sympy import re, im, zeros, IndexedBase
from fast.misc import IJ, find_phase_transformation
from numpy import array as nparray
from numpy import sqrt as npsqrt
from sympy import Expr
from sympy.physics.quantum.qexpr import QExpr
from sympy.core import Mul
from time import time
# from sympy.core.compatibility import u


def define_symbol(name, open_brace, comma, i, j,
                  close_brace, variables, **kwds):
    r"""Define a nice symbol with matrix indices.

    >>> name = "rho"
    >>> from sympy import symbols
    >>> t, x, y, z = symbols("t, x, y, z", positive=True)
    >>> variables = [t, x, y, z]
    >>> open_brace = ""
    >>> comma = ""
    >>> close_brace = ""
    >>> i = 0
    >>> j = 1
    >>> f = define_symbol(name, open_brace, comma, i, j, close_brace,
    ...                   variables, positive=True)
    >>> print f
    rho12(t, x, y, z)

    """
    if variables is None:
        return Symbol(name+open_brace+str(i+1)+comma+str(j+1) +
                      close_brace, **kwds)
    else:
        return Function(name+open_brace+str(i+1)+comma+str(j+1) +
                        close_brace, **kwds)(*variables)


def define_density_matrix(Ne, explicitly_hermitian=False, normalized=False,
                          name="rho", variables=None):
    """Return a symbolic density matrix.

    INPUT:

    - ``Ne`` - The number of atomic states.
    - ``explicitly_hermitian`` - (boolean) whether to make\
 :math:`\\rho_{ij}=\\bar{\\rho}_{ij}` for i<j.
    - ``normalized`` - (boolean): Whether to make\
 :math:`\\rho_{11}=1-\\sum_{i>1} \\rho_{ii}`

    A very simple example:

    >>> define_density_matrix(2)
    Matrix([
    [rho11, rho12],
    [rho21, rho22]])

    The density matrix can be made explicitly hermitian

    >>> define_density_matrix(2, explicitly_hermitian=True)
    Matrix([
    [rho11, conjugate(rho21)],
    [rho21,            rho22]])

    or normalized

    >>> define_density_matrix(2, normalized=True)
    Matrix([
    [-rho22 + 1, rho12],
    [     rho21, rho22]])

    or it can be made an explicit function of given variables

    >>> from sympy import symbols
    >>> t, z = symbols("t, z", positive=True)
    >>> define_density_matrix(2, variables=[t, z])
    Matrix([
    [rho11(t, z), rho12(t, z)],
    [rho21(t, z), rho22(t, z)]])

    """
    if Ne > 9:
        comma = ","
        name = "\\" + name
        open_brace = "_{"
        close_brace = "}"
    else:
        comma = ""
        open_brace = ""
        close_brace = ""

    rho = []
    for i in range(Ne):
        row_rho = []
        for j in range(Ne):
            if i == j:
                row_rho += [define_symbol(name, open_brace, comma, i, j,
                                          close_brace, variables,
                                          real=True)]
            elif i > j:
                row_rho += [define_symbol(name, open_brace, comma, i, j,
                                          close_brace, variables)]
            else:
                if explicitly_hermitian:
                    row_rho += [conjugate(define_symbol(name, open_brace,
                                                        comma, j, i,
                                                        close_brace,
                                                        variables))]
                else:
                    row_rho += [define_symbol(name, open_brace, comma, i, j,
                                              close_brace, variables)]
        rho += [row_rho]

    if normalized:
        rho11 = 1-sum([rho[i][i] for i in range(1, Ne)])
        rho[0][0] = rho11

    rho = Matrix(rho)
    return rho


def define_laser_variables(Nl, real_amplitudes=False, variables=None):
    r"""Return the amplitudes and frequencies of Nl fields.

    >>> E0, omega_laser = define_laser_variables(2)
    >>> E0, omega_laser
    ([E_{01}, E_{02}], [varpi_1, varpi_2])

    The amplitudes are complex by default:

    >>> conjugate(E0[0])
    conjugate(E_{01})

    But they can optionally be made real:

    >>> E0, omega_laser = define_laser_variables(2, real_amplitudes=True)
    >>> conjugate(E0[0])
    E_{01}

    They can also be made explicit functions of given variables:

    >>> from sympy import symbols
    >>> t, z = symbols("t, z", real=True)
    >>> E0, omega_laser = define_laser_variables(2, variables=[t, z])
    >>> E0
    [E_{01}(t, z), E_{02}(t, z)]

    """
    if variables is None:
        E0 = [Symbol(r"E_{0"+str(l+1)+"}", real=real_amplitudes)
              for l in range(Nl)]
    else:
        E0 = [Function(r"E_{0"+str(l+1)+"}", real=real_amplitudes)(*variables)
              for l in range(Nl)]

    omega_laser = [Symbol(r"varpi_"+str(l+1), positive=True)
                   for l in range(Nl)]
    return E0, omega_laser


def polarization_vector(phi, theta, alpha, beta, p,
                        numeric=False, abstract=False):
    """This function returns a unitary vector describing the polarization
    of plane waves.:

    INPUT:

    -  ``phi`` - The spherical coordinates azimuthal angle of the wave vector\
 k.
    -  ``theta`` - The spherical coordinates polar angle of the wave vector k.
    -  ``alpha`` - The rotation of a half-wave plate.
    -  ``beta`` - The rotation of a quarter-wave plate.
    -  ``p`` - either 1 or -1 to indicate whether to return epsilon^(+) or\
 epsilon^(-) respectively.

    If alpha and beta are zero, the result will be linearly polarized light
    along some fast axis. alpha and beta are measured from that fast axis.

    Propagation towards y, linear polarization (for pi transitions):

    >>> from sympy import pi
    >>> polarization_vector(phi=pi/2, theta=pi/2, alpha=pi/2, beta= 0,p=1)
    Matrix([
    [0],
    [0],
    [1]])

    Propagation towards +z, circular polarization (for sigma + transitions):

    >>> polarization_vector(phi=0, theta= 0, alpha=pi/2, beta= pi/8,p=1)
    Matrix([
    [  -sqrt(2)/2],
    [-sqrt(2)*I/2],
    [           0]])

    Propagation towards -z, circular polarization for sigma + transitions:

    >>> polarization_vector(phi=0, theta=pi, alpha=   0, beta=-pi/8,p=1)
    Matrix([
    [  -sqrt(2)/2],
    [-sqrt(2)*I/2],
    [           0]])

    Components + and - are complex conjugates of each other

    >>> from sympy import symbols
    >>> phi, theta, alpha, beta = symbols("phi theta alpha beta", real=True)
    >>> ep = polarization_vector(phi,theta,alpha,beta, 1)
    >>> em = polarization_vector(phi,theta,alpha,beta,-1)
    >>> ep-em.conjugate()
    Matrix([
    [0],
    [0],
    [0]])

    We can also define abstract polarization vectors without explicit \
    components

    >>> polarization_vector(0, 0, 0, 0, 1, abstract=True)
    epsilonp
    >>> polarization_vector(0, 0, 0, 0, -1, abstract=True)
    epsilonm

    """
    if abstract:
        Nl = symbols("N_l", integer=True)
        if p == 1:
            epsilon = Vector3D(IndexedBase("epsilonp", shape=(Nl,)))
        else:
            epsilon = Vector3D(IndexedBase("epsilonm", shape=(Nl,)))
        return epsilon

    epsilon = Matrix([cos(2*beta), p*I*sin(2*beta), 0])

    R1 = Matrix([[cos(2*alpha), -sin(2*alpha), 0],
                [sin(2*alpha), cos(2*alpha), 0],
                [0, 0, 1]])

    R2 = Matrix([[cos(theta), 0, sin(theta)],
                 [0, 1, 0],
                 [-sin(theta), 0, cos(theta)]])

    R3 = Matrix([[cos(phi), -sin(phi), 0],
                 [sin(phi), cos(phi), 0],
                 [0, 0, 1]])

    epsilon = R3*R2*R1*epsilon
    if numeric:
        epsilon = nparray([complex(epsilon[i]) for i in range(3)])
    return epsilon


def cartesian_to_helicity(vector, numeric=False):
    r"""This function takes vectors from the cartesian basis to the helicity basis.
    For instance, we can check what are the vectors of the helicity basis.

    >>> from sympy import pi
    >>> em=polarization_vector(phi=0, theta= 0, alpha=0, beta=-pi/8,p= 1)
    >>> em
    Matrix([
    [   sqrt(2)/2],
    [-sqrt(2)*I/2],
    [           0]])
    >>> cartesian_to_helicity(em)
    Matrix([
    [ 0],
    [ 0],
    [-1]])

    >>> e0=polarization_vector(phi=pi/2, theta=pi/2, alpha=pi/2, beta=0,p=1)
    >>> e0
    Matrix([
    [0],
    [0],
    [1]])
    >>> cartesian_to_helicity(e0)
    Matrix([
    [0],
    [1],
    [0]])

    >>> ep=polarization_vector(phi=0, theta= 0, alpha=pi/2, beta= pi/8,p= 1)
    >>> ep
    Matrix([
    [  -sqrt(2)/2],
    [-sqrt(2)*I/2],
    [           0]])
    >>> cartesian_to_helicity(ep)
    Matrix([
    [-1],
    [ 0],
    [ 0]])

    Note that vectors in the helicity basis are built in a weird way by
    convention:

    .. math::

        \vec{a} = -a_{+1}\vec{e}_{-1} +a_0\vec{e}_0 -a_{-1}\vec{e}_{+1}

    >>> from sympy import symbols
    >>> am,a0,ap = symbols("am a0 ap")
    >>> a=-ap*em +a0*e0 -am*ep
    >>> a
    Matrix([
    [    sqrt(2)*am/2 - sqrt(2)*ap/2],
    [sqrt(2)*I*am/2 + sqrt(2)*I*ap/2],
    [                             a0]])
    >>> cartesian_to_helicity(a).expand()
    Matrix([
    [am],
    [a0],
    [ap]])

    We can also convert a numeric array

    >>> r =[[[0.0, 1.0],
    ...      [1.0, 0.0]],
    ...     [[0.0, -1j],
    ...      [ 1j, 0.0]],
    ...     [[1.0, 0.0],
    ...      [0.0,-1.0]]]

    >>> cartesian_to_helicity(r, numeric=True)
    array([[[ 0.    +0.j,  0.    +0.j],
            [ 1.4142+0.j,  0.    +0.j]],
    <BLANKLINE>
           [[ 1.    +0.j,  0.    +0.j],
            [ 0.    +0.j, -1.    +0.j]],
    <BLANKLINE>
           [[-0.    +0.j, -1.4142+0.j],
            [-0.    +0.j, -0.    +0.j]]])

    """
    if numeric:
        vector = list(vector)
        vector[0] = nparray(vector[0])
        vector[1] = nparray(vector[1])
        vector[2] = nparray(vector[2])

        v = [(vector[0]-1j*vector[1])/npsqrt(2),
             vector[2],
             -(vector[0]+1j*vector[1])/npsqrt(2)]
        v = nparray(v)
    else:
        v = [(vector[0]-I*vector[1])/sqrt(2),
             vector[2],
             -(vector[0]+I*vector[1])/sqrt(2)]

    if type(vector[0]) in [type(Matrix([1, 0])), type(nparray([1, 0]))]:
        return v
    else:
        return Matrix(v)


def helicity_to_cartesian(vector, numeric=False):
    r"""Transform a vector in the helicity basis to the cartesian basis.

    >>> sigmam = [1, 0, 0]
    >>> helicity_to_cartesian(sigmam)
    Matrix([
    [  sqrt(2)/2],
    [sqrt(2)*I/2],
    [          0]])

    The input vector can be a list of matrices

    >>> r = define_r_components(2, helicity=True)
    >>> r[0][0, 1] = 0
    >>> r[1][0, 1] = 0
    >>> r[2][0, 1] = 0
    >>> r
    [Matrix([
    [        0, 0],
    [r_{-1;21}, 0]]), Matrix([
    [       0, 0],
    [r_{0;21}, 0]]), Matrix([
    [        0, 0],
    [r_{+1;21}, 0]])]

    >>> helicity_to_cartesian(r)
    [Matrix([
    [                                 0, 0],
    [sqrt(2)*(-r_{+1;21} + r_{-1;21})/2, 0]]), Matrix([
    [                                  0, 0],
    [sqrt(2)*I*(r_{+1;21} + r_{-1;21})/2, 0]]), Matrix([
    [       0, 0],
    [r_{0;21}, 0]])]

    We can also convert a numeric array

    >>> r =[[[0.0        ,        0.0 ],
    ...      [npsqrt(2.0),        0.0 ]],
    ...     [[1.0        ,        0.0 ],
    ...      [0.0        ,       -1.0 ]],
    ...     [[0.0        ,-npsqrt(2.0)],
    ...      [0.0        ,        0.0 ]]]

    >>> helicity_to_cartesian(r, numeric=True)
    array([[[ 0.+0.j,  1.+0.j],
            [ 1.+0.j,  0.+0.j]],
    <BLANKLINE>
           [[ 0.+0.j, -0.-1.j],
            [ 0.+1.j,  0.+0.j]],
    <BLANKLINE>
           [[ 1.+0.j,  0.+0.j],
            [ 0.+0.j, -1.+0.j]]])

    """
    if numeric:
        vector = list(vector)
        vector[0] = nparray(vector[0])
        vector[1] = nparray(vector[1])
        vector[2] = nparray(vector[2])

        v = [(vector[0]-vector[2])/npsqrt(2),
             1j*(vector[0]+vector[2])/npsqrt(2),
             vector[1]]
        v = nparray(v)
    else:
        v = [(vector[0]-vector[2])/sqrt(2),
             I*(vector[0]+vector[2])/sqrt(2),
             vector[1]]

    if type(vector[0]) in [type(Matrix([1, 0])), type(nparray([1, 0]))]:
        return v
    else:
        return Matrix(v)


def helicity_dot_product(v1, v2):
    r"""Calculate the dot product of two vectors in the helicity basis.

    >>> from sympy import symbols
    >>> u = Matrix(symbols("u_{-1}, u_0, u_{+1}"))
    >>> v = Matrix(symbols("v_{-1}, v_0, v_{+1}"))
    >>> helicity_dot_product(u, v)
    u_0*v_0 - u_{+1}*v_{-1} - u_{-1}*v_{+1}

    We can check that this is in deed the same thing as the usual cartesian
    dot product:

    >>> u = Matrix(symbols("u_x, u_y, u_z"))
    >>> v = Matrix(symbols("v_x, v_y, v_z"))
    >>> u_helicity = cartesian_to_helicity(u)
    >>> v_helicity = cartesian_to_helicity(v)
    >>> helicity_dot_product(u_helicity, v_helicity).expand()
    u_x*v_x + u_y*v_y + u_z*v_z

    The inputs vectors can be a list of matrices representing operators.

    >>> rp = define_r_components(2, helicity=True, p=1,
    ...                          explicitly_hermitian=True)
    >>> rm = define_r_components(2, helicity=True, p=-1,
    ...                          explicitly_hermitian=True)
    >>> from sympy import pi
    >>> em = polarization_vector(0, 0, 0, pi/8, -1)
    >>> ep = polarization_vector(0, 0, 0, pi/8, 1)
    >>> ep = cartesian_to_helicity(ep)
    >>> em = cartesian_to_helicity(em)
    >>> H = helicity_dot_product(ep, rm) + helicity_dot_product(em, rp)
    >>> H
    Matrix([
    [         0, -r_{+1;21}],
    [-r_{+1;21},          0]])

    """
    return -v1[2]*v2[0] + v1[1]*v2[1] - v1[0]*v2[2]


def cartesian_dot_product(v1, v2):
    r"""Calculate the dot product of two vectors in the cartesian basis.

    >>> from sympy import symbols
    >>> u = Matrix(symbols("u_x, u_y, u_z"))
    >>> v = Matrix(symbols("v_x, v_y, v_z"))
    >>> cartesian_dot_product(u, v)
    u_x*v_x + u_y*v_y + u_z*v_z
    """
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]


def define_r_components(Ne, xi=None, explicitly_hermitian=False,
                        helicity=False, real=True, p=None):
    r"""Define the components of the position operators.

    In general, these are representations of the position operators x, y, z

    >>> define_r_components(2)
    [Matrix([
    [     0, x_{12}],
    [x_{21},      0]]), Matrix([
    [     0, y_{12}],
    [y_{21},      0]]), Matrix([
    [     0, z_{12}],
    [z_{21},      0]])]

    We can make these operators explicitly hermitian

    >>> define_r_components(2, explicitly_hermitian=True)
    [Matrix([
    [     0, x_{21}],
    [x_{21},      0]]), Matrix([
    [     0, y_{21}],
    [y_{21},      0]]), Matrix([
    [     0, z_{21}],
    [z_{21},      0]])]

    Make them real

    >>> r = define_r_components(2, real=True, explicitly_hermitian=True)
    >>> print [r[p]-r[p].transpose() for p in range(3)]
    [Matrix([
    [0, 0],
    [0, 0]]), Matrix([
    [0, 0],
    [0, 0]]), Matrix([
    [0, 0],
    [0, 0]])]

    We can get the components of the operator in the helicity basis

    >>> define_r_components(2, helicity=True)
    [Matrix([
    [        0, r_{-1;12}],
    [r_{-1;21},         0]]), Matrix([
    [       0, r_{0;12}],
    [r_{0;21},        0]]), Matrix([
    [        0, r_{+1;12}],
    [r_{+1;21},         0]])]

    And combinations thereof. For instance, let us check that the components
    in the helicity basis produce hermitian operators in the cartesian basis.

    >>> r_helicity = define_r_components(2, helicity=True,
    ...                                  explicitly_hermitian=True)

    [Matrix([
    [        0, -r_{+1;21}],
    [r_{-1;21},          0]]), Matrix([
    [       0, r_{0;21}],
    [r_{0;21},        0]]), Matrix([
    [        0, -r_{-1;21}],
    [r_{+1;21},          0]])]

    >>> r_cartesian = helicity_to_cartesian(r_helicity)
    >>> r_cartesian[0]
    Matrix([
    [                                 0, sqrt(2)*(-r_{+1;21} + r_{-1;21})/2],
    [sqrt(2)*(-r_{+1;21} + r_{-1;21})/2,                                  0]])


    >>> [(r_cartesian[p]-r_cartesian[p].adjoint()).expand() for p in range(3)]
    [Matrix([
    [0, 0],
    [0, 0]]), Matrix([
    [0, 0],
    [0, 0]]), Matrix([
    [0, 0],
    [0, 0]])]

    """
    frequency_sign = p
    if Ne > 9: comma = ","
    else: comma = ""

    if helicity:
        names = ["r_{-1;", "r_{0;", "r_{+1;"]
    else:
        names = ["x", "y", "z"]

    r = []
    if helicity:
        for p in range(3):
            r_comp = []
            for i in range(Ne):
                r_row = []
                for j in range(Ne):
                    if i == j:
                        r_row += [0]
                    elif i > j:
                        r_row += [Symbol(names[p]+str(i+1)+comma+str(j+1)+"}",
                                  real=real)]
                    elif explicitly_hermitian:
                        sign = int((-1)**(p-1))
                        r_row += [sign*conjugate(Symbol(names[2-p]+str(j+1) +
                                                        comma+str(i+1)+"}",
                                                        real=real))]
                    else:
                        r_row += [Symbol(names[p]+str(i+1)+comma+str(j+1)+"}",
                                         real=real)]
                r_comp += [r_row]
            r_comp = Matrix(r_comp)
            r += [r_comp]

    else:
        for p in range(3):
            r_comp = []
            for i in range(Ne):
                r_row = []
                for j in range(Ne):
                    if i == j:
                        r_row += [0]
                    elif i > j:
                        r_row += [Symbol(names[p]+r"_{"+str(i+1) +
                                         comma+str(j+1)+"}", real=real)]
                    elif explicitly_hermitian:
                        r_row += [conjugate(Symbol(names[p]+r"_{"+str(j+1) +
                                                   comma+str(i+1)+"}",
                                                   real=real))]
                    else:
                        r_row += [Symbol(names[p]+r"_{"+str(i+1) +
                                         comma+str(j+1)+"}", real=real)]
                r_comp += [r_row]
            r_comp = Matrix(r_comp)
            r += [r_comp]

    # We select only the upper diagonal or lower diagonal components according
    # to the sign r^(+) or r^(-) provided.
    if frequency_sign == 1:
        r = [Matrix([[r[p][i, j]*delta_lesser(i, j)
                     for j in range(Ne)] for i in range(Ne)])
             for p in range(3)]
    elif frequency_sign == -1:
        r = [Matrix([[r[p][i, j]*delta_greater(i, j)
                     for j in range(Ne)] for i in range(Ne)])
             for p in range(3)]

    if xi is not None:
        Nl = len(xi)
        for p in range(3):
            for i in range(Ne):
                for j in range(Ne):
                    zero = True
                    for l in range(Nl):
                        if xi[l][i, j] != 0:
                            zero = False
                    if zero:
                        r[p][i, j] = 0
    return r


def vector_element(r, i, j):
    r"""Extract an matrix element of a vector operator.

    >>> r = define_r_components(2)
    >>> vector_element(r, 1, 0)
    Matrix([
    [x_{21}],
    [y_{21}],
    [z_{21}]])

    """
    return Matrix([r[p][i, j] for p in range(3)])


def define_frequencies(Ne, explicitly_antisymmetric=False):
    u"""Define all frequencies omega_level, omega, gamma.

    >>> from sympy import pprint
    >>> pprint(define_frequencies(2), use_unicode=True)
    ⎛          ⎡ 0   ω₁₂⎤  ⎡ 0   γ₁₂⎤⎞
    ⎜[ω₁, ω₂], ⎢        ⎥, ⎢        ⎥⎟
    ⎝          ⎣ω₂₁   0 ⎦  ⎣γ₂₁   0 ⎦⎠

    We can make these matrices explicitly antisymmetric.

    >>> pprint(define_frequencies(2, explicitly_antisymmetric=True),
    ...                           use_unicode=True)
    ⎛          ⎡ 0   -ω₂₁⎤  ⎡ 0   -γ₂₁⎤⎞
    ⎜[ω₁, ω₂], ⎢         ⎥, ⎢         ⎥⎟
    ⎝          ⎣ω₂₁   0  ⎦  ⎣γ₂₁   0  ⎦⎠

    """
    omega_level = [Symbol('omega_'+str(i+1), real=True) for i in range(Ne)]

    if Ne > 9:
        opening = "\\"
        comma = ","
        open_brace = "{"
        close_brace = "}"
    else:
        opening = r""
        comma = ""
        open_brace = ""
        close_brace = ""

    omega = []; gamma = []
    for i in range(Ne):
        row_omega = []; row_gamma = []
        for j in range(Ne):
            if i == j:
                om = 0; ga = 0
            elif i > j:
                om = Symbol(opening+r"omega_" +
                            open_brace+str(i+1)+comma+str(j+1) +
                            close_brace, real=True)
                ga = Symbol(opening+r"gamma_" +
                            open_brace+str(i+1)+comma+str(j+1) +
                            close_brace, real=True)
            elif explicitly_antisymmetric:
                om = -Symbol(opening+r"omega_" +
                             open_brace+str(j+1)+comma+str(i+1) +
                             close_brace, real=True)
                ga = -Symbol(opening+r"gamma_" +
                             open_brace+str(j+1)+comma+str(i+1) +
                             close_brace, real=True)
            else:
                om = Symbol(opening+r"omega_" +
                            open_brace+str(i+1)+comma+str(j+1) +
                            close_brace, real=True)
                ga = Symbol(opening+r"gamma_" +
                            open_brace+str(i+1)+comma+str(j+1) +
                            close_brace, real=True)

            row_omega += [om]
            row_gamma += [ga]

        omega += [row_omega]
        gamma += [row_gamma]

    omega = Matrix(omega)
    gamma = Matrix(gamma)

    return omega_level, omega, gamma


def delta_greater(i, j):
    r"""A symbol that 1 if i > j and zero otherwise.

    >>> delta_greater(2, 1)
    1

    >>> delta_greater(1, 2)
    0
    """
    if i > j: return 1
    else: return 0


def delta_lesser(i, j):
    r"""A symbol that 1 if i < j and zero otherwise.

    >>> delta_lesser(2, 1)
    0

    >>> delta_lesser(1, 2)
    1
    """
    if i < j: return 1
    else: return 0


def bra(i, Ne):
    r"""This function returns the transpose of the i-th element of the
    canonical basis of a Hilbert space of dimension Ne (in the form of a
    row vector).

    >>> bra(2,4)
    Matrix([[0, 1, 0, 0]])

    This will return an error if i is not in [1 .. Ne]:

    >>> bra(5,3)
    Traceback (most recent call last):
    ...
    ValueError: i must be in [1 .. Ne].

    """
    if i not in range(1, Ne+1):
        raise ValueError("i must be in [1 .. Ne].")
    return Matrix([KroneckerDelta(i-1, j) for j in range(Ne)]).transpose()


def ket(i, Ne):
    r"""This function returns the i-th element of the canonical basis
    of a Hilbert space of dimension Ne (in the form of a column vector).

    >>> ket(2,4)
    Matrix([
    [0],
    [1],
    [0],
    [0]])

    This will return an error if i is not in [1 .. Ne]:

    >>> ket(5,3)
    Traceback (most recent call last):
    ...
    ValueError: i must be in [1 .. Ne].

    """
    if i not in range(1, Ne+1):
        raise ValueError("i must be in [1 .. Ne].")
    return Matrix([KroneckerDelta(i-1, j) for j in range(Ne)])


def ketbra(i, j, Ne):
    """This function returns the outer product :math:`|i><j|` where\
 :math:`|i>` and :math:`|j>` are elements of the canonical basis of an\
 Ne-dimensional Hilbert space (in matrix form).

    >>> ketbra(2, 3, 3)
    Matrix([
    [0, 0, 0],
    [0, 0, 1],
    [0, 0, 0]])

    """
    return ket(i, Ne)*bra(j, Ne)


def sigma_operator_indices(A):
    r"""If A is an outer-product type operator |a><b| return a, b.

    >>> sig = ket(2, 3)*bra(1, 3)
    >>> sigma_operator_indices(sig)
    (1, 0)
    >>> sigma_operator_indices(sig+sig.adjoint())
    (None, None)

    """
    Ne = A.shape[0]
    band = True
    if sum(A) != 1:
        band = False

    a = None; b = None
    for i in range(Ne):
        for j in range(Ne):
            if A[i, j] == 1:
                a = i; b = j
            elif A[i, j] != 0:
                band = False
    if band:
        return a, b
    else:
        return None, None


def lindblad_operator(A, rho):
    r"""This function returns the action of a Lindblad operator A on a density\
 matrix rho. This is defined as :

    .. math::
        \mathcal{L}(A, \rho) = A \rho A^\dagger -
        (A^\dagger A \rho + \rho A^\dagger A)/2.

    >>> rho=define_density_matrix(3)
    >>> lindblad_operator( ketbra(1,2,3) ,rho )
    Matrix([
    [   rho22, -rho12/2,        0],
    [-rho21/2,   -rho22, -rho23/2],
    [       0, -rho32/2,        0]])

    """
    a, b = sigma_operator_indices(A)
    # print(111, a, b)
    if a is not None and b is not None:
        Ne = A.shape[0]
        L = zeros(Ne, Ne)
        L[a, a] += rho[b, b]

        for j in range(Ne):
            L[b, j] += -rho[b, j]/2

        for i in range(Ne):
            L[i, b] += -rho[i, b]/2

        return L
    else:
        return A*rho*A.adjoint() - (A.adjoint()*A*rho + rho*A.adjoint()*A)/2


def lindblad_terms(gamma, rho, Ne, verbose=1):
    u"""Return the Lindblad terms for decays gamma in matrix form.

    >>> from sympy import pprint
    >>> aux = define_frequencies(4, explicitly_antisymmetric=True)
    >>> omega_level, omega, gamma = aux
    >>> gamma = gamma.subs({gamma[2, 0]:0, gamma[3, 0]:0, gamma[3, 1]:0})
    >>> pprint(gamma, use_unicode=True)
    ⎡ 0   -γ₂₁   0     0  ⎤
    ⎢                     ⎥
    ⎢γ₂₁   0    -γ₃₂   0  ⎥
    ⎢                     ⎥
    ⎢ 0   γ₃₂    0    -γ₄₃⎥
    ⎢                     ⎥
    ⎣ 0    0    γ₄₃    0  ⎦
    >>> rho = define_density_matrix(4)
    >>> pprint(lindblad_terms(gamma, rho, 4), use_unicode=True)
    ⎡                -γ₂₁⋅ρ₁₂             -γ₃₂⋅ρ₁₃             -γ₄₃⋅ρ₁₄      ⎤
    ⎢ γ₂₁⋅ρ₂₂        ─────────            ─────────            ─────────     ⎥
    ⎢                    2                    2                    2         ⎥
    ⎢                                                                        ⎥
    ⎢-γ₂₁⋅ρ₂₁                          γ₂₁⋅ρ₂₃   γ₃₂⋅ρ₂₃    γ₂₁⋅ρ₂₄   γ₄₃⋅ρ₂₄⎥
    ⎢─────────  -γ₂₁⋅ρ₂₂ + γ₃₂⋅ρ₃₃   - ─────── - ───────  - ─────── - ───────⎥
    ⎢    2                                2         2          2         2   ⎥
    ⎢                                                                        ⎥
    ⎢-γ₃₂⋅ρ₃₁     γ₂₁⋅ρ₃₂   γ₃₂⋅ρ₃₂                         γ₃₂⋅ρ₃₄   γ₄₃⋅ρ₃₄⎥
    ⎢─────────  - ─────── - ───────  -γ₃₂⋅ρ₃₃ + γ₄₃⋅ρ₄₄   - ─────── - ───────⎥
    ⎢    2           2         2                               2         2   ⎥
    ⎢                                                                        ⎥
    ⎢-γ₄₃⋅ρ₄₁     γ₂₁⋅ρ₄₂   γ₄₃⋅ρ₄₂    γ₃₂⋅ρ₄₃   γ₄₃⋅ρ₄₃                     ⎥
    ⎢─────────  - ─────── - ───────  - ─────── - ───────       -γ₄₃⋅ρ₄₄      ⎥
    ⎣    2           2         2          2         2                        ⎦

    Notice that there are more terms than simply adding a decay
    gamma_ij*rho_ij/2 for each coherence.

    """
    # We count the necessary Lindblad operators.
    Nterms = 0
    for i in range(Ne):
        for j in range(i):
            if gamma[i, j] != 0:
                Nterms += 1

    L = zeros(Ne)
    counter = 0
    t0 = time()
    for i in range(Ne):
        for j in range(i):
            if gamma[i, j] != 0:
                counter += 1
                sig = ket(j+1, Ne)*bra(i+1, Ne)
                L += gamma[i, j]*lindblad_operator(sig, rho)
                tn = time()
                if tn-t0 > 1:
                    aux = "Calculated up to i={}, j={}, or {}/{} = {:2.2f} %."
                    if verbose > 0:
                        print(aux.format(i, j, counter, Nterms,
                                         float(counter+1)/Nterms*100))
                    t0 = tn
    return L


def define_psi_coefficients(Ne):
    ur"""Define the components of an arbitrary state vector.

    >>> from sympy import pprint
    >>> pprint(define_psi_coefficients(3), use_unicode=True)
    ⎛⎡c₁(t)⎤  ⎡\tilde{c}_{1}(t)⎤  ⎡θ₁⎤⎞
    ⎜⎢     ⎥  ⎢                ⎥  ⎢  ⎥⎟
    ⎜⎢c₂(t)⎥, ⎢\tilde{c}_{2}(t)⎥, ⎢θ₂⎥⎟
    ⎜⎢     ⎥  ⎢                ⎥  ⎢  ⎥⎟
    ⎝⎣c₃(t)⎦  ⎣\tilde{c}_{3}(t)⎦  ⎣θ₃⎦⎠

    """
    t = Symbol("t", real=True)
    c = Matrix([Function("c"+str(i+1))(t) for i in range(Ne)])
    ctilde = Matrix([Function(r"\tilde{c}_{"+str(i+1)+"}")(t)
                     for i in range(Ne)])
    phase = Matrix([Symbol("theta"+str(i+1), real=True) for i in range(Ne)])
    return c, ctilde, phase


def part_symbolic(z, s):
    r"""Extract the real or imaginary part of an expression.

    >>> rho = define_density_matrix(2)
    >>> part_symbolic(rho[1, 1], -1)
    0

    >>> part_symbolic(rho[1, 0], 1)
    re(rho21)

    """
    if s == 1: return re(z)
    else: return im(z)


def define_rho_vector(rho, Ne):
    u"""Define the vectorized density matrix.

    >>> from sympy import pprint
    >>> rho = define_density_matrix(3)
    >>> pprint(define_rho_vector(rho, 3), use_unicode=True)
    ⎡  ρ₂₂  ⎤
    ⎢       ⎥
    ⎢  ρ₃₃  ⎥
    ⎢       ⎥
    ⎢re(ρ₂₁)⎥
    ⎢       ⎥
    ⎢re(ρ₃₁)⎥
    ⎢       ⎥
    ⎢re(ρ₃₂)⎥
    ⎢       ⎥
    ⎢im(ρ₂₁)⎥
    ⎢       ⎥
    ⎢im(ρ₃₁)⎥
    ⎢       ⎥
    ⎣im(ρ₃₂)⎦

    """
    rho_vect = []
    for mu in range(1, Ne**2):
        i, j, s = IJ(mu, Ne)
        i = i-1; j = j-1
        rho_vect += [part_symbolic(rho[i, j], s)]
    return Matrix(rho_vect)


def calculate_A_b(eqs, unfolding, verbose=0):
    u"""Calculate the equations in matrix form.

    >>> from sympy import symbols, pprint, I
    >>> rho = define_density_matrix(2, explicitly_hermitian=True,
    ...                             normalized=True)

    >>> Omega = symbols("Omega")
    >>> delta = symbols("delta", real=True)
    >>> hbar = symbols("hbar", positive=True)
    >>> H = hbar*Matrix([[0, Omega.conjugate()/2], [Omega/2, -delta]])

    >>> Ne = 2
    >>> aux = define_frequencies(Ne, explicitly_antisymmetric=True)
    >>> omega_level, omega, gamma = aux

    >>> eqs = I/hbar*(rho*H-H*rho) + lindblad_terms(gamma, rho, 2)

    >>> from fast import Unfolding
    >>> unfolding = Unfolding(Ne, True, True, True)
    >>> A, b = calculate_A_b(eqs, unfolding)
    >>> pprint(A, use_unicode=True)
    ⎡ -γ₂₁   im(Ω)  -re(Ω)⎤
    ⎢                     ⎥
    ⎢        -γ₂₁         ⎥
    ⎢-im(Ω)  ─────    -δ  ⎥
    ⎢          2          ⎥
    ⎢                     ⎥
    ⎢               -γ₂₁  ⎥
    ⎢re(Ω)     δ    ───── ⎥
    ⎣                 2   ⎦

    >>> pprint(b, use_unicode=True)
    ⎡   0   ⎤
    ⎢       ⎥
    ⎢-im(Ω) ⎥
    ⎢───────⎥
    ⎢   2   ⎥
    ⎢       ⎥
    ⎢ re(Ω) ⎥
    ⎢ ───── ⎥
    ⎣   2   ⎦

    """
    Ne = unfolding.Ne
    Nrho = unfolding.Nrho
    lower_triangular = unfolding.lower_triangular
    rho = define_density_matrix(Ne, explicitly_hermitian=lower_triangular,
                                normalized=unfolding.normalized)

    rho_vect = unfolding(rho)
    if unfolding.real:
        ss_comp = {rho[i, j]: re(rho[i, j])+I*im(rho[i, j])
                   for j in range(Ne) for i in range(Ne)}

    A = []; b = []
    for mu in range(Nrho):
        s, i, j = unfolding.IJ(mu)
        if verbose > 0: print mu
        eq = part_symbolic(eqs[i, j].subs(ss_comp), s)
        eq_new = 0
        row = []
        for nu in range(Nrho):
            variable = rho_vect[nu]
            coefficient = Derivative(eq, variable).doit()
            row += [coefficient]
            eq_new += coefficient*variable

        b += [-(eq-eq_new).expand()]

        A += [row]
    A = Matrix(A); b = Matrix(b)
    return A, b


def phase_transformation(Ne, Nl, r, Lij, omega_laser, phase):
    r"""Obtain a phase transformation to eliminate explicit time dependence.

    >>> Ne = 2

    """
    ph = find_phase_transformation(Ne, Nl, r, Lij)

    return {phase[i]: sum([ph[i][j]*omega_laser[j] for j in range(Nl)])
            for i in range(Ne)}


def calculate_boundaries(Ne, Nl, r, Lij, omega_laser, phase):
    r"""Obtain a phase transformation to eliminate explicit time dependence.

    >>> Ne = 3
    >>> Nl = 2

    >>> r = define_r_components(Ne, helicity=True, explicitly_hermitian=True)
    >>> r = [ri.subs({r[0][2,0]:0,r[1][2,0]:0,r[2][2,0]:0}) for ri in r]

    >>> Lij = [[1,2,[1]],[2,3,[2]]]
    >>> from fast.misc import formatLij
    >>> Lij = formatLij(Lij,Ne)
    >>> E0, omega_laser = define_laser_variables(Nl)
    >>> c, ctilde, phase = define_psi_coefficients(Ne)
    >>> print phase_transformation(Ne, Nl, r, Lij, omega_laser, phase)
    {theta2: varpi_2, theta3: 0, theta1: varpi_1 + varpi_2}

    """
    ph = find_phase_transformation(Ne, Nl, r, Lij)

    return {phase[i]: sum([ph[i][j]*omega_laser[j] for j in range(Nl)])
            for i in range(Ne)}


def hamiltonian(Ep, epsilonp, detuning_knob, rm, omega_level, omega_laser, xi,
                RWA=True, RF=True):
    r"""Return symbolic Hamiltonian.

    >>> from sympy import zeros, pi, pprint, symbols
    >>> Ne = 3
    >>> Nl = 2
    >>> Ep, omega_laser = define_laser_variables(Nl)
    >>> epsilonp = [polarization_vector(0, -pi/2, 0, 0, 1) for l in range(Nl)]
    >>> detuning_knob = symbols("delta1 delta2", real=True)

    >>> xi = [zeros(Ne, Ne) for l in range(Nl)]
    >>> coup = [[(1, 0)], [(2, 0)]]
    >>> for l in range(Nl):
    ...     for pair in coup[l]:
    ...         xi[l][pair[0], pair[1]] = 1
    ...         xi[l][pair[1], pair[0]] = 1

    >>> rm = define_r_components(Ne, xi, explicitly_hermitian=True,
    ...                          helicity=True, p=-1)

    >>> rm = helicity_to_cartesian(rm)
    >>> omega_level, omega, gamma = define_frequencies(Ne, True)

    >>> H = hamiltonian(Ep, epsilonp, detuning_knob, rm, omega_level,
    ...                 omega_laser, xi, RWA=True, RF=False)

    >>> print H[1, 0]
    -E_{01}*e*r_{0;21}*exp(-I*t*varpi_1)/2
    >>> print H[2, 0]
    -E_{02}*e*r_{0;31}*exp(-I*t*varpi_2)/2
    >>> print H[2, 2]
    hbar*omega_3

    """
    # We check what RF is.
    if type(RF) == list:
        theta = RF[:]
        RF = True
    elif type(RF) == Matrix:
        theta = [RF[i, 0] for i in range(RF.shape[0])]
        RF = True
    elif RF:
        # theta should be calculate here!
        s = "We are still missing automatic calculation of phase "
        s += "transformations."
        raise ValueError(s)

    if not RWA and RF:
        s = "The rotating frame does not exist without the rotating wave \
approximation, as far as I know."
        raise ValueError(s)

    # We check that the epsilonp is a list of vectors.
    if not isinstance(epsilonp, list):
        raise ValueError("epsilonp must be a list of polarization vectors.")
    if not isinstance(epsilonp[0], Matrix):
        raise ValueError("epsilonp must be a list of polarization vectors.")

    Ne = len(omega_level)
    Nl = len(omega_laser)

    H = zeros(Ne, Ne)
    hbar, e = symbols("hbar e", positive=True)
    t = symbols("t", real=True)

    for i in range(Ne):
        for j in range(Ne):
            rmij = vector_element(rm, i, j)
            rpij = vector_element(rm, j, i).conjugate()
            for l in range(Nl):
                epsilonpl = epsilonp[l]
                epsilonml = epsilonpl.conjugate()

                if RF:
                    Epl = xi[l][i, j]*Ep[l]
                    Epl *= exp(-I*(theta[j]-theta[i]-t*omega_laser[l]))
                    Eml = xi[l][i, j]*Ep[l].conjugate()
                    Eml *= exp(-I*(theta[j]-theta[i]+t*omega_laser[l]))
                else:
                    Epl = Ep[l]*xi[l][i, j]*exp(-I*omega_laser[l]*t)
                    Eml = Epl.conjugate()

                # The E^(+)r^(-) term
                H[i, j] += -e*Epl/2*cartesian_dot_product(epsilonpl, rmij)
                # The E^(-)r^(+) term
                H[i, j] += -e*Eml/2*cartesian_dot_product(epsilonml, rpij)
                if not RWA:
                    # The E^(+)r^(+) term
                    H[i, j] += -e*Epl/2*cartesian_dot_product(epsilonpl, rpij)
                    # The E^(-)r^(-) term
                    H[i, j] += -e*Eml/2*cartesian_dot_product(epsilonml, rmij)

            if i == j:
                if RF:
                    H[i, j] += hbar*(omega_level[i]+diff(theta[i], t))
                else:
                    H[i, j] += hbar*omega_level[i]

    return H


###############################################################################
# 3D vector class definitions.
class Vector3D(QExpr):
    """Class for 3d vectors."""

    @classmethod
    def default_args(self):
        """Default arguments."""
        return ("v",)

    def _sympystr(self, printer, *args):
        contents = self._print_contents(printer, *args)
        return str(contents)

    def _latex(self, printer, *args):
        contents = self._print_contents_latex(printer, *args)
        # The extra {} brackets are needed to get matplotlib's latex
        # rendered to render this properly.
        return '{%s%s%s}' % (r"\vec{",
                             contents, "}")

    def dot(self, other):
        """Return the dot product."""
        return DotProduct(self, other)

    def cross(self, other):
        """Return the cross product."""
        return CrossProduct(self, other)

    def __getitem__(self, key):
        r"""Return an indexed instance if the vector is an iterable type."""
        arg = self.args[0]
        if isinstance(arg, conjugate):
            arg = arg.args[0]
            if hasattr(arg, '__getitem__'):
                return Vector3D(arg[key].conjugate())
            else:
                raise TypeError("This vector is not iterable.")

        if hasattr(arg, '__getitem__'):
            return Vector3D(arg[key])
        else:
            raise TypeError("This vector is not iterable.")

    def __mul__(self, other):
        r"""Scalar multiplication."""
        if isinstance(other, Vector3D):
            s = "Try methods dot or cross for scalar or vector products."
            raise TypeError(s)
        else:
            return other*self


class DotProduct(Expr):
    """An unevaluated dot product."""

    is_complex = True

    def __new__(cls, v1, v2):
        r"""Create a new dot product."""
        if not isinstance(v1, Vector3D):
            raise TypeError("first object must be an instance of Vector3D.")
        if not isinstance(v2, Vector3D):
            raise TypeError("second object must be an instance of Vector3D.")

        obj = Expr.__new__(cls, v1, v2)
        return obj

    @property
    def v1(self):
        r"""The left-hand side vector."""
        return self.args[0]

    @property
    def v2(self):
        r"""The right-hand side vector."""
        return self.args[1]

    def _sympystr(self, printer, *args):
        sv1 = self.v1._sympystr(printer, *args)
        sv2 = self.v2._sympystr(printer, *args)
        return 'dot(%s, %s)' % (sv1, sv2)

    def _latex(self, printer, *args):
        v1 = printer._print(self.v1, *args)
        v2 = printer._print(self.v2, *args)
        return r'{%s \cdot %s}' % (v1, v2)


class CrossProduct(Expr):
    """An unevaluated dot product."""

    is_complex = True

    def __new__(cls, v1, v2):
        r"""Create a new dot product."""
        if not isinstance(v1, Vector3D):
            raise TypeError("first object must be an instance of Vector3D.")
        if not isinstance(v2, Vector3D):
            raise TypeError("second object must be an instance of Vector3D.")

        obj = Expr.__new__(cls, v1, v2)
        return obj

    @property
    def v1(self):
        r"""The left-hand side vector."""
        return self.args[0]

    @property
    def v2(self):
        r"""The left-hand side vector."""
        return self.args[1]

    def _sympystr(self, printer, *args):
        sv1 = self.v1._sympystr(printer, *args)
        sv2 = self.v2._sympystr(printer, *args)
        return 'cross(%s,%s)' % (sv1, sv2)

    def _latex(self, printer, *args):
        v1 = printer._print(self.v1, *args)
        v2 = printer._print(self.v2, *args)
        return r'{%s \times %s}' % (v1, v2)


def dot(a, b):
    r"""Dot product of two 3d vectors."""
    if isinstance(a, Mul):
        a = a.expand()
        avect = 1
        aivect = -1
        for ai, fact in enumerate(a.args):
            if isinstance(fact, Vector3D):
                avect = fact
                aivect = ai
                break

        acoef = a.args[:aivect] + a.args[aivect+1:]
        acoef = Mul(*acoef)
        return acoef*dot(avect, b)

    if isinstance(b, Mul):
        b = b.expand()
        bvect = 1
        bivect = -1
        for bi, fact in enumerate(b.args):
            if isinstance(fact, Vector3D):
                bvect = fact
                bivect = bi
                break

        bcoef = b.args[:bivect] + b.args[bivect+1:]
        bcoef = Mul(*bcoef)
        return bcoef*dot(a, bvect)

    if isinstance(a, Vector3D) and isinstance(b, Vector3D):
        return DotProduct(a, b)

    if hasattr(a, "shape") and hasattr(b, "shape"):
        return cartesian_dot_product(a, b)

    print a, b, type(a), type(b),
    print isinstance(a, Vector3D), isinstance(b, Vector3D)
    raise NotImplementedError("could not catch these instances in dot!")


def cross(a, b):
    r"""Cross product of two 3d vectors."""
    if isinstance(a, Mul):
        a = a.expand()
        avect = 1
        aivect = -1
        for ai, fact in enumerate(a.args):
            if isinstance(fact, Vector3D):
                avect = fact
                aivect = ai
                break

        acoef = a.args[:aivect] + a.args[aivect+1:]
        acoef = Mul(*acoef)
        return acoef*cross(avect, b)

    if isinstance(b, Mul):
        b = b.expand()
        bvect = 1
        bivect = -1
        for bi, fact in enumerate(b.args):
            if isinstance(fact, Vector3D):
                bvect = fact
                bivect = bi
                break

        bcoef = b.args[:bivect] + b.args[bivect+1:]
        bcoef = Mul(*bcoef)
        return bcoef*cross(a, bvect)

    if isinstance(a, Vector3D) and isinstance(b, Vector3D):
        return CrossProduct(a, b)


if __name__ == "__main__":
    import doctest
    print doctest.testmod(verbose=False)
