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

"""

from sympy import Symbol, Matrix
from sympy import I, conjugate
from sympy import sin, cos, sqrt
from sympy import KroneckerDelta
from sympy import Function, Derivative
from sympy import re, im, zeros
from fast.misc import IJ, find_phase_transformation
from numpy import array as nparray
from numpy import sqrt as npsqrt


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
                          variables=None):
    r"""Return a symbolic density matrix.

    The arguments are

    Ne (integer):
        The number of atomic states.
    explicitly_hermitian (boolean):
        Whether to make $\rho_{ij}=\bar{\rho}_{ij}$ for $i<j$
    normalized (boolean):
        Whether to make $\rho_{11}=1-\sum_{i>1} \rho_{ii}$

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
        name = r"\rho"
        open_brace = "_{"
        close_brace = "}"
    else:
        comma = ""
        name = "rho"
        open_brace = ""
        close_brace = ""

    rho = []
    for i in range(Ne):
        row_rho = []
        for j in range(Ne):
            if i == j:
                row_rho += [define_symbol(name, open_brace, comma, i, j,
                                          close_brace, variables,
                                          positive=True)]
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
    ([E_0^1, E_0^2], [varpi_1, varpi_2])

    The amplitudes are complex by default:
    >>> conjugate(E0[0])
    conjugate(E_0^1)

    But they can optionally be made real:
    >>> E0, omega_laser = define_laser_variables(2, real_amplitudes=True)
    >>> conjugate(E0[0])
    E_0^1

    They can also be made explicit functions of given variables:
    >>> from sympy import symbols
    >>> t, z = symbols("t, z", real=True)
    >>> E0, omega_laser = define_laser_variables(2, variables=[t, z])
    >>> E0
    [E_0^1(t, z), E_0^2(t, z)]

    """
    if variables is None:
        E0 = [Symbol(r"E_0^"+str(l+1), real=real_amplitudes)
              for l in range(Nl)]
    else:
        E0 = [Function(r"E_0^"+str(l+1), real=real_amplitudes)(*variables)
              for l in range(Nl)]

    omega_laser = [Symbol(r"varpi_"+str(l+1), positive=True)
                   for l in range(Nl)]
    return E0, omega_laser


def polarization_vector(phi, theta, alpha, beta, p):
    r"""This function returns a unitary vector describing the polarization
    of plane waves. It recieves as arguments:

    phi   .- The spherical coordinates azimuthal angle of the wave vector k.
    theta .- The spherical coordinates polar angle of the wave vector k.
    alpha .- The rotation of a half-wave plate.
    beta  .- The rotation of a quarter-wave plate.
    p     .- either 1 or -1 to indicate whether to return epsilon^(+) or
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

    """
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

    return R3*R2*R1*epsilon


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
                a = -ap*em +a0*e0 -am*ep

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

    """
    if numeric:
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
    if numeric:
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
    return -v1[2]*v2[0] + v1[1]*v2[1] - v1[0]*v2[2]


def cartesian_dot_product(v1, v2):
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]


def define_r_components(Ne, explicitly_hermitian=False, helicity=False,
                        real=True, p=None):
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

    return r


def vector_element(r, i, j):
    return Matrix([r[p][i, j] for p in range(3)])


def define_frequencies(Ne, explicitly_antisymmetric=False):

    omega_level = [Symbol('omega_'+str(i+1), real=True) for i in range(Ne)]

    if Ne > 9:
        comma = ","
        open_brace = "{"
        close_brace = "}"
    else:
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
                om = Symbol(r"omega_" +
                            open_brace+str(i+1)+comma+str(j+1) +
                            close_brace, real=True)
                ga = Symbol(r"gamma_" +
                            open_brace+str(i+1)+comma+str(j+1) +
                            close_brace, real=True)
            elif explicitly_antisymmetric:
                om = -Symbol(r"omega_" +
                             open_brace+str(j+1)+comma+str(i+1) +
                             close_brace, real=True)
                ga = -Symbol(r"gamma_" +
                             open_brace+str(j+1)+comma+str(i+1) +
                             close_brace, real=True)
            else:
                om = Symbol(r"omega_" +
                            open_brace+str(i+1)+comma+str(j+1) +
                            close_brace, real=True)
                ga = Symbol(r"gamma_" +
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
    if i > j: return 1
    else: return 0


def delta_lesser(i, j):
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
    r"""This function returns the outer product |i><j| where |i> and |j>
    are elements of the canonical basis of an Ne-dimensional Hilbert space
    (in matrix form).

    >>> ketbra(2,3,3)
    Matrix([
    [0, 0, 0],
    [0, 0, 1],
    [0, 0, 0]])

    """
    return ket(i, Ne)*bra(j, Ne)


def lindblad_operator(A, rho):
    r"""This function returns the action of a Lindblad operator A on a density
    matrix rho. This is defined as :
        L(A,rho) = A*rho*A.adjoint()
                 - (A.adjoint()*A*rho + rho*A.adjoint()*A)/2.

    >>> rho=define_density_matrix(3)
    >>> lindblad_operator( ketbra(1,2,3) ,rho )
    Matrix([
    [   rho22, -rho12/2,        0],
    [-rho21/2,   -rho22, -rho23/2],
    [       0, -rho32/2,        0]])

    """
    return A*rho*A.adjoint() - (A.adjoint()*A*rho + rho*A.adjoint()*A)/2


def lindblad_terms(gamma, rho, Ne):
    L = zeros(Ne)
    for i in range(Ne):
        for j in range(i):
            L += gamma[i, j]*lindblad_operator(ket(j+1, Ne)*bra(i+1, Ne), rho)
    return L


def define_psi_coefficients(Ne):
    t = Symbol("t", real=True)
    c = Matrix([Function("c"+str(i+1))(t) for i in range(Ne)])
    ctilde = Matrix([Function(r"\tilde{c}_{"+str(i+1)+"}")(t)
                     for i in range(Ne)])
    phase = Matrix([Symbol("theta"+str(i+1)) for i in range(Ne)])
    return c, ctilde, phase


def part_symbolic(z, s):
    if s == 1: return re(z)
    else: return im(z)


def define_rho_vector(rho, Ne):
    rho_vect = []
    for mu in range(1, Ne**2):
        i, j, s = IJ(mu, Ne)
        i = i-1; j = j-1
        rho_vect += [part_symbolic(rho[i, j], s)]
    return Matrix(rho_vect)


def calculate_A_b(eqs, rho, Ne):
    rho_vect = define_rho_vector(rho, Ne)
    A = []; b = []
    ss_comp = {rho[i, j]: re(rho[i, j])+I*im(rho[i, j])
               for j in range(Ne) for i in range(Ne)}

    for mu in range(1, Ne**2):
        i, j, s = IJ(mu, Ne)
        ii = i-1; jj = j-1
        # print ii,jj,s
        eq = part_symbolic(eqs[ii, jj].subs(ss_comp), s)
        eq_new = 0
        row = []
        for nu in range(1, Ne**2):
            variable = rho_vect[nu-1]
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
    {theta3: 0, theta1: varpi_1 + varpi_2, theta2: varpi_2}

    """
    ph = find_phase_transformation(Ne, Nl, r, Lij)

    return {phase[i]: sum([ph[i][j]*omega_laser[j] for j in range(Nl)])
            for i in range(Ne)}
