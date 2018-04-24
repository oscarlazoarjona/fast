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
from atomic_structure import Atom, State
from sympy import Integer, zeros, eye, Matrix, KroneckerDelta, I, sqrt
from sympy.physics.wigner import clebsch_gordan
from sympy.physics.quantum import TensorProduct
from numpy import array, pi
from scipy.constants import physical_constants, hbar

muB = physical_constants["Bohr magneton"][0]


def lande_g_factors(element, isotope, L=None, J=None, F=None):
    r"""Return the Lande g-factors for a given atom or level.

    >>> element = "Rb"
    >>> isotope = 87
    >>> print lande_g_factors(element, isotope)
    [  9.99993686e-01   2.00231930e+00  -9.95141400e-04]

    The spin-orbit g-factor for a certain J
    >>> print lande_g_factors(element, isotope, L=0, J=1/Integer(2))
    [0.9999936864200584 2.0023193043622 -0.0009951414 2.00231930436220]

    The nuclear-coupled g-factor for a certain F
    >>> print lande_g_factors(element, isotope, L=0, J=1/Integer(2), F=1)
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
    ...     print f_group
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
    ...     print f_group
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


def permj(j1, j2):
    r"""Calculate the allowed total angular momenta.

    >>> from sympy import Integer
    >>> L = 1
    >>> S = 1/Integer(2)
    >>> permj(L, S)
    [1/2, 3/2]

    """
    jmin = abs(j1-j2)
    jmax = j1+j2
    return [jmin + i for i in range(jmax-jmin+1)]


def coupling_matrix_2j(j1, j2):
    ur"""For angular momenta $j_1, j_2$ the unitary transformation from the \
    uncoupled basis into the $j = j_1 \oplus j_2$ coupled basis.

    >>> from sympy import Integer, pprint
    >>> L = 0
    >>> S = 1/Integer(2)
    >>> pprint(coupling_matrix_2j(L, S))
    ⎡1  0⎤
    ⎢    ⎥
    ⎣0  1⎦

    >>> L = 1
    >>> S = 1/Integer(2)
    >>> pprint(coupling_matrix_2j(L, S))
    ⎡   -√6   √3             ⎤
    ⎢0  ────  ──   0    0   0⎥
    ⎢    3    3              ⎥
    ⎢                        ⎥
    ⎢             -√3   √6   ⎥
    ⎢0   0    0   ────  ──  0⎥
    ⎢              3    3    ⎥
    ⎢                        ⎥
    ⎢1   0    0    0    0   0⎥
    ⎢                        ⎥
    ⎢    √3   √6             ⎥
    ⎢0   ──   ──   0    0   0⎥
    ⎢    3    3              ⎥
    ⎢                        ⎥
    ⎢              √6   √3   ⎥
    ⎢0   0    0    ──   ──  0⎥
    ⎢              3    3    ⎥
    ⎢                        ⎥
    ⎣0   0    0    0    0   1⎦

    """
    # We calculate the quantum numbers for the uncoupled basis.
    M1 = [-j1 + i for i in range(2*j1+1)]
    M2 = [-j2 + i for i in range(2*j2+1)]
    j1j2nums = [(j1, m1, j2, m2) for m1 in M1 for m2 in M2]

    # We calculate the quantum numbers for the coupled basis.
    Jper = permj(j1, j2)
    jmjnums = [(J, MJ-J) for J in Jper for MJ in range(2*J+1)]

    # We build the transformation matrix.
    U = zeros((2*j1+1)*(2*j2+1))
    for ii, numj in enumerate(jmjnums):
        j, mj = numj
        for jj, numi in enumerate(j1j2nums):
            j1, m1, j2, m2 = numi
            U[ii, jj] = clebsch_gordan(j1, j2, j, m1, m2, mj)
    return U


def coupling_matrix_3j(j1, j2, j3):
    ur"""For angular momenta $j_1, j_2, j_3$ the unitary transformation from the \
    uncoupled basis into the $j = (j_1 \oplus j_2)\oplus j_3$ coupled basis.

    >>> from sympy import Integer, pprint
    >>> L = 0
    >>> S = 1/Integer(2)
    >>> II = 3/Integer(2)
    >>> pprint(coupling_matrix_3j(L, S, II))
    ⎡                     √3             ⎤
    ⎢0  -1/2   0     0    ──   0    0   0⎥
    ⎢                     2              ⎥
    ⎢                                    ⎥
    ⎢         -√2              √2        ⎥
    ⎢0   0    ────   0     0   ──   0   0⎥
    ⎢          2               2         ⎥
    ⎢                                    ⎥
    ⎢               -√3                  ⎥
    ⎢0   0     0    ────   0   0   1/2  0⎥
    ⎢                2                   ⎥
    ⎢                                    ⎥
    ⎢1   0     0     0     0   0    0   0⎥
    ⎢                                    ⎥
    ⎢    √3                              ⎥
    ⎢0   ──    0     0    1/2  0    0   0⎥
    ⎢    2                               ⎥
    ⎢                                    ⎥
    ⎢          √2              √2        ⎥
    ⎢0   0     ──    0     0   ──   0   0⎥
    ⎢          2               2         ⎥
    ⎢                                    ⎥
    ⎢                              √3    ⎥
    ⎢0   0     0    1/2    0   0   ──   0⎥
    ⎢                              2     ⎥
    ⎢                                    ⎥
    ⎣0   0     0     0     0   0    0   1⎦

    """
    idj3 = eye(2*j3+1)
    Jper = permj(j1, j2)
    U_Jj3_list = [coupling_matrix_2j(J, j3) for J in Jper]

    size = sum([U_Jj3_list[i].shape[0] for i in range(len(Jper))])
    U_Jj3 = zeros(size, size)
    ind0 = 0
    for i, U_Jj3i in enumerate(U_Jj3_list):
        sizeJ = U_Jj3i.shape[0]
        indf = ind0 + sizeJ
        U_Jj3[ind0: indf, ind0: indf] = U_Jj3_list[i]
        ind0 = indf
    return U_Jj3*TensorProduct(coupling_matrix_2j(j1, j2), idj3)


def angular_momentum_matrix(J, ind="z"):
    ur"""Return the angular momentum operator matrix (divided by hbar) for a\
    given J angular momentum.

    INPUT:

    -  ``ind`` - A string ("x", "y", "z", "all") indicating which direction \
    to calculate, or to return them all as :math:`(J_x, J_y, J_z)`.

    OUTPUT:

    - matrix forms of angular momentum operators in the basis \
    :math:`[|J, -J\rangle, \cdot, |J, J\rangle]`.

    >>> from sympy import Integer, pprint
    >>> pprint(angular_momentum_matrix(1/Integer(2)))
    ⎡-1/2   0 ⎤
    ⎢         ⎥
    ⎣ 0    1/2⎦

    >>> pprint(angular_momentum_matrix(1/Integer(2), "all"))
    ⎛            ⎡     ⅈ⎤             ⎞
    ⎜⎡ 0   1/2⎤, ⎢ 0   ─⎥, ⎡-1/2   0 ⎤⎟
    ⎜⎢        ⎥  ⎢     2⎥  ⎢         ⎥⎟
    ⎜⎣1/2   0 ⎦  ⎢      ⎥  ⎣ 0    1/2⎦⎟
    ⎜            ⎢-ⅈ    ⎥             ⎟
    ⎜            ⎢───  0⎥             ⎟
    ⎝            ⎣ 2    ⎦             ⎠

    >>> pprint(angular_momentum_matrix(1, "all"))
    ⎛⎡    √2    ⎤  ⎡         √2⋅ⅈ       ⎤            ⎞
    ⎜⎢0   ──  0 ⎥, ⎢  0      ────    0  ⎥, ⎡-1  0  0⎤⎟
    ⎜⎢    2     ⎥  ⎢          2         ⎥  ⎢        ⎥⎟
    ⎜⎢          ⎥  ⎢                    ⎥  ⎢0   0  0⎥⎟
    ⎜⎢√2      √2⎥  ⎢-√2⋅ⅈ           √2⋅ⅈ⎥  ⎢        ⎥⎟
    ⎜⎢──  0   ──⎥  ⎢──────    0     ────⎥  ⎣0   0  1⎦⎟
    ⎜⎢2       2 ⎥  ⎢  2              2  ⎥            ⎟
    ⎜⎢          ⎥  ⎢                    ⎥            ⎟
    ⎜⎢    √2    ⎥  ⎢        -√2⋅ⅈ       ⎥            ⎟
    ⎜⎢0   ──  0 ⎥  ⎢  0     ──────   0  ⎥            ⎟
    ⎝⎣    2     ⎦  ⎣          2         ⎦            ⎠

    """
    MJ = [-J+i for i in range(2*J+1)]

    if ind == "x":
        JX = Matrix([[sqrt((J-mj)*(J+mj+1))/2*KroneckerDelta(mi-1, mj)
                     for mj in MJ] for mi in MJ])
        JX += Matrix([[sqrt((J+mj)*(J-mj+1))/2*KroneckerDelta(mi+1, mj)
                      for mj in MJ] for mi in MJ])
        return JX
    elif ind == "y":
        JY = Matrix([[-I*sqrt((J-mj)*(J+mj+1))/2*KroneckerDelta(mi-1, mj)
                     for mj in MJ] for mi in MJ])
        JY += Matrix([[+I*sqrt((J+mj)*(J-mj+1))/2*KroneckerDelta(mi+1, mj)
                      for mj in MJ] for mi in MJ])
        return JY
    elif ind == "z":
        JZ = Matrix([[mi*KroneckerDelta(mi, mj) for mj in MJ] for mi in MJ])
        return JZ
    elif ind == "all":
        JX = angular_momentum_matrix(J, "x")
        JY = angular_momentum_matrix(J, "y")
        JZ = angular_momentum_matrix(J, "z")
        return JX, JY, JZ


def orbital_spin_nuclear_matrices(L, S, II, ind="z"):
    ur"""Return the matrix representation of the orbita, electron-spin, and \
    nuclear-spin angular momentum operators \
    :math:`\hat{\vec{L}}, \hat{\vec{L}}, \hat{\vec{L}}` in the coupled basis \
    :math:`[|J, -J\rangle, \cdot, |J, J\rangle]`.

    INPUT:

    -  ``ind`` - A string ("x", "y", "z", "all") indicating which direction \
    to calculate, or to return them all as :math:`(J_x, J_y, J_z)`.

    >>> from sympy import Integer, pprint
    >>> half = 1/Integer(2)
    >>> Lz, Sz, Iz = orbital_spin_nuclear_matrices(0, half, 3*half)
    >>> pprint(Lz)
    ⎡0  0  0  0  0  0  0  0⎤
    ⎢                      ⎥
    ⎢0  0  0  0  0  0  0  0⎥
    ⎢                      ⎥
    ⎢0  0  0  0  0  0  0  0⎥
    ⎢                      ⎥
    ⎢0  0  0  0  0  0  0  0⎥
    ⎢                      ⎥
    ⎢0  0  0  0  0  0  0  0⎥
    ⎢                      ⎥
    ⎢0  0  0  0  0  0  0  0⎥
    ⎢                      ⎥
    ⎢0  0  0  0  0  0  0  0⎥
    ⎢                      ⎥
    ⎣0  0  0  0  0  0  0  0⎦

    >>> pprint(Sz)
    ⎡                       √3                ⎤
    ⎢1/4   0    0     0     ──    0    0    0 ⎥
    ⎢                       4                 ⎥
    ⎢                                         ⎥
    ⎢ 0    0    0     0     0    1/2   0    0 ⎥
    ⎢                                         ⎥
    ⎢                                 √3      ⎥
    ⎢ 0    0   -1/4   0     0     0   ──    0 ⎥
    ⎢                                 4       ⎥
    ⎢                                         ⎥
    ⎢ 0    0    0    -1/2   0     0    0    0 ⎥
    ⎢                                         ⎥
    ⎢√3                                       ⎥
    ⎢──    0    0     0    -1/4   0    0    0 ⎥
    ⎢4                                        ⎥
    ⎢                                         ⎥
    ⎢ 0   1/2   0     0     0     0    0    0 ⎥
    ⎢                                         ⎥
    ⎢           √3                            ⎥
    ⎢ 0    0    ──    0     0     0   1/4   0 ⎥
    ⎢           4                             ⎥
    ⎢                                         ⎥
    ⎣ 0    0    0     0     0     0    0   1/2⎦

    >>> pprint(Iz)
    ⎡                        -√3                  ⎤
    ⎢-5/4   0     0     0    ────   0     0     0 ⎥
    ⎢                         4                   ⎥
    ⎢                                             ⎥
    ⎢ 0     0     0     0     0    -1/2   0     0 ⎥
    ⎢                                             ⎥
    ⎢                                    -√3      ⎥
    ⎢ 0     0    5/4    0     0     0    ────   0 ⎥
    ⎢                                     4       ⎥
    ⎢                                             ⎥
    ⎢ 0     0     0    -3/2   0     0     0     0 ⎥
    ⎢                                             ⎥
    ⎢-√3                                          ⎥
    ⎢────   0     0     0    -3/4   0     0     0 ⎥
    ⎢ 4                                           ⎥
    ⎢                                             ⎥
    ⎢ 0    -1/2   0     0     0     0     0     0 ⎥
    ⎢                                             ⎥
    ⎢            -√3                              ⎥
    ⎢ 0     0    ────   0     0     0    3/4    0 ⎥
    ⎢             4                               ⎥
    ⎢                                             ⎥
    ⎣ 0     0     0     0     0     0     0    3/2⎦

    >>> Lvec, Svec, Ivec = orbital_spin_nuclear_matrices(0, half, 0, "all")
    >>> pprint(Lvec)
    ⎡⎡0  0⎤, ⎡0  0⎤, ⎡0  0⎤⎤
    ⎢⎢    ⎥  ⎢    ⎥  ⎢    ⎥⎥
    ⎣⎣0  0⎦  ⎣0  0⎦  ⎣0  0⎦⎦

    >>> pprint(Svec)
    ⎡            ⎡     ⅈ⎤             ⎤
    ⎢⎡ 0   1/2⎤, ⎢ 0   ─⎥, ⎡-1/2   0 ⎤⎥
    ⎢⎢        ⎥  ⎢     2⎥  ⎢         ⎥⎥
    ⎢⎣1/2   0 ⎦  ⎢      ⎥  ⎣ 0    1/2⎦⎥
    ⎢            ⎢-ⅈ    ⎥             ⎥
    ⎢            ⎢───  0⎥             ⎥
    ⎣            ⎣ 2    ⎦             ⎦

    >>> pprint(Ivec)
    ⎡⎡0  0⎤, ⎡0  0⎤, ⎡0  0⎤⎤
    ⎢⎢    ⎥  ⎢    ⎥  ⎢    ⎥⎥
    ⎣⎣0  0⎦  ⎣0  0⎦  ⎣0  0⎦⎦

    """
    if ind == "all":
        LSIx = orbital_spin_nuclear_matrices(L, S, II, "x")
        LSIy = orbital_spin_nuclear_matrices(L, S, II, "y")
        LSIz = orbital_spin_nuclear_matrices(L, S, II, "z")
        return [[LSIx[i], LSIy[i], LSIz[i]] for i in range(3)]

    L0 = eye(2*L+1)
    S0 = eye(2*S+1)
    I0 = eye(2*II+1)

    Lind = angular_momentum_matrix(L, ind=ind)
    Sind = angular_momentum_matrix(S, ind=ind)
    Iind = angular_momentum_matrix(II, ind=ind)

    Lind = TensorProduct(TensorProduct(Lind, S0), I0)
    Sind = TensorProduct(TensorProduct(L0, Sind), I0)
    Iind = TensorProduct(TensorProduct(L0, S0), Iind)

    U_LSI = coupling_matrix_3j(L, S, II)

    Lind = U_LSI*Lind*U_LSI.adjoint()
    Sind = U_LSI*Sind*U_LSI.adjoint()
    Iind = U_LSI*Iind*U_LSI.adjoint()

    return Lind, Sind, Iind

if __name__ == "__main__":
    import doctest
    print doctest.testmod(verbose=False)
