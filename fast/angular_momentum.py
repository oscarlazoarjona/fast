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
"""This module contains all routines related to angular momentum."""

from sympy import zeros, eye, Matrix, sqrt, KroneckerDelta, I
from sympy.physics.wigner import clebsch_gordan
from sympy.physics.quantum import TensorProduct


def perm_j(j1, j2):
    r"""Calculate the allowed total angular momenta.

    >>> from sympy import Integer
    >>> L = 1
    >>> S = 1/Integer(2)
    >>> perm_j(L, S)
    [1/2, 3/2]

    """
    jmin = abs(j1-j2)
    jmax = j1+j2
    return [jmin + i for i in range(jmax-jmin+1)]


def perm_m(J):
    r"""Calculate the allowed magnetic quantum numbers for a given `J`.

    >>> from sympy import Integer
    >>> S = 1/Integer(2)
    >>> perm_m(S)
    [-1/2, 1/2]

    >>> perm_m(3)
    [-3, -2, -1, 0, 1, 2, 3]

    """
    return [-J+i for i in range(2*J+1)]


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
    Jper = perm_j(j1, j2)
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
    Jper = perm_j(j1, j2)
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


def spherical_tensor(Ji, Jj, K, Q):
    ur"""Return a matrix representation of the spherical tensor with quantum
    numbers $J_i, J_j, K, Q$.

    >>> from sympy import pprint
    >>> pprint(spherical_tensor(1, 1, 1, 0))
    ⎡-√2        ⎤
    ⎢────  0  0 ⎥
    ⎢ 2         ⎥
    ⎢           ⎥
    ⎢ 0    0  0 ⎥
    ⎢           ⎥
    ⎢         √2⎥
    ⎢ 0    0  ──⎥
    ⎣         2 ⎦

    >>> pprint(spherical_tensor(1, 2, 1, -1))
    ⎡      √10          ⎤
    ⎢0  0  ───   0    0 ⎥
    ⎢       10          ⎥
    ⎢                   ⎥
    ⎢           √30     ⎥
    ⎢0  0   0   ───   0 ⎥
    ⎢            10     ⎥
    ⎢                   ⎥
    ⎢                √15⎥
    ⎢0  0   0    0   ───⎥
    ⎣                 5 ⎦

    """
    keti = {(Ji, Mi): Matrix([KroneckerDelta(i, j)
            for j in range(2*Ji+1)])
            for i, Mi in enumerate(perm_m(Ji))}

    braj = {(Jj, Mj): Matrix([KroneckerDelta(i, j)
            for j in range(2*Jj+1)]).adjoint()
            for i, Mj in enumerate(perm_m(Jj))}

    if K not in perm_j(Ji, Jj):
        raise ValueError("K value is not allowed.")
    if Q not in perm_m(K):
        raise ValueError("Q value is not allowed.")

    Ni = 2*Ji+1
    Nj = 2*Jj+1
    T = zeros(Ni, Nj)
    for i, Mi in enumerate(perm_m(Ji)):
        for j, Mj in enumerate(perm_m(Jj)):
            T += (-1)**(Jj-Mj)*clebsch_gordan(Ji, Jj, K, Mi, -Mj, Q) * \
                keti[(Ji, Mi)]*braj[(Jj, Mj)]

    return T


if __name__ == "__main__":
    import doctest
    print doctest.testmod(verbose=False)
