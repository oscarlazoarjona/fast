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

r"""A general module for functions that are needed in more than one module."""

# from colorsys import hls_to_rgb, hsv_to_rgb
# from scipy.optimize import curve_fit
# from matplotlib import pyplot
import numpy as np
from time import time
import os
from math import sqrt, exp
from sympy import solve, Symbol, diff, pprint
from sympy import zeros as symzeros


use_netcdf = False
if use_netcdf:
    from netCDF4 import Dataset

sage_included = 'sage' in globals().keys()


def fast_linear_system(M, b=None, real=False):
    u"""Return a fast function that solves the equation dx/dt = Mx - b.

    Parameters
    ----------
    ``M`` : numpy array with shape (N, N)
        The coefficient matrix of the system.

    ``b`` : numpy array with shape (N)
        The independent terms vector.

    ``real`` bool:
        Whether to return a real array.

    Returns
    -------
    function
        A function that takes as arguments a numpy array of time values, and an
        initial vector.

    Examples
    --------
    We will for a simple Rabi oscillation for two different vectorizations.
    We first define a few things to prove that our proposed solution is
    correct.

    >>> from fast import (define_density_matrix, Unfolding, calculate_A_b,
    ...     sharp, flat, ket)
    >>> from sympy import (init_printing, symbols, Matrix, sqrt, zeros, diff,
    ...     I, Eq, exp, simplify, pprint)
    >>> import numpy as np
    >>> Ne = 2
    >>> t = symbols("t", real=True)
    >>> rho1 = define_density_matrix(Ne, explicitly_hermitian=False,
    ...     normalized=False, variables=[t])
    >>> unfolding1 = Unfolding(Ne, real=False, lower_triangular=False,
    ...     normalized=False)

    >>> rho2 = define_density_matrix(Ne, explicitly_hermitian=True,
    ...     normalized=True)
    >>> unfolding2 = Unfolding(Ne, real=True, lower_triangular=True,
    ...     normalized=True)

    We define a Hamiltonian operator.
    >>> psi1 = Matrix([1, 1])/sqrt(2)
    >>> psi2 = Matrix([1, -1])/sqrt(2)
    >>> v = [psi1, psi2]
    >>> h = [2, 4]
    >>> H = sum([h[i]*v[i]*v[i].adjoint() for i in range(2)], zeros(2))
    >>> pprint(H, use_unicode=True)
    ⎡3   -1⎤
    ⎢      ⎥
    ⎣-1  3 ⎦

    The equations are
    >>> eqs1_lhs = diff(rho1, t)
    >>> eqs1_rhs = (I*(rho1*H - H*rho1)).expand()
    >>> eqs2_rhs = (I*(rho2*H - H*rho2)).expand()

    On the left-hand side
    >>> pprint(eqs1_lhs)
    ⎡d           d         ⎤
    ⎢──(ρ₁₁(t))  ──(ρ₁₂(t))⎥
    ⎢dt          dt        ⎥
    ⎢                      ⎥
    ⎢d           d         ⎥
    ⎢──(ρ₂₁(t))  ──(ρ₂₂(t))⎥
    ⎣dt          dt        ⎦

    And on the right-hand side
    >>> pprint(eqs1_rhs, use_unicode=True)
    ⎡-ⅈ⋅ρ₁₂(t) + ⅈ⋅ρ₂₁(t)  -ⅈ⋅ρ₁₁(t) + ⅈ⋅ρ₂₂(t)⎤
    ⎢                                          ⎥
    ⎣ⅈ⋅ρ₁₁(t) - ⅈ⋅ρ₂₂(t)   ⅈ⋅ρ₁₂(t) - ⅈ⋅ρ₂₁(t) ⎦


    With the vectorization
    >>> pprint(unfolding1(rho1), use_unicode=True)
    ⎡ρ₁₁(t)⎤
    ⎢      ⎥
    ⎢ρ₂₂(t)⎥
    ⎢      ⎥
    ⎢ρ₂₁(t)⎥
    ⎢      ⎥
    ⎣ρ₁₂(t)⎦

    The equations are `d rho/dt = M1 rho` with
    >>> M1 = sharp(I*H, unfolding1) - flat(I*H, unfolding1)
    >>> pprint(M1, use_unicode=True)
    ⎡0   0   ⅈ   -ⅈ⎤
    ⎢              ⎥
    ⎢0   0   -ⅈ  ⅈ ⎥
    ⎢              ⎥
    ⎢ⅈ   -ⅈ  0   0 ⎥
    ⎢              ⎥
    ⎣-ⅈ  ⅈ   0   0 ⎦

    And with vectorization
    >>> pprint(unfolding2(rho2), use_unicode=True)
    ⎡  ρ₂₂  ⎤
    ⎢       ⎥
    ⎢re(ρ₂₁)⎥
    ⎢       ⎥
    ⎣im(ρ₂₁)⎦

    The equations are `d rho/dt = M2 rho - b2` with
    >>> M2, b2 = calculate_A_b(eqs2_rhs, unfolding2)
    >>> pprint([M2, b2], use_unicode=True)
    ⎡⎡0   0  2⎤  ⎡0 ⎤⎤
    ⎢⎢        ⎥  ⎢  ⎥⎥
    ⎢⎢0   0  0⎥, ⎢0 ⎥⎥
    ⎢⎢        ⎥  ⎢  ⎥⎥
    ⎣⎣-2  0  0⎦  ⎣-1⎦⎦

    And a solution is
    >>> U = sum([exp(-I*t*h[i])*v[i]*v[i].adjoint()
    ...     for i in range(2)], zeros(2))
    >>> psi0 = 3*ket(1, Ne) + 4*ket(2, Ne)
    >>> psi0 = psi0.normalized()
    >>> rho0 = psi0*psi0.adjoint()
    >>> rhot = simplify(U*rho0*U.adjoint())
    >>> pprint(rhot, use_unicode=True)
    ⎡   7⋅cos(2⋅t)   1    7⋅ⅈ⋅sin(2⋅t)   12⎤
    ⎢ - ────────── + ─    ──────────── + ──⎥
    ⎢       50       2         50        25⎥
    ⎢                                      ⎥
    ⎢  7⋅ⅈ⋅sin(2⋅t)   12   7⋅cos(2⋅t)   1  ⎥
    ⎢- ──────────── + ──   ────────── + ─  ⎥
    ⎣       50        25       50       2  ⎦

    We substitute in the equations to prove that the solution is correct.
    >>> Eq(simplify(diff(unfolding1(rhot), t).expand()),
    ...     simplify(M1*unfolding1(rhot)).expand())
    True

    We will now test that the generated solvers return this solution. Numeric
    versions of the dynamical matrices are

    >>> M1_num = np.array([[complex(M1[i, j]) for j in range(Ne**2)] for i in
    ...     range(Ne**2)])
    >>> print(M1_num)
    [[0.+0.j 0.+0.j 0.+1.j 0.-1.j]
     [0.+0.j 0.+0.j 0.-1.j 0.+1.j]
     [0.+1.j 0.-1.j 0.+0.j 0.+0.j]
     [0.-1.j 0.+1.j 0.+0.j 0.+0.j]]
    >>> M2_num = np.array([[float(M2[i, j]) for j in range(Ne**2-1)] for i in
    ...     range(Ne**2-1)])
    >>> b2_num = np.array([float(b2[i])  for i in range(Ne**2-1)])
    >>> print(M2_num)
    [[ 0.  0.  2.]
     [ 0.  0.  0.]
     [-2.  0.  0.]]
    >>> print(b2_num)
    [ 0.  0. -1.]

    The initial conditions:
    >>> Nt = 101
    >>> tnum = np.linspace(0, 2*np.pi, Nt)

    >>> rho10_num = unfolding1(rho0)
    >>> rho10_num = np.array([complex(rho10_num[i]) for i in range(Ne**2)])

    >>> rho20_num = unfolding2(rho0)
    >>> rho20_num = np.array([float(rho20_num[i]) for i in range(Ne**2-1)])

    We create the solvers!
    >>> linear_system1 = fast_linear_system(M1_num)
    >>> linear_system2 = fast_linear_system(M2_num, b2_num, real=True)

    We call the solvers
    >>> rho1_num = linear_system1(tnum, rho10_num)
    >>> rho2_num = linear_system2(tnum, rho20_num)

    We extract the physically meaningful parts:
    >>> rho11re_1 = np.real(rho1_num[0, :])
    >>> rho22re_1 = np.real(rho1_num[1, :])
    >>> rho21re_1 = np.real(rho1_num[2, :])
    >>> rho21im_1 = np.imag(rho1_num[2, :])

    >>> rho22re_2 = np.real(rho2_num[0, :])
    >>> rho11re_2 = 1-rho22re_2
    >>> rho21re_2 = np.real(rho2_num[1, :])
    >>> rho21im_2 = np.real(rho2_num[2, :])

    We calculate the same parts from our analytic solution before.
    >>> rho22re_a = 7*np.cos(2*tnum)/50.0 + 0.5
    >>> rho11re_a = 1-rho22re_a
    >>> rho21re_a = 12.0/25
    >>> rho21im_a = -7*np.sin(2*tnum)/50.0

    We check that the solutions are the same.

    >>> err11 = np.abs(rho11re_a-rho11re_1)
    >>> err21 = np.abs(rho22re_a-rho22re_1)
    >>> err31 = np.abs(rho21re_a-rho21re_1)
    >>> err41 = np.abs(rho21im_a-rho21im_1)

    >>> err12 = np.abs(rho11re_a-rho11re_2)
    >>> err22 = np.abs(rho22re_a-rho22re_2)
    >>> err32 = np.abs(rho21re_a-rho21re_2)
    >>> err42 = np.abs(rho21im_a-rho21im_2)

    >>> assert np.amax(err11) < 1e-10
    >>> assert np.amax(err21) < 1e-10
    >>> assert np.amax(err31) < 1e-10
    >>> assert np.amax(err41) < 1e-10
    >>> assert np.amax(err12) < 1e-10
    >>> assert np.amax(err22) < 1e-10
    >>> assert np.amax(err32) < 1e-10
    >>> assert np.amax(err42) < 1e-10


    """
    NN = M.shape[0]
    lam, P = np.linalg.eig(M)
    Pinv = np.linalg.inv(P)

    if b is not None:
        d = np.dot(Pinv, b)
        for i in range(NN):
            if lam[i] != 0:
                d[i] = d[i]/lam[i]
            elif d[i] == 0:
                pass
            else:
                mes = "Having (P^{-1} b)_i != 0 and lambda_i = 0 leads to"
                mes += " infinite growth."
                raise ValueError(mes)

    def linear_system(t, x0):
        r"""Solve the linear system at time `t` from initial condition `x0`.
        Parameters
        ----------
        ``t`` : numpy array with shape (Nt)
            The time values to evaluate the system.

        ``x0`` : numpy array with shape (N)
            The initial condition.

        Returns
        -------
        numpy array with shape (N, Nt)
            An array with the solution at times `t`.

        """
        Nt = t.shape[0]
        yt = np.zeros((NN, Nt), complex)
        y0 = np.dot(Pinv, x0)
        if b is not None:
            y0 += -d

        for i in range(NN):
            yt[i, :] = y0[i]*np.exp(lam[i]*t)
            if b is not None:
                yt[i, :] += d[i]

        if real:
            return np.real(np.dot(P, yt))
        else:
            return np.dot(P, yt)

    return linear_system


def fprint(expr, print_ascii=False):
    r"""This function chooses whether to use ascii characters to represent
    a symbolic expression in the notebook or to use sympy's pprint.

    >>> from sympy import cos
    >>> omega=Symbol("omega")
    >>> fprint(cos(omega),print_ascii=True)
    cos(omega)


    """
    if print_ascii:
        pprint(expr, use_unicode=False, num_columns=120)
    else:
        return expr


def format_double(num):
    r"""Format a Python float into a Fortran double."""
    num = str(num)
    if 'e' in num:
        return num.replace('e', 'd')
    else:
        return num+'d0'


def Mu(i, j, s, N, excluded_mu=[]):
    """This function calculates the global index mu for the element i, j.

    It returns the index for the real part if s=0 and the one for the
    imaginary part if s=1.
    """
    if i == j:
        if s == -1:
            if i == 1:
                return 0
            else:
                mes = 'There is no population rhoii with i='+str(i)+'.'
                raise ValueError(mes)
        mu = i-1
    elif i > j and s == 1:
        mu = i - j + sum([N-k for k in range(1, j)])+N-1
    elif i > j and s == -1:
        mu = i - j + sum([N-k for k in range(1, j)])+N-1 + N*(N-1)/2
    else:
        mes = 'i='+str(i)+', j='+str(j)
        mes += ' Equations for i<j are not calculated.'+str(s)
        raise ValueError(mes)

    if excluded_mu != []:
        mu = mu-len([ii for ii in excluded_mu if ii < mu])
    return mu


def IJ(mu, N):
    """Return i, j, s for any given mu."""
    if mu == 0: return 1, 1, 1

    if mu not in range(0, N**2):
        raise ValueError('mu has an invalid value mu='+str(mu)+'.')

    if 1 <= mu <= N-1:
        return mu+1, mu+1, 1
    else:
        m = N-1
        M = N*(N-1)/2
        for jj in range(1, N):
            for ii in range(jj+1, N+1):
                m += 1
                if m == mu or m+M == mu:
                    if mu > N*(N+1)/2 - 1:
                        return ii, jj, -1
                    else:
                        return ii, jj, 1


def read_result(path, name, i=None, j=None, s=0, N=None,
                excluded_mu=[], use_netcdf=True, clone=None):
    r"""This function reads the results stored in path under file name.dat
    returning them as a list of N^2 lists of the form
    [frequency, rho22, rho33, ... rho_N,N-1]. Alternatively it can return only
    two lists [frequency, rho_i,j,s] where s must be 1 for the real part and
    -1 for the imaginary part.
    """
    if clone is not None:
        clone = '_'+str(clone)
    else:
        clone = ''

    if use_netcdf:
        ncfile = Dataset(path+name+clone+'.nc', 'r')
        matrix = ncfile.variables['matrix'][:]
        vector = ncfile.variables['vector'][:]
        Nrho, n_points = len(matrix), len(matrix[0])
        ncfile.close()

        if i is not None:
            mu = Mu(i, j, s, N, excluded_mu)
            if sage_included:
                return [(vector[k], matrix[mu-1][k]) for k in range(n_points)]
            else:
                return vector, matrix[mu-1]
        else:
            return [vector] + [matrix[k] for k in range(Nrho)]
        # ln=[for i in range(m)]
    else:
        r = file(path+name+clone+'.dat', 'r')
        l = r.readlines()
        r.close()

        try:
            ln = [[float(num) for num in li.split()] for li in l]
        except:
            print(num)
            print([l[0]])
            print([l[-2]])
            print([l[-1]])
            raise ValueError

        if i is not None:
            mu = Mu(i, j, s, N, excluded_mu)
            if sage_included:
                return [(li[0], li[mu]) for li in ln]
            else:
                x = [li[0] for li in ln]
                y = [li[mu] for li in ln]
                return x, y
        else:
            dat = [[ln[ii][0] for ii in range(len(ln))]]
            dat += [[ln[ii][mu+1] for ii in range(len(ln))]
                    for mu in range(N**2-1)]
            return dat


def dft(s, max_freq=False):
    r"""Calculate a discrete Fourier transform."""
    f = [i[1] for i in s]
    T = s[-1][0]
    Om = 1/T

    fn = np.fft.fft(f)
    fn = [sqrt(np.real(i)**2 + np.imag(i)**2) for i in fn]
    m = int(len(fn)/2)
    fn = fn[m:]+fn[:m]

    omm = Om*m

    sf = [(i*Om-omm, fn[i]) for i in range(len(s))]
    if max_freq:
        m = max(fn)
        im = fn.index(m)
        return sf[im][0]
    return sf


def formatLij(Lij0, Ne):
    """This function transforms a list of laser conections of the form
    [i,j,[l1,l2,...]] between states i and j by lasers l1,l2,... into a
    Ne x Ne matrix whose elements are the lasers connecting the corresponding
    indices.
    """
    # We create Lij as a matrix of lists of laser indices
    global Lij
    Lij = []
    for i in range(Ne):
        fila = []
        for j in range(Ne):
            band = False
            for tri in Lij0:
                if [i+1, j+1] == tri[:2]:
                    band = True
                    break
                elif [j+1, i+1] == tri[:2]:
                    band = True
                    break
            if band:
                fila += [tri[2]]
            else:
                fila += [[]]
        Lij += [fila]
    return Lij


def find_phase_transformation(Ne, Nl, r, Lij, verbose=0,
                              return_equations=False, **kwds):
    """This function returns a phase transformation specified as a list of
    lenght Ne whose elements correspond to each theta_i. Each element is a
    list of length Nl which specifies the coefficients multiplying each
    optical frequency omega^l. So for instance [[1,1],[1,0],[0,0]] means

    theta_1=omega^1+omega^2
    theta_2=omega^1
    theta_3=0.
    """
    # We first define the needed variables
    _omega_laser = [Symbol('omega_laser'+str(l+1)) for l in range(Nl)]
    _theta = [Symbol('theta'+str(i+1)) for i in range(Ne)]

    # We find all the equations that the specified problem has to fulfil.
    eqs = []
    for i in range(Ne):
        for j in range(i+1, Ne):
            if type(r[0]) == list:
                if (r[0][i][j] != 0) or (r[1][i][j] != 0) or (r[2][i][j] != 0):
                    for l in range(Nl):
                        if l+1 in Lij[i][j]:
                            eqs += [_omega_laser[l] + _theta[j] - _theta[i]]
            else:
                if (r[0][i, j] != 0) or (r[1][i, j] != 0) or (r[2][i, j] != 0):
                    for l in range(Nl):
                        if l+1 in Lij[i][j]:
                            eqs += [_omega_laser[l] + _theta[j] - _theta[i]]

    if return_equations:
        return eqs

    sol = solve(eqs, _theta, dict=True)[0]
    for i in range(Ne):
        if _theta[i] not in sol.keys():
            sol.update({_theta[i]: _theta[i]})

    sol_simple = {_theta[i]: sol[_theta[i]]-sol[_theta[-1]] for i in range(Ne)}

    sol = []
    for i in range(Ne):
        soli = []
        for l in range(Nl):
            soli += [diff(sol_simple[_theta[i]], _omega_laser[l])]
        sol += [soli]

    return sol


def calculate_iI_correspondence(omega):
    r"""Get the correspondance between degenerate and nondegenerate schemes."""
    Ne = len(omega[0])
    om = omega[0][0]
    correspondence = []
    I = 0
    for i in range(Ne):
        if omega[i][0] != om:
            om = omega[i][0]
            I += 1
        correspondence += [(i+1, I+1)]
    Nnd = I+1

    def I_nd(i):
        return correspondence[i-1][1]

    def i_d(I):
        for i in range(Ne):
            if correspondence[i][1] == I:
                return correspondence[i][0]

    return i_d, I_nd, Nnd


def part(z, s):
    r"""Get the real or imaginary part of a complex number."""
    if sage_included:
        if s == 1: return np.real(z)
        elif s == -1: return np.imag(z)
        elif s == 0:
            return z
    else:
        if s == 1: return z.real
        elif s == -1: return z.imag
        elif s == 0: return z


def detuning_combinations(lists):
    r"""This function recieves a list of length Nl with the number of
    transitions each laser induces. It returns the cartesian product of all
    these posibilities as a list of all possible combinations.
    """
    Nl = len(lists)
    comb = [[i] for i in range(lists[0])]
    for l in range(1, Nl):
        combn = []
        for c0 in comb:
            for cl in range(lists[l]):
                combn += [c0[:]+[cl]]
        comb = combn[:]
    return comb


def laser_detunings(Lij, Nl, i_d, I_nd, Nnd):
    r"""This function returns the list of transitions i,j that each laser
    produces as lists of length Ne, whose elements are all zero except for the
    ith element =1 and the jth element = -1.

    Also, it returns detuningsij, which contains the same information
    but as pairs if ij. The indices in it start from 0.
    """
    Ne = len(Lij)
    detunings = [[] for i in range(Nl)]
    detuningsij = [[] for i in range(Nl)]
    detuning = [0 for i in range(Nnd)]
    for i in range(1, Ne):
        for j in range(i):
            for l in Lij[i][j]:
                det = detuning[:]
                det[I_nd(i+1)-1] += 1; det[I_nd(j+1)-1] -= 1
                if det not in detunings[l-1]:
                    detunings[l-1] += [det]
                    detuningsij[l-1] += [(I_nd(i+1)-1, I_nd(j+1)-1)]

    return detunings, detuningsij


def Theta(i, j, theta, omega_rescaled, omega_min,
          detunings, detuningsij, combinations, detuning_indices,
          Lij, i_d, I_nd, Nnd,
          states=None, verbose=1, other_the=None):
    r"""This function returns code for Theta_i j as defined in the equation
    labeled Theta. in terms of detunings. It recieves indexes i,j starting
    from 1.
    """
    if i == j:
        return ''
    elif j > i:
        raise ValueError('i should never be less than j.')

    Ne = len(omega_rescaled[0]); Nl = len(theta[0])
    if other_the is not None:
        the = other_the[:]
    else:
        the = [theta[j-1][l]-theta[i-1][l] for l in range(Nl)]

    ##########################################################################
    # This part is about finding detunings that reduce to -omega_ij
    if the == [0 for l in range(Nl)]:
        if omega_rescaled[i-1][j-1] == 0:
            return ''
        # We need a combination of detunings that yields -omega_i,j.
        # According to equation labeled "detuning-exception1"
        # -omega_ij= delta^l_ik - delta^l_jk
        # so we seek a k that is lower than i and j
        # and an l such that l is in L_ik and in L_jk

        ii = min(i, j)
        band1 = False
        band2 = False
        for k in range(1, ii):
            for l in range(1, Nl+1):
                if l in Lij[i-1][k-1] and l in Lij[j-1][k-1]:
                    s = 'I found that -omega_'+str(i)+','+str(j)
                    s += '= delta^'+str(l)+'_'+str(i)+','+str(k)
                    s += '-delta^'+str(l)+'_'+str(j)+','+str(k)
                    s += '\n and also l='+str(l)
                    s += ' is in  L'+str(i)+','+str(k)+'='+str(Lij[i-1][k-1])
                    s += ' and in  L'+str(j)+','+str(k)+'='+str(Lij[j-1][k-1])
                    s += '\n'

                    band1 = True
                    break
            if band1:
                break

        if not band1:
            # If no such a k exists then we follow te equation labeled
            # "detuning-exception2"
            # -omega_ij= delta^l_kj - delta^l_ki
            # to look for a k greater than i and j
            # and an l such that l is in L_ik and in L_jk
            for k in range(ii+1, Ne+1):
                for l in range(1, Nl+1):
                    if l in Lij[i-1][k-1] and l in Lij[j-1][k-1]:
                        s = 'I found that -omega_'+str(i)+','+str(j)
                        s += '= delta^'+str(l)+'_'+str(k)+','+str(j)
                        s += '-delta^'+str(l)+'_'+str(k)+','+str(i)
                        s += '\n and also l='+str(l)
                        s += ' is in  L'+str(i)+','+str(k)+'='
                        s += str(Lij[i-1][k-1])
                        s += ' and in  L'+str(j)+','+str(k)+'='
                        s += str(Lij[j-1][k-1])+'\n'
                        band2 = True
                        break

                if band2:
                    break

        if band1:
            The = ''
            # We need to find which detunings are delta^l_i,k and delta^l_j,k.
            # Since they are detunings of laser l, they must have a number
            # greater than those with smaller l:
            acum = sum([detuning_indices[ll] for ll in range(l-1)])
            # We test the indices of all detunings of laser l to find the ones
            # we need.
            for kk in range(detuning_indices[l-1]):
                if detuningsij[l-1][kk][0]+1 == I_nd(i) and \
                   detuningsij[l-1][kk][1]+1 == I_nd(k):
                    The += '+detuning('+str(acum+kk+1)+')'
                if detuningsij[l-1][kk][0]+1 == I_nd(j) and \
                   detuningsij[l-1][kk][1]+1 == I_nd(k):
                    The += '-detuning('+str(acum+kk+1)+')'
            return The
        elif band2:
            The = ''
            # We need to find which detunings are delta^l_k,j and delta^l_k,i.
            # Since they are detunings of laser l, they must have a number
            # greater than those with smaller l:
            acum = sum([detuning_indices[ll] for ll in range(l-1)])
            # We test the indices of all detunings of laser l to find the ones
            # we need.
            for kk in range(detuning_indices[l-1]):
                if detuningsij[l-1][kk][0]+1 == I_nd(k) and \
                   detuningsij[l-1][kk][1]+1 == I_nd(j):
                    The += '+detuning('+str(acum+kk+1)+')'
                if detuningsij[l-1][kk][0]+1 == I_nd(k) and \
                   detuningsij[l-1][kk][1]+1 == I_nd(i):
                    The += '-detuning('+str(acum+kk+1)+')'
            return The
        else:
            if verbose > 1:
                s = 'WARNING: Optical frequencies will be used '
                s += 'instead for -omega_'+str(i)+","+str(j)
                s += '=' + str(omega_rescaled[i-1][j-1]) + '\n'
                print(s)
            return format_double(-omega_rescaled[i-1][j-1])
    ###########################################################################
    # This part is about finding detunings that reduce to
    # (theta_j -theta_i -omega_ij)

    # We establish to which omegaij the detunings should reduce to
    omegaij = [0 for k in range(Nnd)]
    omegaij[I_nd(i)-1] = 1; omegaij[I_nd(j)-1] = -1

    # We search for a combination of detunings that reduces to omegaij
    # That is a linear combination of detunings with coefficients the
    # that reduces to omegaij. Here
    detuning_failure = True
    for comb in combinations:
        suma = [sum([the[l]*detunings[l][comb[l]][iii] for l in range(Nl)])
                for iii in range(Nnd)]
        if suma == omegaij:
            detuning_failure = False
            break

    # We stop if no combination reduces to the needed expression.
    if detuning_failure:
        if verbose > 1:
            s = 'We will see if it is possible to express Theta_'+str(i)
            s += ','+str(j)+'=theta_'+str(j)+'-theta_'+str(i)+'-omega_'
            s += str(i)+','+str(j)
            print(s)
            s = 'in terms of other indices a,b such that omega_ab=omega_ij'
            s += ' and transition i -> j is allowed by Lij.'
            print(s)

        band3 = False
        for a in range(Ne):
            for b in range(a):
                if omega_rescaled[a][b] == omega_rescaled[i-1][j-1] and \
                   Lij[a][b] != []:
                    band3 = True
                    break
            if band3: break
        if verbose > 1:
            print(band3, a, b, i, j)

        a = a+1; b = b+1
        if band3 and (i != a or j != b):
            if verbose > 1:
                print(omega_rescaled[i-1][j-1], omega_rescaled[a-1][b-1])
                print(the)
                print('This was possible for omega_'+str(a)+','+str(b))

            return Theta(a, b, theta, omega_rescaled, omega_min, detunings,
                         detuningsij, combinations, detuning_indices, Lij,
                         other_the=the, verbose=verbose, states=states)
        else:
            if verbose > 0:
                aux = 'WARNING: It was impossible to express frequencies',
                aux += 'the and atomic transition -omega_'
                aux += str(i)+" "+str(j)+'='
                aux += str(-omega_rescaled[i-1][j-1],)
                aux += 'in terms of detunings from the'
                aux += 'transitions given by Lij'
                print(aux)

            if verbose > 0: print('I Will use optical frequencies instead.')
            The = ''
            # We give the optical frequencies
            for l in range(Nl):
                a = the[l]
                if a == 1:
                    The += '+'+format_double(omega_min[l])
                    The += '+detuning_knob('+str(l+1)+')'
                elif a == -1:
                    The += '-('+format_double(omega_min[l])
                    The += '+detuning_knob('+str(l+1)+'))'
                elif a == 0:
                    The += ''
                elif a > 0:
                    The += '+'+str(a)+'*'+format_double(omega_min[l])
                    The += '+detuning_knob('+str(l+1)+')'
                else:
                    The += str(a)+'*('+format_double(omega_min[l])
                    The += '+detuning_knob('+str(l+1)+'))'

            # We substract omega_ij
            The += format_double(-omega_rescaled[i-1][j-1])
            if verbose > 1: print(The)
            if verbose > 1: print()
            return The

    # For each optical frequency in the, we write the corresponding detuning.
    # This way of assigining a global index ll to the detunings ammounts to
    #   ll=   number_of_previous_detunings
    #       + number_of_detuning_ordered_by_row_and_from_left_to_right_column
    The = ''
    acum = 0
    ###########################################################################
    for l in range(Nl):
        if the[l] == 1:
            The += '+detuning('+str(acum+comb[l]+1)+')'
        elif the[l] == -1:
            The += '-detuning('+str(acum+comb[l]+1)+')'
        elif the[l] == 0:
            pass
        elif the[l] > 0:
            The += '+'+str(the[l])+'*detuning('+str(acum+comb[l]+1)+')'
        else:
            The += str(the[l])+'*detuning('+str(acum+comb[l]+1)+')'
        acum += detuning_indices[l]

    if The[:1] == '+':
        The = The[1:]

    ###########################################################################
    return The


def dot_product(laserl, sign, r, i, j):
    """This function calculates the dot product epsilon^(l(+-)) . vec(r_ij)."""
    if sign == 1:
        dp = sum([laserl.Yp[1-p]*r[p+1][i-1][j-1]*(-1)**p
                  for p in range(-1, 2)])
    elif sign == -1:
        dp = sum([laserl.Ym[1-p]*r[p+1][i-1][j-1]*(-1)**p
                  for p in range(-1, 2)])

    if not sage_included:
        return complex(dp)
    else:
        return dp


def find_omega_min(omega, Nl, detuningsij, i_d, I_nd):
    r"""This function returns a list of length Nl containing the mininmal frequency
    that each laser excites.
    """
    omega_min = []
    omega_min_indices = []
    for l in range(Nl):
        omegas = sorted([(omega[i_d(p[0]+1)-1][i_d(p[1]+1)-1], p)
                         for p in detuningsij[l]])
        omega_min += [omegas[0][0]]
        omega_min_indices += [omegas[0][1]]
    return omega_min, omega_min_indices


def write_equations_code(path, name, laser, omega, gamma, r, Lij, states=None,
                         excluded_mu=[], verbose=1):
    r"""Write code for the equations."""
    Ne = len(omega[0])
    Nl = len(laser)
    N_excluded_mu = len(excluded_mu)

    if states is None: states = range(1, Ne+1)

    omega_rescaled = omega[:]

    # We determine whether it is possible to eliminate explicit time-dependance
    theta = find_phase_transformation(Ne, Nl, r, Lij)

    # We construct the correspondence i <-> I between degenerate and
    # non-degenerate indices.
    i_d, I_nd, Nnd = calculate_iI_correspondence(omega)
    # We find the non-degenerate detunings
    detunings, detuningsij = laser_detunings(Lij, Nl, i_d, I_nd, Nnd)
    aux = find_omega_min(omega_rescaled, Nl, detuningsij, i_d, I_nd)
    omega_min, omega_min_indices = aux

    detuning_indices = [len(detunings[i]) for i in range(Nl)]
    Nd = sum([len(detunings[l]) for l in range(Nl)])
    combinations = detuning_combinations(detuning_indices)

    ##########################################################################

    # We add the code to caculate the detunings for each laser.
    code0 = ''
    code0 += '	!We calculate the detunings.\n'
    # We find the minimal frequency corresponding to each laser.
    code0 += '	!The list of detunings has the following meaning:\n'
    conta = 0
    for ll in range(Nl):
        for kk in range(len(detuningsij[ll])):
            conta += 1
            ii, jj = detuningsij[ll][kk]

            if states is None:
                code0 += '	!detuning('+str(conta)+')= delta^'+str(ll+1)+'_'
                code0 += str(i_d(ii+1))+','+str(i_d(jj+1))+'\n'
            else:
                state_i = str(states[i_d(ii+1)-1])[5:]
                state_j = str(states[i_d(jj+1)-1])[5:]
                code0 += '	!detuning('+str(conta)+')= delta^'+str(ll+1)+'_'
                code0 += state_i+','+state_j+'\n'

    det_index = 1
    for l in range(Nl):
        # omega0 = omega_min[l]
        i_min, j_min = omega_min_indices[l]
        for ii, jj in detuningsij[l]:
            code0 += '	detuning('+str(det_index)+')='
            # code0+=format_double(omega0-omega_rescaled[ii][jj])+'+(detuning_knob('+str(l+1)+'))\n'
            code0 += 'detuning_knob('+str(l+1)+') '
            code0 += '-('+format_double(
                omega_rescaled[i_d(ii+1)-1][i_d(i_min+1)-1])+')'
            code0 += '-('+format_double(
                omega_rescaled[i_d(j_min+1)-1][i_d(jj+1)-1])+')\n'
            det_index += 1
    code0 += '\n'

    ###########################################################################
    # We add here the code for the equations
    ###########################################################################
    code = ''
    ###########################################################################
    # We need to check that the resulting matrix A doesn't have any
    # row or column made of zeroes. If there is such a row or column
    # it means that for that particular i(mu),j(mu) or i(nu),j(nu)
    # neither the hamiltonian part of the correspondin equation nor the
    # phase transformation part nor the decay part were nonzero
    # (at least for the given light polarizations).
    # In this cases it is necessary to re-calculate all of the problem
    # without the corresponding mu components of A.
    # Also the right hand side will be checked for non zero elements.
    # Now we start the check for these exceptions.

    row_check = [False for mu in range(Ne**2-1-N_excluded_mu)]
    col_check = [False for nu in range(Ne**2-1-N_excluded_mu)]
    rhs_check = [False for nu in range(Ne**2-1-N_excluded_mu)]
    ##########################################################################
    # We give the code to calculate the independent vector.
    code += '	!We calculate the independent vector.\n'
    for i in range(2, Ne+1-N_excluded_mu):
        for s in [1, -1]:
            nu = Mu(i, 1, s, Ne, excluded_mu)
            rhs_check[nu-1] = True

            for l in Lij[i-1][0]:
                dp = dot_product(laser[l-1], 1, r, i, 1)
                dp = s*part(dp, -s)
                if dp != 0:
                    code += '	B('+str(nu)+',1)=B('+str(nu)+',1) +E0('+str(l)
                    code += ')*('+format_double(dp)+')\n'

    code += '\n'
    code += '	B=B/2.0d0\n\n'  # +str(1/sqrt(2.0))+'d0\n\n'

    ###########################################################################
    # We give the code to calculate the equations for populations.
    code += '	!We calculate the equations for populations.\n'
    for i in range(2, Ne+1):
        mu = Mu(i, i, 1, Ne, excluded_mu)
        for k in range(1, Ne+1):
            if k < i:
                for l in Lij[k-1][i-1]:
                    dp1 = dot_product(laser[l-1], -1, r, k, i)
                    dp2 = dot_product(laser[l-1], 1, r, i, k)

                    real_coef = part(dp1-dp2, -1)
                    imag_coef = part(dp1+dp2, 1)

                    if real_coef != 0:
                        nu = Mu(i, k, 1, Ne, excluded_mu)
                        code += '	A('+str(mu)+','+str(nu)+')=A('+str(mu)+','
                        code += str(nu)+')'
                        code += '+E0('+str(l)+')*('+format_double(real_coef)
                        code += ')\n'
                        row_check[mu-1] = True; col_check[nu-1] = True

                    if imag_coef != 0:
                        nu = Mu(i, k, -1, Ne, excluded_mu)
                        code += '	A('+str(mu)+','+str(nu)+')=A('+str(mu)+','
                        code += str(nu)+')'
                        code += '+E0('+str(l)+')*('+format_double(imag_coef)
                        code += ')\n'
                        row_check[mu-1] = True; col_check[nu-1] = True

            if k > i:
                for l in Lij[k-1][i-1]:
                    dp1 = dot_product(laser[l-1], -1, r, i, k)
                    dp2 = dot_product(laser[l-1], 1, r, k, i)

                    real_coef = -part(dp1-dp2, -1)
                    imag_coef = -part(dp1+dp2, 1)

                    if real_coef != 0:
                        nu = Mu(k, i, 1, Ne, excluded_mu)
                        code += '	A('+str(mu)+','+str(nu)+')=A('+str(mu)+','
                        code += str(nu)+')'
                        code += '+E0('+str(l)+')*('+format_double(real_coef)
                        code += ')\n'
                        row_check[mu-1] = True; col_check[nu-1] = True

                    if imag_coef != 0:
                        nu = Mu(k, i, -1, Ne, excluded_mu)
                        code += '	A('+str(mu)+','+str(nu)+')=A('+str(mu)+','
                        code += str(nu)+')'
                        code += '+E0('+str(l)+')*('+format_double(imag_coef)
                        code += ')\n'
                        row_check[mu-1] = True; col_check[nu-1] = True

    code += '\n'
    ###########################################################################
    # We give the code to calculate the equations for coherences
    # given in equations with label "stationary-coherences"

    code += '	!The code to calculate the equations for coherences.\n'

    for i in range(2, Ne+1):
        for j in range(1, i):
            for s in [1, -1]:
                mu = Mu(i, j, s, Ne, excluded_mu)
                for k in range(1, Ne+1):
                    for l in Lij[k-1][j-1]:
                        if k < i:
                            if k < j:
                                dp = s*dot_product(laser[l-1], -1, r, k, j)
                                dp1 = part(dp, -s)
                                dp2 = part(s*dp, +s)
                                nu = Mu(i, k, 1, Ne, excluded_mu)
                                if dp1 != 0:
                                    code += '	A('+str(mu)+','+str(nu)+')=A('
                                    code += str(mu)+','+str(nu)+')'
                                    code += '+E0('+str(l)+')*('
                                    code += format_double(dp1)+')\n'
                                    row_check[mu-1] = True
                                    col_check[nu-1] = True
                                nu = Mu(i, k, -1, Ne, excluded_mu)
                                if dp2 != 0:
                                    code += '	A('+str(mu)+','+str(nu)+')=A('
                                    code += str(mu)+','+str(nu)+')'
                                    code += '+E0('+str(l)+')*('
                                    code += format_double(dp2)+')\n'
                                    row_check[mu-1] = True
                                    col_check[nu-1] = True
                            elif k > j:
                                dp = s*dot_product(laser[l-1], +1, r, k, j)
                                dp1 = part(dp, -s)
                                dp2 = part(s*dp, s)
                                nu = Mu(i, k, 1, Ne, excluded_mu)
                                if dp1 != 0:
                                    code += '	A('+str(mu)+','+str(nu)+')=A('
                                    code += str(mu)+','+str(nu)+')'
                                    code += '+E0('+str(l)+')*('
                                    code += format_double(dp1)+')\n'
                                    row_check[mu-1] = True
                                    col_check[nu-1] = True
                                nu = Mu(i, k, -1, Ne, excluded_mu)
                                if dp2 != 0:
                                    code += '	A('+str(mu)+','+str(nu)+')=A('
                                    code += str(mu)+','+str(nu)+')'
                                    code += '+E0('+str(l)+')*('
                                    code += format_double(dp2)+')\n'
                                    row_check[mu-1] = True
                                    col_check[nu-1] = True
                        elif k > i:
                            dp = s*dot_product(laser[l-1], 1, r, k, j)
                            dp1 = part(dp, -s)
                            dp2 = part(-s*dp, s)
                            nu = Mu(k, i, 1, Ne, excluded_mu)
                            if dp1 != 0:
                                code += '	A('+str(mu)+','+str(nu)
                                code += ')=A('+str(mu)+','+str(nu)+')'
                                code += '+E0('+str(l)+')*('
                                code += format_double(dp1)+')\n'
                                row_check[mu-1] = True; col_check[nu-1] = True
                            nu = Mu(k, i, -1, Ne, excluded_mu)
                            if dp2 != 0:
                                code += '	A('+str(mu)+','+str(nu)+')=A('
                                code += str(mu)+','+str(nu)+')'
                                code += '+E0('+str(l)+')*('+format_double(dp2)
                                code += ')\n'
                                row_check[mu-1] = True; col_check[nu-1] = True
                    for l in Lij[i-1][k-1]:
                        if k > j:
                            if k > i:
                                dp = -s*dot_product(laser[l-1], -1, r, i, k)
                                dp1 = part(dp, -s)
                                dp2 = part(s*dp, s)
                                nu = Mu(k, j, 1, Ne, excluded_mu)
                                if dp1 != 0:
                                    code += '	A('+str(mu)+','+str(nu)
                                    code += ')=A('+str(mu)+','+str(nu)+')'
                                    code += '+E0('+str(l)+')*('
                                    code += format_double(dp1)+')\n'
                                    row_check[mu-1] = True
                                    col_check[nu-1] = True
                                nu = Mu(k, j, -1, Ne, excluded_mu)
                                if dp2 != 0:
                                    code += '	A('+str(mu)+','+str(nu)
                                    code += ')=A('+str(mu)+','+str(nu)+')'
                                    code += '+E0('+str(l)+')*('
                                    code += format_double(dp2)+')\n'
                                    row_check[mu-1] = True
                                    col_check[nu-1] = True
                            elif k < i:
                                dp = -s*dot_product(laser[l-1], 1, r, i, k)
                                dp1 = part(dp, -s)
                                dp2 = part(s*dp, s)
                                nu = Mu(k, j, 1, Ne, excluded_mu)
                                if dp1 != 0:
                                    code += '	A('+str(mu)+','+str(nu)
                                    code += ')=A('+str(mu)+','+str(nu)+')'
                                    code += '+E0('+str(l)+')*('
                                    code += format_double(dp1)+')\n'
                                    row_check[mu-1] = True
                                    col_check[nu-1] = True
                                nu = Mu(k, j, -1, Ne, excluded_mu)
                                if dp2 != 0:
                                    code += '	A('+str(mu)+','+str(nu)
                                    code += ')=A('+str(mu)+','+str(nu)+')'
                                    code += '+E0('+str(l)+')*('
                                    code += format_double(dp2)+')\n'
                                    row_check[mu-1] = True
                                    col_check[nu-1] = True
                        elif k < j:
                            dp = -s*dot_product(laser[l-1], 1, r, i, k)
                            dp1 = part(dp, -s)
                            dp2 = part(-s*dp, s)
                            nu = Mu(j, k, 1, Ne, excluded_mu)
                            if dp1 != 0:
                                code += '	A('+str(mu)+','+str(nu)+')=A('
                                code += str(mu)+','+str(nu)+')'
                                code += '+E0('+str(l)+')*('
                                code += format_double(dp1)+')\n'
                                row_check[mu-1] = True; col_check[nu-1] = True
                            nu = Mu(j, k, -1, Ne, excluded_mu)
                            if dp2 != 0:
                                code += '	A('+str(mu)+','+str(nu)+')=A('
                                code += str(mu)+','+str(nu)+')'
                                code += '+E0('+str(l)+')*('+format_double(dp2)
                                code += ')\n'
                                row_check[mu-1] = True; col_check[nu-1] = True
                for l in Lij[i-1][j-1]:
                    dp = s*part(dot_product(laser[l-1], 1, r, i, j), -s)
                    nu = Mu(i, i, 1, Ne, excluded_mu)
                    if dp != 0:
                        code += '	A('+str(mu)+','+str(nu)+')=A('
                        code += str(mu)+','+str(nu)+')'
                        code += '+E0('+str(l)+')*('+format_double(dp)+')\n'
                        row_check[mu-1] = True; col_check[nu-1] = True
                        nu = Mu(j, j, 1, Ne, excluded_mu)
                        if nu == 0:
                            for n in range(1, Ne):
                                code += '	A('+str(mu)+','+str(n)+')=A('
                                code += str(mu)+','+str(n)+')'
                                code += '+E0('+str(l)+')*('+format_double(dp)
                                code += ')\n'
                                row_check[mu-1] = True; col_check[n-1] = True
                        else:
                            code += '	A('+str(mu)+','+str(nu)+')=A('
                            code += str(mu)+','+str(nu)+')'
                            code += '+E0('+str(l)+')*('+format_double(-dp)
                            code += ')\n'
                            row_check[mu-1] = True; col_check[nu-1] = True

    code += '\n'
    code += '	A=A/2.0d0\n\n'  # +str(1/sqrt(2.0))+'d0\n\n'
    ###########################################################################
    # We add the terms associated with the phase transformation.
    code += '	!We calculate the terms associated with the phase '
    code += 'transformation.\n'

    # for i in range(2,10):
    for i in range(2, Ne+1):
        for j in range(1, i):
            extra = Theta(i, j, theta, omega_rescaled, omega_min, detunings,
                          detuningsij, combinations, detuning_indices, Lij,
                          i_d, I_nd, Nnd, verbose=verbose, states=states)

            if extra != '':
                for s in [1, -1]:
                    mu = Mu(i, j, s, Ne, excluded_mu)
                    nu = Mu(i, j, -s, Ne, excluded_mu)

                    code += '    A('+str(mu)+','+str(nu)+')=A('+str(mu)+','
                    code += str(nu)+')'
                    if s == 1:
                        code += '-('+str(extra)+')\n'
                    elif s == -1:
                        code += '+('+str(extra)+')\n'

                    row_check[mu-1] = True; col_check[nu-1] = True

    code += '\n'
    ###########################################################################
    # We add the terms associated with spontaneous decay.
    code += '	!We calculate the terms associated with spontaneous decay.\n'
    # First for populations.
    for i in range(2, Ne+1):
        mu = Mu(i, i, 1, Ne, excluded_mu)
        for k in range(1, Ne+1):
            gams = 0
            if k < i:
                gams += gamma[i-1][k-1]
            elif k > i:
                nu = Mu(k, k, 1, Ne, excluded_mu)
                ga = gamma[i-1][k-1]
                if ga != 0:
                    code += '    A('+str(mu)+','+str(nu)+')=A('+str(mu)+','
                    code += str(nu)+')'
                    code += '-('+format_double(ga)+')\n'
                    row_check[mu-1] = True; col_check[nu-1] = True
            if gams != 0:
                code += '    A('+str(mu)+','+str(mu)+')=A('+str(mu)+','
                code += str(mu)+')'
                code += '-('+format_double(gams)+')\n'
                row_check[mu-1] = True; col_check[mu-1] = True

    # And now for coherences
    for i in range(1, Ne+1):
        for j in range(1, i):
            gams = gamma[i-1][j-1]/2
            if gams != 0:
                for a in range(i+1, Ne+1):
                    mu = Mu(a, i, 1, Ne, excluded_mu)
                    code += '    A('+str(mu)+','+str(mu)+')=A('+str(mu)+','
                    code += str(mu)+')'
                    code += '-('+format_double(gams)+')\n'
                    row_check[mu-1] = True; col_check[mu-1] = True
                    mu = Mu(a, i, -1, Ne, excluded_mu)
                    code += '    A('+str(mu)+','+str(mu)+')=A('+str(mu)+','
                    code += str(mu)+')'
                    code += '-('+format_double(gams)+')\n'
                    row_check[mu-1] = True; col_check[mu-1] = True

                for b in range(1, i):
                    mu = Mu(i, b, 1, Ne, excluded_mu)
                    code += '    A('+str(mu)+','+str(mu)+')=A('+str(mu)+','
                    code += str(mu)+')'
                    code += '-('+format_double(gams)+')\n'
                    row_check[mu-1] = True; col_check[mu-1] = True
                    mu = Mu(i, b, -1, Ne, excluded_mu)
                    code += '    A('+str(mu)+','+str(mu)+')=A('+str(mu)+','
                    code += str(mu)+')'
                    code += '-('+format_double(gams)+')\n'
                    row_check[mu-1] = True; col_check[mu-1] = True

    code += '\n'
    code = code0+code
    Nd = sum([len(detunings[l]) for l in range(Nl)])
    res = [code, Nd, row_check, col_check, rhs_check, Ne, N_excluded_mu,
           states, omega_min, detuningsij, omega_rescaled]
    return res


def compile_code(path, name, optimization_flag=' -O3',
                 lapack=False, parallel=True, clone=None, verbose=0):
    r"""Compile fortran code."""
    from config import use_netcdf
    t0 = time()
    parallel_flag = ''; end_flags = ''
    if lapack: end_flags += ' -llapack'
    if use_netcdf: end_flags += ' -lnetcdff -lnetcdf'
    if parallel: parallel_flag += ' -fopenmp'

    # We establish the name of the clone.
    if clone is not None:
        clone = '_'+str(clone)
    else:
        clone = ''

    # We read the code and add the "clone" to every file name.
    f = file(path+name+'.f90', 'r')
    code = f.read()
    f.close()
    code = code.replace('.dat', clone+'.dat')
    code = code.replace('.nc', clone+'.nc')

    # We save the code in a clone file.
    f = file(path+name+clone+'.f90', 'w')
    f.write(code)
    f.close()

    from fast.config import fast_path
    # com='gfortran -I '+fast_path+' '+parallel_flag+optimization_flag+
    # ' '+path+name+clone+'.f90 -o '+path+name+clone+end_flags
    com = 'gfortran -I '
    com += fast_path+" "
    com += parallel_flag
    com += optimization_flag+' '
    com += path+name+clone+'.f90 -o '+path+name+clone
    com += end_flags
    if verbose > 0:
        print(com)
    exit_code = os.system(com)
    if exit_code != 0:
        s = 'command: '+com+' returned exit_code '+str(exit_code)
        raise RuntimeError(s)

    if clone != '':
        os.system('rm '+path+name+clone+'.f90')
    return time()-t0


def convolve_with_gaussian(x, f, sigma):
    r"""Return the convolution of an array with a gaussian of width sigma."""
    # We will calculate with data from a zero-centered normalized gaussian
    # distribution such that the steps between data are the same as the
    # original signal.
    N = len(f); step = (x[1]-x[0])
    a = (x[-1]-x[0])/2.0

    x_gaussian = [-a + i*step for i in range(N)]
    gaussian = [0.398942280401433*exp(-xi**2/(2*sigma**2))/sigma
                for xi in x_gaussian]

    # We calculate N points of the convolution.
    fg = np.convolve(f, gaussian, mode='same')

    # We correct for border effects.
    stepcor = step*float(N)/float(N+1)
    xfg = [x[0] + i*stepcor for i in range(len(fg))]

    # We give the convolution it's correct units.
    fg = [fg[i]*step for i in range(len(fg))]
    return xfg, fg


def block_diagonal_matrix(matrices, type=None):
    """Build a block-diagonal matrix out of a given list of matrices.

    The type of matrix is chosen according to the input matrices.

    >>> import numpy as np
    >>> import sympy as sy
    >>> lis = [sy.ones(2), 2*sy.ones(3), 3*sy.ones(4)]
    >>> sy.pprint(block_diagonal_matrix(lis))
    ⎡1  1  0  0  0  0  0  0  0⎤
    ⎢                         ⎥
    ⎢1  1  0  0  0  0  0  0  0⎥
    ⎢                         ⎥
    ⎢0  0  2  2  2  0  0  0  0⎥
    ⎢                         ⎥
    ⎢0  0  2  2  2  0  0  0  0⎥
    ⎢                         ⎥
    ⎢0  0  2  2  2  0  0  0  0⎥
    ⎢                         ⎥
    ⎢0  0  0  0  0  3  3  3  3⎥
    ⎢                         ⎥
    ⎢0  0  0  0  0  3  3  3  3⎥
    ⎢                         ⎥
    ⎢0  0  0  0  0  3  3  3  3⎥
    ⎢                         ⎥
    ⎣0  0  0  0  0  3  3  3  3⎦

    >>> lis = [np.ones((2, 2)), 2*np.ones((3, 3)), 3*np.ones((4, 4))]
    >>> print(block_diagonal_matrix(lis))
    [[1. 1. 0. 0. 0. 0. 0. 0. 0.]
     [1. 1. 0. 0. 0. 0. 0. 0. 0.]
     [0. 0. 2. 2. 2. 0. 0. 0. 0.]
     [0. 0. 2. 2. 2. 0. 0. 0. 0.]
     [0. 0. 2. 2. 2. 0. 0. 0. 0.]
     [0. 0. 0. 0. 0. 3. 3. 3. 3.]
     [0. 0. 0. 0. 0. 3. 3. 3. 3.]
     [0. 0. 0. 0. 0. 3. 3. 3. 3.]
     [0. 0. 0. 0. 0. 3. 3. 3. 3.]]

    """
    if type is None:
        type = np.float64
    sizes = [Ai.shape[0] for Ai in matrices]
    size = sum(sizes)
    symbolic = hasattr(matrices[0][0], "subs")
    if symbolic:
        A = symzeros(size, size)
    else:
        A = np.zeros((size, size), type)
    ini = 0; fin = 0
    for i, sizei in enumerate(sizes):
        fin += sizei
        A[ini: fin, ini: fin] = matrices[i]
        ini += sizei

    return A


if __name__ == "__main__":
    import doctest
    print(doctest.testmod(verbose=False))
