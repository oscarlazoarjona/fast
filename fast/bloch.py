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
r"""This file will contain the new write_bloch function."""

from sympy import Symbol, diff
from fast.symbolic import cartesian_dot_product, define_frequencies
from fast.symbolic import define_laser_variables
from scipy.constants import physical_constants

import numpy as np
import fast
import sympy

from numpy import array as nparray
from numpy import sqrt as npsqrt
from sympy import Matrix, sqrt, I

hbar_num = physical_constants["Planck constant over 2 pi"][0]
e_num = physical_constants["elementary charge"][0]
a0 = physical_constants["Bohr radius"][0]


def phase_transformation(Ne, Nl, rm, xi, return_equations=False):
    """This function returns a phase transformation specified as a list of
    lenght Ne whose elements correspond to each theta_i. Each element is a list
    of length Nl which specifies the coefficients multiplying each optical
    frequency omega^l. So for instance [[1,1],[1,0],[0,0]] means

    theta_1=omega^1+omega^2
    theta_2=omega^1
    theta_3=0.
    """
    # We first define the needed variables
    E0, omega_laser = define_laser_variables(Nl)
    theta = [Symbol('theta'+str(i+1)) for i in range(Ne)]

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
    sol = sympy.solve(eqs, theta, dict=True)[0]
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
    r"""Find the smallest transition frequency for each field."""
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
    r"""Get the indices of the detunings of all fields.

    They are returned in the form
        [[(i1, j1), (i2, j2)], ...,[(i1, j1)]].
    One list of pairs of indices for each field.
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
    r"""Get the code to calculate the simplified detunings."""
    code_det = ""

    for l in range(Nl):
        for pair in pairs[l]:
            iu, ju = pair
            code_det += "    delta"+str(l+1)
            code_det += "_"+str(iu+1)
            code_det += "_"+str(ju+1)
            code_det += " = detuning_knob["+str(l)+"]"
            corr = -omega_levelu[iu]+omega_levelu[iu0[l]]
            corr = omega_levelu[ju0[l]]-omega_levelu[ju]
            if corr != 0:
                code_det += " + ("+str(corr)+")"
            code_det += "\n"
    return code_det


def fast_hamiltonian(Ep, epsilonp, detuning_knob, rm, omega_level, xi, theta,
                     file_name=None):
    r"""Return a fast function that returns a Hamiltonian as a numerical array.
    lowest frequency transition. They can either be symbolic expressions or

    All quantities should be in SI units.

    The arguments Ep, epsilonp, and detuning_knob represent the electric field
    amplitudes, field polarizations, and the detunings of each field from the
    numerical values.

    Eventually, the way this will work is that the returned function will take
    these Ep, epsilonp, and detuning_knob as numerical arguments if symbolic
    expressions were used, or use the given numerical values by default. At the
    moment however, only variable detuning_knob is supported.

    The argument rm should be numerical values of the below-diagonal components
    of the position operator in the cartesian basis:
        rm = [ x_ij, y_ij, z_ij ] for 0 < j < i

    The argument omega_level should be numerical values of the energy levels.
    theta should be a phase transformation returned by the phase_transformation
    function. xi should be an array of ones and zeros such that xi[l, i, j]
    represents whether the |i> -> |j> transition is driven by field l.

    """
    if True:
        # We find out the number of fields and states.
        Nl = len(Ep)
        Ne = np.array(rm[0]).shape[0]
        #######################################################################
        # We determine which arguments are constants.
        try:
            Ep = np.array([complex(Ep[l]) for l in range(Nl)])
            constant_Ep = True
        except:
            constant_Ep = False

        try:
            epsilonp = [np.array([complex(epsilonp[l][i]) for i in range(3)])
                        for l in range(Nl)]
            constant_epsilonp = True
        except:
            constant_epsilonp = False

        try:
            detuning_knob = np.array([float(detuning_knob[l])
                                      for l in range(Nl)])
            constant_detuning_knob = True
        except:
            constant_detuning_knob = False

        if not constant_Ep:
            raise NotImplementedError("Ep must be a constant.")
        if not constant_epsilonp:
            raise NotImplementedError("epsilonp must be constant.")
        if constant_detuning_knob:
            raise NotImplementedError("detuning_knob must be variable.")

        # We convert rm to a numpy array
        rm = np.array([[[complex(rm[k][i, j])
                       for j in range(Ne)] for i in range(Ne)]
                       for k in range(3)])

        #######################################################################
        # We establish the arguments.
        code = ""
        code += "def hamiltonian("
        if not constant_Ep: code += "Ep, "
        if not constant_epsilonp: code += "epsilonp, "
        if not constant_detuning_knob: code += "detuning_knob, "
        if code[-2:] == ", ":
            code = code[:-2] + "):\n"

        code += "    H = np.zeros(("+str(Ne)+", "+str(Ne)+"), complex)\n\n"
        #######################################################################
        # We get the code for the below-diagonal elements.
        code += "    # We calculate the below-diagonal elements.\n"

        for i in range(Ne):
            for j in range(i):
                for l in range(Nl):
                    if xi[l, i, j] == 1.0:
                        # We get the below-diagonal terms.
                        code += "    H["+str(i)+", "+str(j)+"] = "
                        # We get the code for Ep.
                        if constant_Ep:
                            code += str(Ep[l])
                        else:
                            code += "Ep["+str(l)+"]"
                        # We get the code for epsilonp dot rm
                        if constant_epsilonp:
                            rmij = rm[:, i, j]
                            dp = cartesian_dot_product(epsilonp[l], rmij)
                            code += "*("+str(dp)+")"
                        else:
                            code += "cartesian_dot_product(epsilonp[l], rm)"

                        code += "\n"

        #######################################################################
        # We get the code for the above-diagonal elements.
        code += r"""
    # We calculate the above-diagonal elements.
    for i in range("""+str(Ne)+"""):
        for j in range(i+1, """+str(Ne)+"""):
            H[i, j] = H[j, i].conjugate()\n\n"""
    ###########################################################################
    # We get the code for the diagonal elements.
    code += "    # We calculate the diagonal elements.\n"
    # 1 We build the degeneration simplification and its inverse (to avoid
    # large combinatorics).
    aux = define_simplification(omega_level, xi, Nl)
    u, invu, omega_levelu, Neu, xiu = aux
    # For each field we find the smallest transition frequency, and its
    # simplified indices.
    omega_min, iu0, ju0 = find_omega_min(omega_levelu, Neu, Nl, xiu)
    # We get the code to calculate the non degenerate detunings.
    pairs = detunings_indices(Neu, Nl, xiu)
    code_det = detunings_code(Neu, Nl, pairs, omega_levelu, iu0, ju0)
    code += code_det
    #####################################
    # We find the coefficients a_l that multiply omega_laser_l in
    # H_ii = omega_level_iu + theta_iu = \sum_i a_i varpi_i + other terms
    _omega_levelu, omega, gamma = define_frequencies(Neu)
    E0, omega_laser = define_laser_variables(Nl)
    # _omega_level, omega, gamma = define_frequencies(Ne)
    print pairs
    combs = [[pairs[0][k]] for k in range(len(pairs[0]))]
    print combs

    for i in range(Ne):
        _Hii = theta[i] + _omega_levelu[u(i)]
        a = [diff(_Hii, omega_laser[l]) for l in range(Nl)]
        # We look for a combination of the detunings obtained with the function
        # detunings_code. So we build all combinations.
        a = str(a)

    # We define symbolic omegas.
    #####################################
    # 5 We get the code to calculate the non degenerate detunings.

    # print omega_min, iu0, ju0

    # print xi
    # print xiu
    #
    # print Neu
    # print omega_level
    # print omega_levelu
    # for iu in range(Neu):
    #     print iu, invu(iu)
    #####################################
    # We get the code.
    # code += "    # We calculate the diagonal elements.\n"
    # for i in range(Ne):
    #     if not constant_detuning_knob:
    #         code += "    H["+str(i)+", "+str(i)+"] = "
    #         code += str(theta[i]+_omega_level[i])+"\n"

    ###########################################################################
    # print Ne, Nl
    code += "\n"
    code += "    return H\n"

    if file_name is not None:
        f = file(file_name, "w")
        f.write(code)
        f.close()

    hamiltonian = 0
    # print code
    exec code
    return hamiltonian


###############################################################################
# Ne = 2; Nl = 1
# E0, omega_laser = fast.symbolic.define_laser_variables(Nl)
# E0 = [np.sqrt(2.0)]
#
# epsilonp = fast.symbolic.polarization_vector(0, 0, 0, sympy.pi/8, 1)
#
# detuning_knob = [sympy.symbols("delta1")]
#
# # sympy.pprint(E0)
# # sympy.pprint(epsilonp)
# # sympy.pprint(detuning_knob)
#
# rm = fast.symbolic.define_r_components(Ne, explicitly_hermitian=True, p=-1)
# rm[0][1, 0] = 1.0
# rm[1][1, 0] = 1.0
# rm[2][1, 0] = 1.0
# # sympy.pprint(rm)
#
# omega_level = [100.0*i for i in range(Ne)]
# # print omega_level
# xi = np.zeros((Nl, Ne, Ne))
# xi[0, 1, 0] = 1.0
# xi[0, 0, 1] = 1.0
#
# phase = phase_transformation(Ne, Nl, rm, xi, return_equations=False)
#
# hamiltonian2 = fast_hamiltonian(E0, [epsilonp], detuning_knob, rm,
#                                 omega_level, xi, phase, "code.py")
#
# detuning_knob = [1.0]
# print hamiltonian2(detuning_knob)
###############################################################################
# A real example.
element = "Rb"; isotope = 87; N = 5
fine_states = [fast.State(element, isotope, N, 0, 1/fast.Integer(2)),
               fast.State(element, isotope, N, 1, 3/fast.Integer(2))]
magnetic_states = fast.make_list_of_states(fine_states, "magnetic")

Ne = len(magnetic_states)
Nl = 1
E0 = [1.0]
epsilonp = [[0.0, 0.0, 1.0]]
detuning_knob = [0.0]


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

    """
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


r = fast.calculate_matrices(magnetic_states)
# r = [np.array(r[i]) for i in range(3)]
print r[0][0][0], type(r[0][0][0])
r = helicity_to_cartesian(r, numeric=True)
# print len(r), r[0].shape
print r[0].shape
# print r[0][0, 0], type(r[0][0, 0])
print

# omega_level = [1.0, 1.0, 1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 4.0]
# print omega_level
# simp = define_simplification(omega_level)
#
# for i in range(len(omega_level)):
#     print i, simp(i)
