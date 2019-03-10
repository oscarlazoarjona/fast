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

r"""This is a template."""
import numpy as np
from sympy import symbols
from scipy.constants import physical_constants

from fast.bloch import (Unfolding, fast_bloch_equations,
                        detunings_code, define_simplification, find_omega_min,
                        detunings_indices, detunings_combinations,
                        detunings_rewrite)
from fast.symbolic import define_frequencies, define_laser_variables


c_num = physical_constants["speed of light in vacuum"][0]
k_B_num = physical_constants["Boltzmann constant"][0]


class Inhomogeneity(object):
    r"""A class describing an ensemble of atoms driven different dynamics."""

    def __init__(self, domain, distribution, terms):
        r"""A class describing an ensemble of atoms driven different dynamics.

        To describe Doppler broadening from a Maxwell-Boltzmann velocity
        distribution.

        >>> from fast.atomic_structure import speed_average
        >>> from fast import PlaneWave, Atom
        >>> from fast.bloch import phase_transformation
        >>> import numpy as np
        >>> from sympy import Matrix

        >>> element = "Rb"
        >>> isotope = 85
        >>> T = 273.15+20
        >>> v_av = speed_average(T, element, isotope)
        >>> v = np.linspace(-4*v_av, 4*v_av, 11)
        >>> mass = Atom(element, isotope).mass
        >>> f = fast_maxwell_boltzmann(mass)
        >>> distribution = f([v], T)
        >>> print distribution
        [  3.34624364e-12   5.12403858e-09   1.53771999e-06   9.04382537e-05
           1.04240650e-03   2.35468108e-03   1.04240650e-03   9.04382537e-05
           1.53771999e-06   5.12403858e-09   3.34624364e-12]

        We obtain a fast function for the Doppler terms.

        >>> Ne = 2
        >>> Nl = 1
        >>> unfolding = Unfolding(Ne, True, True, True)

        >>> a0 = physical_constants["Bohr radius"][0]
        >>> rm = [Matrix([[0, 0], [a0, 0]]),
        ...       Matrix([[0, 0], [0, 0]]),
        ...       Matrix([[0, 0], [0, 0]])]
        >>> xi = np.array([[[0, 1], [1, 0]]])
        >>> theta = phase_transformation(Ne, Nl, rm, xi)

        >>> omega_level = [1, 2.4e15]

        >>> v_sym = symbols("vx vy vz")
        >>> detuning_knob = [symbols("delta1")]
        >>> laser = PlaneWave(0, 0, 0, 0, 1)
        >>> k = [laser.k]

        >>> doppler_terms = fast_doppler_terms(v_sym, detuning_knob, k,
        ...                                    omega_level,
        ...                                    xi, theta, unfolding,
        ...                                    matrix_form=True)

        >>> doppler_effect = Inhomogeneity([v], distribution, doppler_terms)

        """
        self.domain = domain
        self.distribution = distribution
        self.terms = terms
        self.shape = distribution.shape

    def __repr__(self):
        r"""Return a string representation."""
        return "An inhomogeneity of shape "+str(self.shape)

    def average(self, rho):
        r"""Return the average density matrix of an inhomogeneous ensemble."""
        def marginal(f, rho):
            remaining = len(f.shape)
            if remaining == 0:
                return rho
            rho = sum([f[i]*rho[i] for i in range(rho.shape[0])])
            f = np.sum(f, 0)
            return marginal(f, rho)

        return marginal(self.distribution, rho)


class DopplerBroadening(Inhomogeneity):
    r"""An object representing an ensemble of atom at different velocities."""

    def __init__(self, shape, stds, T, mass, detuning_knob, k,
                 omega_level, xi, theta, unfolding, axes=["x", "y", "z"],
                 matrix_form=False, file_name=None, return_code=False):
        r"""An object representing Doppler broadening.

        >>> from fast.atomic_structure import speed_average
        >>> from fast import PlaneWave, Atom
        >>> from fast.bloch import phase_transformation
        >>> import numpy as np
        >>> from sympy import Matrix

        >>> element = "Rb"
        >>> isotope = 85
        >>> T = 273.15+20
        >>> v_av = speed_average(T, element, isotope)
        >>> v = np.linspace(-4*v_av, 4*v_av, 11)
        >>> mass = Atom(element, isotope).mass
        >>> f = fast_maxwell_boltzmann(mass)
        >>> distribution = f(v, T)
        >>> print distribution
        [  3.34624364e-12   5.12403858e-09   1.53771999e-06   9.04382537e-05
           1.04240650e-03   2.35468108e-03   1.04240650e-03   9.04382537e-05
           1.53771999e-06   5.12403858e-09   3.34624364e-12]

        We obtain a fast function for the Doppler terms.

        >>> Ne = 2
        >>> Nl = 1
        >>> unfolding = Unfolding(Ne, True, True, True)

        >>> a0 = physical_constants["Bohr radius"][0]
        >>> rm = [Matrix([[0, 0], [a0, 0]]),
        ...       Matrix([[0, 0], [0, 0]]),
        ...       Matrix([[0, 0], [0, 0]])]
        >>> xi = np.array([[[0, 1], [1, 0]]])
        >>> theta = phase_transformation(Ne, Nl, rm, xi)

        >>> omega_level = [1, 2.4e15]

        >>> v_sym = symbols("vx vy vz")
        >>> detuning_knob = [symbols("delta1")]
        >>> laser = PlaneWave(0, 0, 0, 0, 1)
        >>> k = [laser.k]

        >>> shape = [11]
        >>> stds = [[-4, 4]]

        >>> doppler_effect = DopplerBroadening(shape, stds, T, mass,
        ...                                    detuning_knob, k,
        ...                                    omega_level, xi, theta,
        ...                                    unfolding)
        >>> print doppler_effect.domain
        [array([-677.70074607, -542.16059685, -406.62044764, -271.08029843,
               -135.54014921,    0.        ,  135.54014921,  271.08029843,
                406.62044764,  542.16059685,  677.70074607])]

        >>> print doppler_effect.distribution
        [  1.07064870e-04   1.90728284e-03   1.79157396e-02   8.87372390e-02
           2.31754734e-01   3.19155879e-01   2.31754734e-01   8.87372390e-02
           1.79157396e-02   1.90728284e-03   1.07064870e-04]
        >>> print sum(doppler_effect.distribution)
        1.0

        >>> print doppler_effect.v_sig
        169.425186516

        """
        ######################################################################
        # We obtain the domain of the velocity distribution.
        if hasattr(shape, "__getitem__"):
            dimension = len(shape)
        else:
            dimension = 1
        v_symb = symbols("v1:"+str(dimension+1))
        v_sig = np.sqrt(k_B_num*T/mass)

        domain = [np.linspace(stds[i][0]*v_sig, stds[i][1]*v_sig, shape[i])
                  for i in range(dimension)]
        domain = np.meshgrid(*domain)
        if dimension in [2, 3]:
            for i in range(dimension):
                domain[i] = np.swapaxes(domain[i], 0, 1)

        ######################################################################
        # We obtain a function to calculate the distribution.
        f = fast_maxwell_boltzmann(mass)
        distribution = f(domain, T)
        # We renormalize from probability density to probability.
        not_one = sum(distribution.flatten())
        distribution = distribution/not_one

        ######################################################################
        doppler_terms = fast_doppler_terms(v_symb, detuning_knob, k,
                                           omega_level, xi, theta, unfolding,
                                           axes=axes,
                                           matrix_form=matrix_form,
                                           file_name=file_name)

        Inhomogeneity.__init__(self, domain, distribution, doppler_terms)
        self.T = T
        self.v_sig = v_sig
        self.stds = stds
        self.mass = mass
        self.detuning_knob = detuning_knob
        self.k = k
        self.omega_level = omega_level
        self.xi = xi
        self.theta = theta
        self.unfolding = unfolding
        self.matrix_form = matrix_form
        self.axes = axes

    def reset(self, T):
        r"""Recalculate the doppler broadening for a given temperature."""
        self.__init__(self.shape, self.stds, T,
                      self.mass, self.detuning_knob, self.k,
                      self.omega_level, self.xi, self.theta, self.unfolding,
                      self.axes,
                      self.matrix_form)


def fast_doppler_terms(v, detuning_knob, k, omega_level, xi, theta,
                       unfolding, axes=["x", "y", "z"], matrix_form=False,
                       file_name=None,
                       return_code=False):
    r"""Return a fast function that returns the Doppler terms.

    >>> from sympy import Matrix, symbols
    >>> from scipy.constants import physical_constants
    >>> from fast import PlaneWave
    >>> from fast.bloch import phase_transformation
    >>> Ne = 2
    >>> Nl = 1
    >>> unfolding = Unfolding(Ne, True, True, True)

    >>> a0 = physical_constants["Bohr radius"][0]
    >>> rm = [Matrix([[0, 0], [a0, 0]]),
    ...       Matrix([[0, 0], [0, 0]]),
    ...       Matrix([[0, 0], [0, 0]])]
    >>> xi = np.array([[[0, 1], [1, 0]]])
    >>> theta = phase_transformation(Ne, Nl, rm, xi)

    >>> omega_level = [1, 2.4e15]
    >>> detuning_knob = [symbols("delta1")]

    >>> v = symbols("vx vy vz")
    >>> laser = PlaneWave(0, 0, 0, 0, 1)
    >>> k = [laser.k]

    >>> doppler_terms = fast_doppler_terms(v, detuning_knob, k, omega_level,
    ...                                    xi, theta, unfolding,
    ...                                    matrix_form=True)

    >>> detuning_knobs = [0]
    >>> A, b = doppler_terms([0, 0, 167], detuning_knobs)
    >>> print A/2/np.pi*1e-9
    [[ 0.          0.          0.        ]
     [ 0.          0.          0.21277821]
     [ 0.         -0.21277821  0.        ]]

    """
    # We unpack variables.
    if True:
        Ne = unfolding.Ne
        Nrho = unfolding.Nrho
        Nl = xi.shape[0]
        IJ = unfolding.IJ
        Mu = unfolding.Mu
    # We determine which arguments are constants.
    if True:
        try:
            detuning_knob = np.array([float(detuning_knob[l])
                                      for l in range(Nl)])
            variable_detuning_knob = False
        except:
            variable_detuning_knob = True
    # We establish the arguments of the output function.
    if True:
        code = ""
        code += "def doppler_terms(v, "
        if not matrix_form: code += "rho, "
        if variable_detuning_knob: code += "detuning_knob, "
        if code[-2:] == ", ":
            code = code[:-2]
        code += "):\n"

        code += '    r"""A fast calculation of the Doppler terms."""\n'
    # We initialize the output and auxiliaries.
    if True:
        # We introduce the factor that multiplies all terms.
        if unfolding.real:
            code += "    fact = 1.0\n\n"
        else:
            code += "    fact = 1.0j\n\n"

        if matrix_form:
            code += "    A = np.zeros(("+str(Nrho)+", "+str(Nrho)
            if not unfolding.real:
                code += "), complex)\n\n"
            else:
                code += "))\n\n"
            if unfolding.normalized:
                code += "    b = np.zeros(("+str(Nrho)
                if not unfolding.real:
                    code += "), complex)\n\n"
                else:
                    code += "))\n\n"
        else:
            code += "    rhs = np.zeros(("+str(Nrho)
            if not unfolding.real:
                code += "), complex)\n\n"
            else:
                code += "))\n\n"

    # We build the degeneration simplification and is inverse (to avoid
    # large combinatorics).
    aux = define_simplification(omega_level, xi, Nl)
    u, invu, omega_levelu, Neu, xiu = aux
    # For each field we find the smallest transition frequency, and its
    # simplified indices.
    omega_min, iu0, ju0 = find_omega_min(omega_levelu, Neu, Nl, xiu)
    #####################################
    # We get the code to calculate the non degenerate detunings.
    pairs = detunings_indices(Neu, Nl, xiu)
    if not variable_detuning_knob:
        code += "    detuning_knob = np.zeros("+str(Nl)+")\n"
        for l in range(Nl):
            code += "    detuning_knob["+str(l)+"] = " +\
                str(detuning_knob[l])+"\n"
    # We put in the code to calculate the Doppler shift
    if True:
        # Delta omega = omega_laser /c k.v
        code += "    # The Doppler shift\n"
        dimension = len(v)
        code += "    c = %s\n" % c_num
        for l in range(Nl):
            lp1 = l+1
            code += "    omega_laser%s = %s" % (lp1, omega_min[l])
            code += "+detuning_knob[%s]\n" % l
            code += "    detuning_knob[%s] = 0\n" % l

            axeindices = {"x": 0, "y": 1, "z": 2}
            axeindices = [axeindices[axe] for axe in axes]
            for ii in axeindices[:dimension]:
                code += "    detuning_knob[%s] += -%s" % (l, k[l][ii])
                if dimension == 1:
                    code += "*v/c*omega_laser%s\n" % lp1
                else:
                    code += "*v[%s]/c*omega_laser%s\n" % (ii, lp1)

    code += "\n"

    code_det = detunings_code(Neu, Nl, pairs, omega_levelu, iu0, ju0)
    code += code_det
    code += "\n"
    #####################################
    # There is a term
    # I * Theta_ij * rho_ij = I * (omega_level_j - omega_level_i
    #                              theta_j - theta_i)
    # for all i != j.
    # This term can be re expressed as
    # re(Theta_ij*rho_ij) = - Theta_ij * im(rho_ij)
    # im(Theta_ij*rho_ij) = + Theta_ij * re(rho_ij)
    _omega_level, omega, gamma = define_frequencies(Ne)
    _omega_levelu, omega, gamma = define_frequencies(Neu)
    E0, omega_laser = define_laser_variables(Nl)
    # We build all combinations.
    combs = detunings_combinations(pairs)
    # We add all terms.
    for mu in range(Nrho):
        s, i, j = IJ(mu)
        if i != j:
            _Thetaij = _omega_levelu[u(j)] - _omega_levelu[u(i)]
            _Thetaij += theta[j] - theta[i]
            aux = (_Thetaij, combs, omega_laser,
                   _omega_levelu, omega_levelu, iu0, ju0)
            assign = detunings_rewrite(*aux)
            if assign != "":
                if s == 0:
                    nu = mu
                elif s == 1:
                    assign = "-(%s)" % assign
                    nu = Mu(-s, i, j)
                elif s == -1:
                    assign = "+(%s)" % assign
                    nu = Mu(-s, i, j)

                if matrix_form:
                    term_code = "    A[%s, %s] = %s\n" % (mu, nu, assign)
                else:
                    term_code = "    rhs[%s] = (%s)*rho[%s]\n" % (mu,
                                                                  assign, nu)
            else:
                term_code = ""
            code += term_code
    #####################################
    # We finish the code.
    if True:
        if matrix_form:
            if unfolding.normalized:
                code += "    A *= fact\n"
                code += "    b *= fact\n"
                code += "    return A, b\n"
            else:
                code += "    A *= fact\n"
                code += "    return A\n"
        else:
            code += "    rhs *= fact\n"
            code += "    return rhs\n"
    # We write the code to file if provided, and execute it.
    if True:
        if file_name is not None:
            f = file(file_name+".py", "w")
            f.write(code)
            f.close()

        doppler_terms = code
        if not return_code:
            exec doppler_terms
    return doppler_terms


def fast_maxwell_boltzmann(mass, file_name=None,
                           return_code=False):
    r"""Return a function that returns values of a Maxwell-Boltzmann
    distribution.

    >>> from fast import Atom
    >>> mass = Atom("Rb", 87).mass
    >>> f = fast_maxwell_boltzmann(mass)
    >>> print f(0, 273.15+20)
    0.00238221482739

    >>> import numpy as np
    >>> v = np.linspace(-600, 600, 101)
    >>> dist = f(v, 273.15+20)
    >>> dv = v[1]-v[0]
    >>> print sum(dist)*dv
    0.999704711134

    """
    # We get the mass of the atom.
    code = ""
    code = "def maxwell_boltzmann(v, T):\n"
    code += '    r"""A fast calculation of the'
    code += ' Maxwell-Boltzmann distribution."""\n'
    code += "    if hasattr(v, 'shape'):\n"
    code += "        d = 1\n"
    code += "        m = %s\n" % mass
    code += "        f = np.sqrt(m/2/np.pi/k_B_num/T)**d\n"
    code += "        f = f * np.exp(-m*v**2/2/k_B_num/T)\n"
    code += "        return f\n"
    code += "    elif hasattr(v, '__len__'):\n"
    code += "        d = len(v)\n"
    code += "        m = %s\n" % mass
    code += "        f = np.sqrt(m/2/np.pi/k_B_num/T)**d\n"
    code += "        vsquare = sum([v[i]**2 for i in range(d)])\n"
    code += "        f = f * np.exp(-m*vsquare/2/k_B_num/T)\n"
    code += "        return f\n"
    code += "    else:\n"
    code += "        d = 1\n"
    code += "        m = %s\n" % mass
    code += "        f = np.sqrt(m/2/np.pi/k_B_num/T)**d\n"
    code += "        f = f * np.exp(-m*v**2/2/k_B_num/T)\n"
    code += "        return f\n"

    # We write the code to file if provided, and execute it.
    if file_name is not None:
        f = file(file_name+".py", "w")
        f.write(code)
        f.close()

    maxwell_boltzmann = code
    if not return_code:
        exec maxwell_boltzmann

    return maxwell_boltzmann


def fast_inhomo_bloch_equations(Ep, epsilonp, detuning_knob, T, gamma,
                                omega_level, rm, xi, theta,
                                unfolding, inhomogeneity, matrix_form=False,
                                file_name=None, return_code=False):
    r"""Return a fast function that returns the numeric right-hand sides of \
    inhomogeneous Bloch equations.

    We test a basic two-level system.

    >>> import numpy as np
    >>> from scipy.constants import physical_constants
    >>> from sympy import Matrix, symbols
    >>> from fast.electric_field import electric_field_amplitude_top
    >>> from fast.electric_field import PlaneWave
    >>> from fast.symbolic import (define_laser_variables,
    ...                            polarization_vector)
    >>> from fast.atomic_structure import Atom
    >>> from fast.bloch import phase_transformation

    >>> Ne = 2
    >>> Nl = 1
    >>> a0 = physical_constants["Bohr radius"][0]
    >>> rm = [np.array([[0, 0], [a0, 0]]),
    ...       np.array([[0, 0], [0, 0]]),
    ...       np.array([[0, 0], [0, 0]])]
    >>> xi = np.array([[[0, 1], [1, 0]]])
    >>> omega_level = [0, 2.4e15]
    >>> gamma21 = 2*np.pi*6e6
    >>> gamma = np.array([[0, -gamma21], [gamma21, 0]])
    >>> theta = phase_transformation(Ne, Nl, rm, xi)

    We define symbolic variables to be used as token arguments.
    >>> Ep, omega_laser = define_laser_variables(Nl)
    >>> laser = PlaneWave(0, 0, 0, 0)
    >>> epsilonp = [laser.epsilonp]
    >>> k = [laser.k]
    >>> detuning_knob = [symbols("delta1", real=True)]

    A map to unfold the density matrix.
    >>> unfolding = Unfolding(Ne, True, True, True)

    We define the Doppler broadening.

    >>> shape = [9]
    >>> stds = [[-4, 4]]
    >>> T = 273.15+20
    >>> mass = Atom("Rb", 87).mass
    >>> aux = (shape, stds, T, mass, detuning_knob, k,
    ...        omega_level, xi, theta, unfolding, ["z", "x", "y"],
    ...        True)

    >>> doppler_effect = DopplerBroadening(*aux)
    >>> doppler_effect.domain
    [array([-669.86784872, -502.40088654, -334.93392436, -167.46696218,
              0.        ,  167.46696218,  334.93392436,  502.40088654,
            669.86784872])]

    We obtain a function to calculate the Bloch equations.

    >>> T_symb = symbols("T", positive=True)
    >>> aux = (Ep, epsilonp, detuning_knob, T_symb, gamma,
    ...        omega_level, rm, xi, theta, unfolding, doppler_effect,
    ...        True)
    >>> bloch_equations = fast_inhomo_bloch_equations(*aux)

    We calculate an example.

    >>> detuning_knobs = [0]
    >>> Eps = electric_field_amplitude_top(0, 1e-3, 1, "SI")
    >>> Eps *= np.exp(1j*np.pi)
    >>> Eps = [Eps]
    >>> A, b = bloch_equations(Eps, detuning_knobs, T)
    >>> print A[:, 2, 1]*1e-6/2/np.pi
    [ 853.49268666  640.12094531  426.74705849  213.37341005    0.
     -213.37317167 -426.74610495 -640.11879984 -853.49125636]
    >>> print b*1e-6
    [[ 0.  0.  0.]
     [ 0.  0.  0.]
     [ 0.  0.  0.]
     [ 0.  0.  0.]
     [ 0.  0.  0.]
     [ 0.  0.  0.]
     [ 0.  0.  0.]
     [ 0.  0.  0.]
     [ 0.  0.  0.]]

    """
    if not unfolding.lower_triangular:
        mes = "It is very inefficient to solve using all components of the "
        mes += "density matrix. Better set lower_triangular=True in Unfolding."
        raise NotImplementedError(mes)
    if matrix_form and (not unfolding.real) and (unfolding.lower_triangular):
        mes = "It is not possible to express the equations in matrix form "
        mes += "for complex lower triangular components only."
        raise ValueError(mes)
    Nl = len(Ep)
    Nrho = unfolding.Nrho
    # We determine which arguments are constants.
    if True:
        try:
            Ep = np.array([complex(Ep[l]) for l in range(Nl)])
            variable_Ep = False
        except:
            variable_Ep = True

        try:
            epsilonp = [np.array([complex(epsilonp[l][i]) for i in range(3)])
                        for l in range(Nl)]
            variable_epsilonp = False
        except:
            variable_epsilonp = True
        try:
            detuning_knob = np.array([float(detuning_knob[l])
                                      for l in range(Nl)])
            variable_detuning_knob = False
        except:
            variable_detuning_knob = True
        try:
            T = float(T)
            variable_T = False
        except:
            variable_T = True
    # We obtain code for the homogeneous terms.
    if True:
        if file_name is not None:
            file_name_bloch = file_name+"_bloch"
        else:
            file_name_bloch = file_name
        aux = (Ep, epsilonp, detuning_knob, gamma, omega_level, rm, xi, theta,
               unfolding, matrix_form, file_name_bloch, True)
        bloch_equations = fast_bloch_equations(*aux)
        code = bloch_equations+"\n\n"
    # We establish the arguments of the output function.
    if True:
        code += "def inhomo_bloch_equations("
        code_args = ""
        if not matrix_form: code_args += "rho, "
        if variable_Ep: code_args += "Ep, "
        if variable_epsilonp: code_args += "epsilonp, "
        if variable_detuning_knob: code_args += "detuning_knob, "
        code += code_args
        if variable_T: code += "T, "
        code += "inhomogeneity=inhomogeneity, "
        code += "bloch_equations=bloch_equations):\n"
        code += '    r"""A fast calculation of inhomogeneous '
        code += 'Bloch equations."""\n'
    # We initialize the output and auxiliaries.
    if True:
        # We introduce the factor that multiplies all terms.
        sha = str(inhomogeneity.shape)[1:-1]+" "
        if matrix_form:
            code += "    A = np.zeros(("+sha+str(Nrho)+", "+str(Nrho)
            if not unfolding.real:
                code += "), complex)\n\n"
            else:
                code += "))\n\n"
            if unfolding.normalized:
                code += "    b = np.zeros(("+sha+str(Nrho)
                if not unfolding.real:
                    code += "), complex)\n\n"
                else:
                    code += "))\n\n"
        else:
            code += "    rhs = np.zeros(("+sha+str(Nrho)
            if not unfolding.real:
                code += "), complex)\n\n"
            else:
                code += "))\n\n"
    # We calculate the equations for each ensemble.
    if True:
        if variable_T: code += "    inhomogeneity.reset(T)\n"
        if code_args[-2:] == ", ": code_args = code_args[:-2]
        code += "    homogeneous = bloch_equations("+code_args+")\n\n"
        code += "    terms = inhomogeneity.terms\n"
        code += "    shape = inhomogeneity.shape\n"
        code += "    domain = inhomogeneity.domain\n"

        shape = inhomogeneity.shape
        dimension = len(shape)
        if dimension == 1:
            code += "    for i in range(shape[0]):\n"
            code += "        result = terms(domain[0][i], detuning_knob)\n"
            if matrix_form:
                if unfolding.normalized:
                    code += "        A[i] = homogeneous[0]+result[0]\n"
                    code += "        b[i] = homogeneous[1]+result[1]\n"
                else:
                    code += "        A[i] = homogeneous+result\n"
            else:
                code += "        rhs[i] = homogeneous+result\n"
    # We finish the code.
    if True:
        # code = rabi_code + "\n\n" + code
        if matrix_form:
            if unfolding.normalized:
                code += "    return A, b\n"
            else:
                code += "    return A\n"
        else:
            code += "    return rhs\n"
    # We write the code to file if provided, and execute it.
    if True:
        if file_name is not None:
            f = file(file_name+".py", "w")
            f.write(code)
            f.close()

        inhomo_bloch_equations = code
        if not return_code:
            exec inhomo_bloch_equations
    return inhomo_bloch_equations


def fast_inhomo_time_evolution(Ep, epsilonp, detuning_knob, gamma,
                               omega_level, rm, xi, theta, inhomogeneity,
                               semi_analytic=True,
                               file_name=None, return_code=False):
    r"""Return a fast function to calculate the time evolution of a state.

    We test a basic two-level system.

    >>> from fast.bloch import phase_transformation
    >>> from fast import PlaneWave, electric_field_amplitude_top, Atom
    >>> Ne = 2
    >>> Nl = 1
    >>> a0 = physical_constants["Bohr radius"][0]
    >>> rm = [np.array([[0, 0], [a0, 0]]),
    ...       np.array([[0, 0], [0, 0]]),
    ...       np.array([[0, 0], [0, 0]])]
    >>> xi = np.array([[[0, 1], [1, 0]]])
    >>> omega_level = [0, 2.4e15]
    >>> gamma21 = 2*np.pi*6e6
    >>> gamma = np.array([[0, -gamma21], [gamma21, 0]])
    >>> theta = phase_transformation(Ne, Nl, rm, xi)

    >>> Ep, omega_laser = define_laser_variables(Nl)
    >>> laser = PlaneWave(0, 0, 0, 0)
    >>> epsilonp = [laser.epsilonp]
    >>> k = [laser.k]
    >>> detuning_knob = [symbols("delta1", real=True)]

    A map to unfold the density matrix.
    >>> unfolding = Unfolding(Ne, True, True, True)

    We obtain a function to calculate time evolution.
    >>> aux = (Ep, epsilonp, detuning_knob, gamma,
    ...        omega_level, rm, xi, theta)
    >>> time_evolution = fast_time_evolution(*aux)

    >>> detuning_knobs = [0]
    >>> Eps = electric_field_amplitude_top(1e-3, 1e-3, 1, "SI")
    >>> Eps *= np.exp(1j*np.pi)
    >>> Eps = [Eps]

    >>> t = np.linspace(0, 1e-6, 11)
    >>> rho0 = np.array([[1, 0], [0, 0]])
    >>> rho0 = unfolding(rho0)

    >>> rhov0 = time_evolution(t, rho0, Eps, detuning_knobs)
    >>> print rhov0
    [[  0.00000000e+00   0.00000000e+00   0.00000000e+00]
     [  2.65923407e-01  -2.10635455e-17  -3.75764598e-01]
     [  2.61242286e-01  -6.34832515e-18  -3.52144243e-01]
     [  2.60791600e-01  -5.15628720e-18  -3.53251746e-01]
     [  2.60854215e-01  -4.93108886e-18  -3.53220826e-01]
     [  2.60849370e-01  -4.89778155e-18  -3.53220121e-01]
     [  2.60849648e-01  -4.89278818e-18  -3.53220302e-01]
     [  2.60849636e-01  -4.89202035e-18  -3.53220286e-01]
     [  2.60849637e-01  -4.89190454e-18  -3.53220287e-01]
     [  2.60849637e-01  -4.89188691e-18  -3.53220286e-01]
     [  2.60849637e-01  -4.89188423e-18  -3.53220286e-01]]

    We define the Doppler broadening.
    >>> Nvz = 3
    >>> shape = [Nvz]
    >>> stds = [[-4, 4]]
    >>> T = 273.15+20
    >>> mass = Atom("Rb", 87).mass
    >>> aux = (shape, stds, T, mass, detuning_knob, k,
    ...    omega_level, xi, theta, unfolding, ["z", "x", "y"],
    ...    True)

    >>> doppler_effect = DopplerBroadening(*aux)

    >>> aux = (Ep, epsilonp, detuning_knob, gamma,
    ...        omega_level, rm, xi, theta, doppler_effect,
    ...        True)

    >>> time_evolution = fast_inhomo_time_evolution(*aux)
    >>> rhot = time_evolution(t, rho0, Eps, detuning_knobs)
    >>> rhov0_inhomo = np.swapaxes(rhot, 0, 1)[Nvz/2]
    >>> print rhov0_inhomo-rhov0
    [[ 0.  0.  0.]
     [ 0.  0.  0.]
     [ 0.  0.  0.]
     [ 0.  0.  0.]
     [ 0.  0.  0.]
     [ 0.  0.  0.]
     [ 0.  0.  0.]
     [ 0.  0.  0.]
     [ 0.  0.  0.]
     [ 0.  0.  0.]
     [ 0.  0.  0.]]

    >>> print np.swapaxes(rhot, 0, 1)
    [[[  0.00000000e+00   0.00000000e+00   0.00000000e+00]
      [  8.10014861e-06   2.82707002e-03  -3.28172407e-04]
      [  6.83586416e-06   2.61409347e-03   4.78186293e-05]
      [  6.69318379e-06   2.58707370e-03  -1.19691879e-05]
      [  6.74368244e-06   2.59682364e-03  -9.92679502e-06]
      [  6.73782860e-06   2.59569997e-03  -8.91437415e-06]
      [  6.73771085e-06   2.59567649e-03  -9.14281734e-06]
      [  6.73786684e-06   2.59570660e-03  -9.12523444e-06]
      [  6.73784161e-06   2.59570175e-03  -9.12311787e-06]
      [  6.73784253e-06   2.59570193e-03  -9.12390236e-06]
      [  6.73784295e-06   2.59570201e-03  -9.12381062e-06]]
    <BLANKLINE>
     [[  0.00000000e+00   0.00000000e+00   0.00000000e+00]
      [  2.65923407e-01  -2.10635455e-17  -3.75764598e-01]
      [  2.61242286e-01  -6.34832515e-18  -3.52144243e-01]
      [  2.60791600e-01  -5.15628720e-18  -3.53251746e-01]
      [  2.60854215e-01  -4.93108886e-18  -3.53220826e-01]
      [  2.60849370e-01  -4.89778155e-18  -3.53220121e-01]
      [  2.60849648e-01  -4.89278818e-18  -3.53220302e-01]
      [  2.60849636e-01  -4.89202035e-18  -3.53220286e-01]
      [  2.60849637e-01  -4.89190454e-18  -3.53220287e-01]
      [  2.60849637e-01  -4.89188691e-18  -3.53220286e-01]
      [  2.60849637e-01  -4.89188423e-18  -3.53220286e-01]]
    <BLANKLINE>
     [[  0.00000000e+00   0.00000000e+00   0.00000000e+00]
      [  8.10014861e-06  -2.82707002e-03  -3.28172407e-04]
      [  6.83586416e-06  -2.61409347e-03   4.78186293e-05]
      [  6.69318379e-06  -2.58707370e-03  -1.19691879e-05]
      [  6.74368244e-06  -2.59682364e-03  -9.92679502e-06]
      [  6.73782860e-06  -2.59569997e-03  -8.91437415e-06]
      [  6.73771085e-06  -2.59567649e-03  -9.14281734e-06]
      [  6.73786684e-06  -2.59570660e-03  -9.12523444e-06]
      [  6.73784161e-06  -2.59570175e-03  -9.12311787e-06]
      [  6.73784253e-06  -2.59570193e-03  -9.12390236e-06]
      [  6.73784295e-06  -2.59570201e-03  -9.12381062e-06]]]

    """
    if not semi_analytic:
        s = "The numeric integrator has not been implemented."
        raise NotImplementedError(s)
    # We unpack variables.
    if True:
        Ne = len(omega_level)
        Nl = xi.shape[0]
        unfolding = Unfolding(Ne, True, True, True)
        Nrho = unfolding.Nrho
    # We determine which arguments are constants.
    if True:
        try:
            Ep = np.array([complex(Ep[l]) for l in range(Nl)])
            variable_Ep = False
        except:
            variable_Ep = True

        try:
            epsilonp = [np.array([complex(epsilonp[l][i]) for i in range(3)])
                        for l in range(Nl)]
            variable_epsilonp = False
        except:
            variable_epsilonp = True
        try:
            detuning_knob = np.array([float(detuning_knob[l])
                                      for l in range(Nl)])
            variable_detuning_knob = False
        except:
            variable_detuning_knob = True
    # We obtain code for the three parts.
    if True:
        args = (Ep, epsilonp, detuning_knob, gamma,
                omega_level, rm, xi, theta,
                unfolding, True, None, True)

        bloch_equations = fast_bloch_equations(*args)

        code = bloch_equations+"\n\n"

        if ((not variable_Ep) and
           (not variable_epsilonp) and
           (not variable_detuning_knob)):
            # We can call bloch_equations here!
            code += "bloch_equations = bloch_equations()\n"
    # We establish the arguments of the output function.
    if True:
        code += "def inhomo_time_evolution(t, rho0, "
        if variable_Ep: code += "Ep, "
        if variable_epsilonp: code += "epsilonp, "
        if variable_detuning_knob: code += "detuning_knob, "
        code += "bloch_equations=bloch_equations, "
        code += "inhomogeneity=inhomogeneity):\n"
        code += '    r"""A fast calculation of '
        code += 'inhomogeneous time evolution."""\n'
    # We call the Bloch equations.
    if True:
        shape = inhomogeneity.shape
        sha = str(shape)[1:-1]+" "
        code += "    Nt = t.shape[0]\n"
        code += "    terms = inhomogeneity.terms\n"
        code += "    shape = inhomogeneity.shape\n"
        code += "    domain = inhomogeneity.domain\n"
        code += "    rho = np.zeros((Nt, %s%s))\n" % (sha, Nrho)
        ###################################################################
        # We put in the arguments of bloch_equations.
        code += "    A_homo, b_homo = bloch_equations"
        if ((not variable_Ep) and
           (not variable_epsilonp) and
           (not variable_detuning_knob)):

            code += "\n"
        else:
            code += "("
            if variable_Ep: code += "Ep, "
            if variable_epsilonp: code += "epsilonp, "
            if variable_detuning_knob: code += "detuning_knob, "
            if code[-2:] == ", ": code = code[:-2]
            code += ")\n"
        ###################################################################
        dimension = len(shape)
        if dimension == 1:
            indent = "    "*2
            code += "    for i in range(shape[0]):\n"
            code += "        A_inhomo, b_inhomo = terms("
            code += "domain[0][i]"
            code += ", detuning_knob)\n"
            code += "        A = A_homo+A_inhomo\n"
            code += "        b = b_homo+b_inhomo\n"
        elif dimension == 2:
            indent = "    "*3
            code += "    for i in range(shape[0]):\n"
            code += "        for j in range(shape[1]):\n"
            code += "            A_inhomo, b_inhomo = terms("
            code += "[domain[0][i, j], domain[1][i, j]]"
            code += ", detuning_knob)\n"
            code += "            A = A_homo+A_inhomo\n"
            code += "            b = b_homo+b_inhomo\n"
        elif dimension == 3:
            indent = "    "*4
            code += "    for i in range(shape[0]):\n"
            code += "        for j in range(shape[1]):\n"
            code += "            for k in range(shape[2]):\n"
            code += "                A_inhomo, b_inhomo = terms("
            code += "[domain[0][i, j, k], "
            code += "domain[1][i, j, k], "
            code += "domain[2][i, j, k]]"
            code += ", detuning_knob)\n"
            code += "                A = A_homo+A_inhomo\n"
            code += "                b = b_homo+b_inhomo\n"

        code += indent+"lam, S = np.linalg.eig(A)\n"

        # code += indent+"print\n"

        code += indent+"Sinv = np.linalg.inv(S)\n"
        code += indent+"d = np.dot(Sinv, b)/lam\n"
        code += indent+"r = np.dot(Sinv, rho0) - d\n"
        code += indent+"rho_steady = np.dot(S, d)\n"

        # code += indent+"if i == 1: print A, list(lam)\n"

        code += indent+"for ii, ti in enumerate(t):\n"
        code += indent+"    rho_prime = r*np.exp(lam*ti)\n"
        code += indent+"    rho[ii, "
        if dimension == 1: code += "i, "
        if dimension == 2: code += "i, j, "
        if dimension == 3: code += "i, j, k, "
        code += ":] = np.real(np.dot(S, rho_prime) "
        code += "+ rho_steady)\n"

        code += "    return rho\n"
    # We write the code to file if provided, and execute it.
    if True:
        if file_name is not None:
            f = file(file_name+".py", "w")
            f.write(code)
            f.close()

        inhomo_time_evolution = code
        if not return_code:
            exec inhomo_time_evolution
    return inhomo_time_evolution


def fast_inhomo_sweep_time_evolution(Ep, epsilonp, gamma,
                                     omega_level, rm, xi, theta,
                                     inhomogeneity,
                                     semi_analytic=True,
                                     file_name=None, return_code=False):
    r"""Return a spectrum of time evolutions of the density matrix.

    We test a basic two-level system.

    >>> from fast.bloch import phase_transformation
    >>> from fast import PlaneWave, electric_field_amplitude_top, Atom
    >>> Ne = 2
    >>> Nl = 1
    >>> a0 = physical_constants["Bohr radius"][0]
    >>> rm = [np.array([[0, 0], [a0, 0]]),
    ...       np.array([[0, 0], [0, 0]]),
    ...       np.array([[0, 0], [0, 0]])]
    >>> xi = np.array([[[0, 1], [1, 0]]])
    >>> omega_level = [0, 2.4e15]
    >>> gamma21 = 2*np.pi*6e6
    >>> gamma = np.array([[0, -gamma21], [gamma21, 0]])
    >>> theta = phase_transformation(Ne, Nl, rm, xi)

    >>> Ep, omega_laser = define_laser_variables(Nl)
    >>> laser = PlaneWave(0, 0, 0, 0)
    >>> epsilonp = [laser.epsilonp]
    >>> k = [laser.k]
    >>> detuning_knob = [symbols("delta1", real=True)]

    A map to unfold the density matrix.
    >>> unfolding = Unfolding(Ne, True, True, True)

    >>> Eps = electric_field_amplitude_top(1e-3, 1e-3, 1, "SI")
    >>> Eps = [Eps]

    >>> t = np.linspace(0, 1e-6, 11)
    >>> rho0 = np.array([[1, 0], [0, 0]])
    >>> rho0 = unfolding(rho0)

    We define the Doppler broadening.
    >>> Nvz = 15
    >>> shape = [Nvz]
    >>> stds = [[-4, 4]]
    >>> T = 273.15+20
    >>> mass = Atom("Rb", 87).mass
    >>> aux = (shape, stds, T, mass, detuning_knob, k,
    ...    omega_level, xi, theta, unfolding, ["z", "x", "y"],
    ...    True)

    >>> doppler_effect = DopplerBroadening(*aux)

    We get a function for the frequency sweep of time evolution.
    >>> aux = (Ep, epsilonp, gamma,
    ...        omega_level, rm, xi, theta,
    ...        doppler_effect,
    ...        True,
    ...        "eqs")

    >>> inhomo_time_evolution = fast_inhomo_sweep_time_evolution(*aux)

    >>> amp = 1000e6*2*np.pi
    >>> Ndelta = 101
    >>> detuning_knobs = [[-amp, amp, Ndelta]]

    >>> deltas, rhot = inhomo_time_evolution(t, rho0, Eps, detuning_knobs)
    >>> print rhot.shape
    (101, 11, 15, 3)

    """
    # We unpack variables.
    if True:
        Nl = xi.shape[0]
    # We determine which arguments are constants.
    if True:
        try:
            Ep = np.array([complex(Ep[l]) for l in range(Nl)])
            variable_Ep = False
        except:
            variable_Ep = True

        try:
            epsilonp = [np.array([complex(epsilonp[l][i]) for i in range(3)])
                        for l in range(Nl)]
            variable_epsilonp = False
        except:
            variable_epsilonp = True
    # We obtain code for the time evolution.
    if True:
        detuning_knob = symbols("delta1:"+str(Nl))
        args = (Ep, epsilonp, detuning_knob, gamma, omega_level, rm, xi, theta,
                file_name, True)

        args = (Ep, epsilonp, detuning_knob, gamma, omega_level, rm, xi,
                theta, inhomogeneity, True, file_name, True)
        inhomo_time_evolution = fast_inhomo_time_evolution(*args)
        code = inhomo_time_evolution+"\n\n"
    # We establish the arguments of the output function.
    if True:
        code += "def inhomo_sweep_time_evolution(t, rho0, "
        if variable_Ep: code += "Ep, "
        if variable_epsilonp: code += "epsilonp, "
        code += "detuning_knob, "
        code += "inhomo_time_evolution=inhomo_time_evolution):\n"
        code += '    r"""A fast frequency sweep of the steady state."""\n'
    # Code to determine the sweep range.
    if True:
        code += """    sweepN = -1\n"""
        code += """    for i, delta in enumerate(detuning_knob):\n"""
        code += """        if hasattr(delta, "__getitem__"):\n"""
        code += """            sweepN = i\n"""
        code += """            delta0 = delta[0]\n"""
        code += """            deltaf = delta[1]\n"""
        code += """            Ndelta = delta[2]\n"""
        code += """            break\n\n"""
        code += """    if sweepN == -1:\n"""
        code += """        s = 'One of the detuning knobs '\n"""
        code += """        s += 'must be of the form '\n"""
        code += """        s += '(start, stop, Nsteps)'\n"""
        code += """        raise ValueError(s)\n\n"""
        code += """    deltas = np.linspace(delta0, deltaf, Ndelta)\n\n"""
    # We call time_evolution.
    if True:
        code += "    args = [[t, rho0, "
        if variable_Ep: code += "Ep, "
        if variable_epsilonp: code += "epsilonp, "
        code += """list(detuning_knob[:sweepN]) +\n"""
        code += """            [deltas[i]] +\n"""
        code += """            list(detuning_knob[sweepN+1:])]\n"""
        code += """          for i in range(Ndelta)]\n\n"""
        code += "    rho = np.array([inhomo_time_evolution(*argsi)\n"
        code += "                   for argsi in args])\n\n"
    # We finish the code.
    if True:
        code += "    return deltas, rho\n"
    # We write the code to file if provided, and execute it.
    if True:
        if file_name is not None:
            f = file(file_name+".py", "w")
            f.write(code)
            f.close()

        inhomo_sweep_time_evolution = code
        if not return_code:
            exec inhomo_sweep_time_evolution
    return inhomo_sweep_time_evolution


if __name__ == "__main__":
    import doctest
    print doctest.testmod(verbose=False)
