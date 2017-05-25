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

r"""This file contains all the information about the atoms.

This inclures routines to calculate the necessary matrices.

>>> g=State("Cs", 133, 6, 0, 1/Integer(2))
>>> make_list_of_states([g], "hyperfine")
[133Cs 6S_1/2^3, 133Cs 6S_1/2^4]

"""

from sympy.core.numbers import Rational as Integer
from math import sqrt, pi
from sympy.physics.wigner import wigner_3j, wigner_6j

# Physical constants (SI units):
from scipy.constants import physical_constants
from scipy.constants import pi as Pi  # Every one's favourite constant
#                                                              units
c = physical_constants["speed of light in vacuum"][0]      # m/s
a0 = physical_constants["Bohr radius"][0]                  # m
hbar = physical_constants["Planck constant over 2 pi"][0]  # J s
mu0 = physical_constants["mag. constant"][0]               # N / A^2
epsilon0 = physical_constants["electric constant"][0]      # F / m
e = physical_constants["elementary charge"][0]             # C
me = physical_constants["electron mass"][0]                # kg
k_B = physical_constants["Boltzmann constant"][0]          # J / K
uma = physical_constants["unified atomic mass unit"][0]    # kg


TmeltRb = 39.30+273.15  # K. [6]
TboilRb = 688+273.15    # K. [6]
TmeltCs = 28.5+273.15   # K. [3]
TboilCs = 671+273.15    # K. [3]

S = 'S'
P = 'P'
D = 'D'
F = 'F'
G = 'G'
H = 'H'
I = 'I'


class Atom(object):
    r"""This class implements specific atoms and their properties.

    The atoms can be identified by element and optionally by isotope.

    >>> atom=Atom("Rb")
    >>> print atom
    Rb
    """

    def __init__(self, element, isotope=None):
        r"""Initialise specific atoms and their properties.

        The atoms can be identified by element and optionally by isotope.

        >>> atom=Atom("Rb")
        >>> print atom
        Rb

        We can get the atomic number
        >>> atom.Z
        37

        We can get the available isotopes
        >>> atom.isotopes
        [85, 87]

        The atomic radius (in meters)
        >>> print atom.radius
        2.35e-10

        The melting temperature
        >>> atom.Tboil
        961.15

        The boiling temperature
        >>> atom.Tboil
        961.15

        If we also specify an isotope, new information becomes available
        >>> atom=Atom("Rb",85)
        >>> print atom
        85Rb

        The number of neutrons
        >>> atom.neutrons
        48

        The abundance of the isotope
        >>> print atom.abundance
        0.7217

        The mass of the atom (in kilograms)
        >>> print atom.mass
        1.40999341816e-25

        The nuclear spin of the atom
        >>> atom.nuclear_spin
        5/2

        The frequency of the ionization limit (in Hz)
        >>> print atom.ionization_frequency
        1.01002474142e+13

        """
        m_Rb85 = 84.9117897379*uma  # Rb85  mass in kg [4]
        m_Rb87 = 86.9091805310*uma  # Rb87  mass in kg [4]
        m_Cs133 = 132.9054519610*uma  # Cs133 mass in kg [4]

        abundance_Rb85 = 0.7217  # [5]
        abundance_Rb87 = 0.2783  # [5]
        abundance_Cs133 = 1.0    # [5]

        # Atomic radii from [7]
        # This is the database of implemented atoms.
        # nuclear spins from [5]
        # ionization frequencies from [3,6,8]
        #          element, isotope, atomic number,   mass, abundance      ,
        #          Tmelt  , Tboil  , radius (m) , nuclear spin,
        #          ionization frequencies

        database = [["Rb",      85,         37,  m_Rb85, abundance_Rb85,
                    TmeltRb, TboilRb, 2.35e-10, 5/Integer(2),
                    c*33690.79890],
                    ["Rb",      87,         37,  m_Rb87, abundance_Rb87,
                    TmeltRb, TboilRb, 2.35e-10, 3/Integer(2),
                    c*33690.80480],
                    ["Cs",     133,         55, m_Cs133, abundance_Cs133,
                    TmeltCs, TboilCs, 2.60e-10, 7/Integer(2),
                    c*31406.46766]]

        # We scan the database
        valid_input = False
        isotopes = []
        for item in database:
            if element == item[0] and isotope == item[1]:
                valid_input = True
                self.element = item[0]
                self.isotope = item[1]
                self.Z = item[2]
                self.mass = item[3]
                self.abundance = item[4]
                self.Tmelt = item[5]
                self.Tboil = item[6]
                self.radius = item[7]
                self.nuclear_spin = item[8]
                self.ionization_frequency = item[9]

                self.neutrons = self.isotope-self.Z
                break

            # If an isotope is not provided we return an object with a reduced
            # number of properties.
            if isotope is None:
                if element == item[0]:
                    isotopes += [item[1]]
                    valid_input = True
                    self.element = item[0]
                    self.Z = item[2]
                    self.Tmelt = item[5]
                    self.Tboil = item[6]
                    self.radius = item[7]
                    self.abundance = 1.0
        if isotope is None:
            self.isotope = None
            self.isotopes = isotopes

        if not valid_input:
            s = "The isotope "+str(isotope)+str(element)
            s += " is not in the database."
            raise ValueError, s

    def __repr__(self):
        r"""The representation routine for atoms.

        >>> Atom("Rb",85).__repr__()
        '85Rb'

        """
        if self.isotope is not None:
            return str(self.isotope)+self.element
        else:
            return self.element

    def __str__(self):
        r"""The string routine for atoms.

        >>> Atom("Cs",133).__str__()
        '133Cs'

        """
        return self.__repr__()

    def _latex_(self):
        r"""The LaTeX routine for atoms.

        >>> Atom("Rb",87)._latex_()
        '^{87}\\mathrm{Rb}'

        """
        return '^{'+str(self.isotope)+'}\\mathrm{'+self.element+'}'

    def states(self,
               Nmax=50, omega_min=None, omega_max=None, return_missing=False):
        r"""Find all states of available in an atom.

        This function returns all available states up to the fine structure
        (ordered by energy) such that the principal quantum number is N<=Nmax.
        Nmax is 50 by default.

        >>> atom=Atom("Rb",85)
        >>> states=atom.states()
        >>> print states
        [85Rb 5S_1/2, 85Rb 5P_1/2, 85Rb 5P_3/2, 85Rb 4D_5/2, 85Rb 4D_3/2, 85Rb 6S_1/2, 85Rb 6P_1/2, 85Rb 6P_3/2, 85Rb 5D_3/2, 85Rb 5D_5/2, 85Rb 7S_1/2, 85Rb 7P_1/2, 85Rb 7P_3/2, 85Rb 6D_3/2, 85Rb 7D_3/2, 85Rb 14S_1/2, 85Rb 15S_1/2, 85Rb 16S_1/2, 85Rb 17S_1/2, 85Rb 18S_1/2, 85Rb 19S_1/2, 85Rb 20S_1/2, 85Rb 21S_1/2, 85Rb 22S_1/2, 85Rb 23S_1/2, 85Rb 24S_1/2, 85Rb 25S_1/2, 85Rb 26S_1/2, 85Rb 27S_1/2, 85Rb 28S_1/2, 85Rb 29S_1/2, 85Rb 30S_1/2, 85Rb 31S_1/2, 85Rb 32S_1/2, 85Rb 33S_1/2, 85Rb 34S_1/2, 85Rb 35S_1/2, 85Rb 36S_1/2, 85Rb 37S_1/2, 85Rb 38S_1/2, 85Rb 39S_1/2, 85Rb 40S_1/2, 85Rb 41S_1/2, 85Rb 42S_1/2, 85Rb 43S_1/2, 85Rb 44S_1/2, 85Rb 45S_1/2, 85Rb 46S_1/2, 85Rb 47S_1/2, 85Rb 48S_1/2, 85Rb 49S_1/2, 85Rb 50S_1/2]

        If an omega_max is provided any state with an energy higher than
        hbar*omega will not be returned. If an omega_min is provided any state
        with an energy lower than hbar*omega will not be returned.

        >>> atom.states(omega_min=1.00845e15*2*pi, omega_max=1.0086e+15*2*pi)
        [85Rb 49S_1/2, 85Rb 50S_1/2]

        If return_missing=True then the function will return a 2-tuple composed
        by a list of the states available, and a list of the valid states not
        available.

        >>> available,not_available=atom.states(Nmax=5,return_missing=True)
        >>> print available
        [85Rb 5S_1/2, 85Rb 5P_1/2, 85Rb 5P_3/2, 85Rb 4D_5/2, 85Rb 4D_3/2, 85Rb 5D_3/2, 85Rb 5D_5/2]
        >>> print not_available
        [('Rb', 85, 1, 0, 1/2), ('Rb', 85, 2, 0, 1/2), ('Rb', 85, 2, 1, 1/2), ('Rb', 85, 2, 1, 3/2), ('Rb', 85, 3, 0, 1/2), ('Rb', 85, 3, 1, 1/2), ('Rb', 85, 3, 1, 3/2), ('Rb', 85, 3, 2, 3/2), ('Rb', 85, 3, 2, 5/2), ('Rb', 85, 4, 0, 1/2), ('Rb', 85, 4, 1, 1/2), ('Rb', 85, 4, 1, 3/2), ('Rb', 85, 4, 3, 5/2), ('Rb', 85, 4, 3, 7/2), ('Rb', 85, 5, 3, 5/2), ('Rb', 85, 5, 3, 7/2), ('Rb', 85, 5, 4, 7/2), ('Rb', 85, 5, 4, 9/2)]


        """
        # We generate all possible quantum numbers for N<=Nmax.

        S = 1/Integer(2)  # The spin of the electron.
        available = []
        not_available = []
        for N in range(1, Nmax+1):
            for L in range(N):
                Jmin = abs(L-S)
                Jmax = L+S
                Jpos = [Jmin+i for i in range(Jmax-Jmin+1)]
                for J in Jpos:
                    try:
                        state = State(self.element, self.isotope, N, L, J)
                        available += [state]
                    except:
                        not_available += [(self.element, self.isotope,
                                           N, L, J)]

        if omega_min is not None:
            available = [s for s in available if s.omega >= omega_min]
        if omega_max is not None:
            available = [s for s in available if s.omega <= omega_max]

        # We sort the states by energy.
        available = [(s.omega, s) for s in available]
        available = sorted(available)
        available = [s[1] for s in available]

        if return_missing:
            return available, not_available
        else:
            return available

    def transitions(self, omega_min=None, omega_max=None):
        r"""Find all allowed transitions.

        This function finds all allowed transitions (by electric-dipole
        selection rules) in the atom.

        >>> atom=Atom("Rb",85)
        >>> transitions=atom.transitions()
        >>> print len(transitions)
        270

        Arguments omega_min and omega_max can be used to make filter out the
        results.

        >>> from scipy.constants import c
        >>> wavelength_min=770e-9
        >>> wavelength_max=790e-9
        >>> omega_min=2*pi*c/wavelength_max
        >>> omega_max=2*pi*c/wavelength_min

        >>> easy_transitions=atom.transitions(omega_min=omega_min,omega_max=omega_max)
        >>> for ti in easy_transitions:
        ...     print abs(ti.wavelength)*1e9, ti
        780.241476935 85Rb 5S_1/2 -----> 85Rb 5P_3/2
        776.157015322 85Rb 5P_3/2 -----> 85Rb 5D_3/2
        775.978619616 85Rb 5P_3/2 -----> 85Rb 5D_5/2

        """
        states = self.states()
        transitions = states
        transitions = []
        for i in range(len(states)):
            si = states[i]
            for j in range(i):
                sj = states[j]
                t = Transition(sj, si, verbose=0)
                if t.allowed:
                    transitions += [t]

        if omega_min is not None:
            transitions = [ti for ti in transitions
                           if abs(ti.omega) >= omega_min]
        if omega_max is not None:
            transitions = [ti for ti in transitions
                           if abs(ti.omega) <= omega_max]

        return transitions

    def find_decays(self, fine_state):
        r"""Find all possible decays from a given fine state.

        This function finds all the states to which a given fine structure
        state can decay (through electric dipole selection rules).

        >>> atom=Atom("Cs",133)
        >>> e=State("Cs",133,6,"P",3/Integer(2))
        >>> atom.find_decays(e)
        [133Cs 6P_3/2, 133Cs 6S_1/2]

        >>> s=State("Cs",133,6,"D",5/Integer(2))
        >>> atom.find_decays(s)
        [133Cs 6D_5/2, 133Cs 6P_3/2, 133Cs 7P_3/2, 133Cs 6S_1/2, 133Cs 5D_3/2, 133Cs 5D_5/2, 133Cs 7S_1/2, 133Cs 6P_1/2]

        """
        def decays_from(fine_state, transitions):
            states_below = []
            for t in transitions:
                if t.e2 == fine_state:
                    states_below += [t.e1]
            return states_below

        def iterate(states, transitions):
            new_and_old_states = states[:]
            for state in states:
                new_states = decays_from(state, transitions)
                for new_state in new_states:
                    if new_state not in new_and_old_states:
                        new_and_old_states += [new_state]

            if states == new_and_old_states:
                return states
            else:
                return iterate(new_and_old_states, transitions)

        transitions = self.transitions()
        states = iterate([fine_state], transitions)

        return states


class State(object):
    r"""This class implements specific eigenstates of the atomic hamiltonian.

    The states must be identified by element, isotope and quantum numbers
    in one of these forms:

    1 .- |n, l, j>
    >>> State("Rb",85,5,0,1/Integer(2))
    85Rb 5S_1/2

    2 .- |n, l, j, f>
    >>> State("Cs", 133, 6, 0, 1/Integer(2), 3)
    133Cs 6S_1/2^3

    3 .- |n, l, j, f, m>
    >>> State("Rb", 87, 5, 0, 1/Integer(2), 2, -1)
    87Rb 5S_1/2^2,-1

    States have many useful properties:
    >>> g2=State("Cs", 133, 6, 0, 1/Integer(2), 4, 4)
    >>> e =State("Cs", 133, 6, 1, 3/Integer(2), 5, 5)

    The defining properties of the state:
    >>> g2.element, g2.isotope, g2.n, g2.l, g2.j, g2.f, g2.m
    ('Cs', 133, 6, 0, 1/2, 4, 4)

    The atomic number of the element:
    >>> g2.Z
    55

    The number of neutrons:
    >>> g2.neutrons
    78

    The mass of the atom (in kilograms):
    >>> print g2.mass
    2.2069469161e-25

    The absolute frequency of the state relative to the groundstate energy
    (in Hz):
    >>> print g2.nu
    4021776399.38

    The angular frequency of the state (in rad/s):
    >>> print g2.omega
    25269566381.3

    The hyperfine constants Ahfs, Bhfs, Chfs:
    >>> print e.Ahfs, e.Bhfs, e.Chfs
    50288250.0 -494000.0 560.0

    A latex representation:
    >>> print g2._latex_()
    ^{133}\mathrm{Cs}\ 6S_{1/2}^{4,4}
    """

    def __init__(self, element, isotope, n, l=None, j=None, f=None, m=None):
        r"""Initialize states.

        >>> State("Rb",85,5,0,1/Integer(2))
        85Rb 5S_1/2
        """
        # We declare the atom to which this state belongs.
        atom = Atom(element, isotope)
        self.atom = atom
        self.element = element
        self.isotope = isotope

        # We draw a few properties from the atom.
        self.Z = atom.Z
        self.neutrons = atom.neutrons
        self.abundance = atom.abundance
        self.mass = atom.mass

        if l is None:
            s = "The orbital angular momentum quantum number l "
            s += "must be specified."
            raise NotImplementedError, s
        if j is None:
            s = "The total angular momentum quantum number j=l+s "
            s += "must be specified."
            raise NotImplementedError, s

        # We go from chemical notation to quantum numbers.
        if str(l) in ['S', 's']: l = 0
        if str(l) in ['P', 'p']: l = 1
        if str(l) in ['D', 'd']: l = 2
        if str(l) in ['F', 'f']: l = 3
        if str(l) in ['G', 'g']: l = 4
        if str(l) in ['H', 'h']: l = 5
        if str(l) in ['I', 'i']: l = 6

        # We check the value of l.
        lperm = range(0, n)
        if l not in lperm:
            raise ValueError, 'l = ' + str(l) + ' is not allowed.'

        # We check the value of j.
        jmin = abs(l - Integer(1) / Integer(2))
        nj = int(2 * min(l, Integer(1) / Integer(2)) + 1)
        jperm = [jmin] + [jmin + ii for ii in range(1, nj)]
        if j not in jperm:
            raise ValueError, 'j = ' + str(j) + ' is not allowed.'

        self.n = n
        self.l = l
        self.j = j
        self.f = f
        self.m = m
        self.quantum_numbers = [isotope, n, l, j]

        # We calculate the energy of the state.
        ################################################################
        # the l values of letters in
        S = 0
        P = 1
        D = 2
        F = 3
        G = 4
        H = 5
        I = 6
        from scipy.constants import c
        # All tables are given in (cm^-1).
        i = atom.nuclear_spin
        if element == "Rb":
            if isotope==85:
                #        N, L,     K       , E (cm^-1),      A (cm^-1)      B (cm^-1)    C (cm^-1)
                nivfin=[[ 5, S, 1/Integer(2), 0.0000000  , 0.033753721     , 0.0       , 0.0],
                        [ 5, P, 1/Integer(2), 12578.950  , 0.004026        , 0.0       , 0.0],
                        [ 5, P, 3/Integer(2), 12816.545  , 0.0008352       , 0.0008676 , 0.0],

                        [ 4, D, 5/Integer(2), 19355.203  , 0.000168        , 0.0       , 0.0],
                        [ 4, D, 3/Integer(2), 19355.649  , 0.000244        , 0.0       , 0.0],

                        [ 6, S, 1/Integer(2), 20132.460  , 239.18e6/(c*100), 0.0       , 0.0],

                        [ 6, P, 1/Integer(2), 23715.081  , 0.001305        , 0.0       , 0.0],
                        [ 6, P, 3/Integer(2), 23792.591  , 0.0002723       , 0.0002732 , 0.0],

                        [ 5, D, 3/Integer(2), 25700.536  , 0.000140834     , 0.00006373, 0.0],
                        [ 5, D, 5/Integer(2), 25703.498  ,-0.00007309      , 0.0000894 , 0.0],

                        [ 7, S, 1/Integer(2), 26311.437  , 0.003159        , 0.0       , 0.0],

                        [ 7, P, 1/Integer(2), 27835.02   , 0.000590        , 0.0       , 0.0],
                        [ 7, P, 3/Integer(2), 27870.11   , 0.0001238       , 0.000123  , 0.0],

                        [ 6, D, 3/Integer(2), 28687.127  , 0.000077        , 0.000054  , 0.0],
                        [ 7, D, 3/Integer(2), 30280.113  , 0.0000472       , 0.000010  , 0.0],
                        [14, S, 1/Integer(2), 32761.60069, 0.0             , 0.0       , 0.0],

                        [15, S, 1/Integer(2), 32911.63322, 0.0             , 0.0       , 0.0],
                        [16, S, 1/Integer(2), 33028.05246, 0.0             , 0.0       , 0.0],
                        [17, S, 1/Integer(2), 33120.19975, 0.0             , 0.0       , 0.0],
                        [18, S, 1/Integer(2), 33194.38247, 0.0             , 0.0       , 0.0],
                        [19, S, 1/Integer(2), 33254.98464, 0.0             , 0.0       , 0.0],
                        [20, S, 1/Integer(2), 33305.12986, 0.0             , 0.0       , 0.0],
                        [21, S, 1/Integer(2), 33347.09244, 0.0             , 0.0       , 0.0],
                        [22, S, 1/Integer(2), 33382.56103, 0.0             , 0.0       , 0.0],
                        [23, S, 1/Integer(2), 33412.80991, 0.0             , 0.0       , 0.0],
                        [24, S, 1/Integer(2), 33438.81426, 0.0             , 0.0       , 0.0],
                        [25, S, 1/Integer(2), 33461.33387, 0.0             , 0.0       , 0.0],

                        [26, S, 1/Integer(2), 33480.96379, 0.0             , 0.0       , 0.0],
                        [27, S, 1/Integer(2), 33498.17857, 0.0             , 0.0       , 0.0],
                        [28, S, 1/Integer(2), 33513.35849, 0.0             , 0.0       , 0.0],
                        [29, S, 1/Integer(2), 33526.81233, 0.0             , 0.0       , 0.0],
                        [30, S, 1/Integer(2), 33538.79178, 0.0             , 0.0       , 0.0],
                        [31, S, 1/Integer(2), 33549.50518, 0.0             , 0.0       , 0.0],
                        [32, S, 1/Integer(2), 33559.12441, 0.0             , 0.0       , 0.0],
                        [33, S, 1/Integer(2), 33567.79401, 0.0             , 0.0       , 0.0],
                        [34, S, 1/Integer(2), 33575.63464, 0.0             , 0.0       , 0.0],
                        [35, S, 1/Integer(2), 33582.74865, 0.0             , 0.0       , 0.0],

                        [36, S, 1/Integer(2), 33589.22335, 0.0             , 0.0       , 0.0],
                        [37, S, 1/Integer(2), 33595.13297, 0.0             , 0.0       , 0.0],
                        [38, S, 1/Integer(2), 33600.54153, 0.0             , 0.0       , 0.0],
                        [39, S, 1/Integer(2), 33605.50403, 0.0             , 0.0       , 0.0],
                        [40, S, 1/Integer(2), 33610.06832, 0.0             , 0.0       , 0.0],
                        [41, S, 1/Integer(2), 33614.27566, 0.0             , 0.0       , 0.0],
                        [42, S, 1/Integer(2), 33618.16244, 0.0             , 0.0       , 0.0],
                        [43, S, 1/Integer(2), 33621.76072, 0.0             , 0.0       , 0.0],
                        [44, S, 1/Integer(2), 33625.09801, 0.0             , 0.0       , 0.0],
                        [45, S, 1/Integer(2), 33628.19886, 0.0             , 0.0       , 0.0],

                        [46, S, 1/Integer(2), 33631.08524, 0.0             , 0.0       , 0.0],
                        [47, S, 1/Integer(2), 33633.77664, 0.0             , 0.0       , 0.0],
                        [48, S, 1/Integer(2), 33636.29016, 0.0             , 0.0       , 0.0],
                        [49, S, 1/Integer(2), 33638.64086, 0.0             , 0.0       , 0.0],
                        [50, S, 1/Integer(2), 33640.84283, 0.0             , 0.0       , 0.0],
                        ]

            elif isotope == 87:
                #        N, L,     K       , E (cm^-1),      A (cm^-1)      B (cm^-1)      C (cm^-1)
                nivfin=[[5, S, 1/Integer(2), 0.0000000, 0.113990236053642, 0.0            , 0.0],
                        [5, P, 1/Integer(2), 12578.950, 0.01365          , 0.0            , 0.0],
                        [5, P, 3/Integer(2), 12816.545, 0.0028259        , 0.00041684     , 0.0],

                        [4, D, 5/Integer(2), 19355.203,-16.801e6/(c*100) , 3.645e6/(c*100), 0.0],# Not Sansonetti.
                        [4, D, 3/Integer(2), 19355.649, 24.750e6/(c*100) , 2.190e6/(c*100), 0.0],# Not Sansonetti.

                        [6, S, 1/Integer(2), 20132.460, 807.66e6/(c*100) , 0.0            , 0.0],# Not Sansonetti.

                        [6, P, 1/Integer(2), 23715.081, 0.0044217        , 0.0            , 0.0],
                        [6, P, 3/Integer(2), 23792.591, 0.000924         , 0.000132       , 0.0],

                        [5, D, 3/Integer(2), 25700.536, 0.00048134       , 0.00003109     , 0.0],
                        [5, D, 5/Integer(2), 25703.498, -0.00024886      , 0.00004241     , 0.0],

                        [7, S, 1/Integer(2), 26311.437, 0.010664         , 0.0            , 0.0],

                        [7, P, 1/Integer(2), 27835.020, 0.001999         , 0.0            , 0.0],
                        [7, P, 3/Integer(2), 27870.110, 0.0004193        , 0.00005700     , 0.0]
                        ]

            else:
                s = "The isotope "+str(isotope)+str(element)+" is not in the database."
                raise ValueError, s
            # We rewrite the table in Hz
            nivfin=[ nivfin[ii][:3] + [nivfin[ii][3]*c*100 ] +[nivfin[ii][4]*c*100 ] +[nivfin[ii][5]*c*100 ]
                                    + [nivfin[ii][6]*c*100 ] for ii in range(len(nivfin)) ]

        elif element == "Cs":
            if isotope == 133:
                # Reference [1], others not used yet [2]:
                #        N, L,     K       , E (cm^-1),       A (MHz)      B (MHz)   C (MHz)

                nivfin=[[ 6, S, 1/Integer(2),     0         , 2298.1579425, 0     , 1.0     ],# This is exact.

                        [ 6, P, 1/Integer(2), 11178.26815870,  291.9309   , 0.0   , 0.0     ],
                        [ 6, P, 3/Integer(2), 11732.3071041 ,  50.28825   ,-0.4940, 0.000560],# C: 0.000560 Steck

                        [ 5, D, 3/Integer(2), 14499.2568    ,  48.78      , 0.1   , 0.0     ],
                        [ 5, D, 5/Integer(2), 14596.84232   , -21.24      , 0.2   , 0.0     ],

                        [ 7, S, 1/Integer(2), 18535.5286    , 545.90      , 0.0   , 0.0     ],

                        [ 7, P, 1/Integer(2), 21765.348     ,  94.35      , 0.0   , 0.0     ],
                        [ 7, P, 3/Integer(2), 21946.397     ,  16.609     , 0.0   , 0.0     ],

                        [ 6, D, 3/Integer(2), 22588.8210    ,  16.34      ,-0.1   , 0.0     ],
                        [ 6, D, 5/Integer(2), 22631.6863    ,  -4.66      , 0.9   , 0.0     ],

                        [ 8, S, 1/Integer(2), 24317.149400  , 219.12      , 0.0   , 0.0     ],
                        [ 4, F, 7/Integer(2), 24472.0455    ,   0.0       , 0.0   , 0.0     ],
                        [ 4, F, 5/Integer(2), 24472.2269    ,   0.0       , 0.0   , 0.0     ],
                        [ 8, P, 1/Integer(2), 25708.85473   ,  42.97      , 0.0   , 0.0     ],
                        [ 8, P, 3/Integer(2), 25791.508     ,   7.626     , 0.0   , 0.0     ],
                        [ 7, D, 3/Integer(2), 26047.8342    ,   7.4       , 0.0   , 0.0     ],
                        [ 7, D, 5/Integer(2), 26068.7730    ,  -1.7       , 0.0   , 0.0     ],

                        [ 9, S, 1/Integer(2), 26910.6627    , 110.1       , 0.0   , 0.0     ],
                        [ 5, F, 7/Integer(2), 26971.1535    ,   0.0       , 0.0   , 0.0     ],
                        [ 5, F, 5/Integer(2), 26971.3030    ,   0.0       , 0.0   , 0.0     ],
                        [ 5, G, 7/Integer(2), 27008.0541    ,   0.0       , 0.0   , 0.0     ],
                        [ 5, G, 9/Integer(2), 27008.0569    ,   0.0       , 0.0   , 0.0     ],
                        [ 9, P, 1/Integer(2), 27636.9966    ,   0.0       , 0.0   , 0.0     ],
                        [ 9, P, 3/Integer(2), 27681.6782    ,  23.19      , 0.0   , 0.0     ],
                        [ 8, D, 3/Integer(2), 27811.2400    ,   4.129     , 0.0   , 0.0     ],
                        [ 8, D, 5/Integer(2), 27822.8802    ,   3.95      , 0.0   , 0.0     ],
                        [10, S, 1/Integer(2), 28300.2287    ,  -0.85      , 0.0   , 0.0     ],
                        [ 6, F, 7/Integer(2), 28329.4075    ,  63.2       , 0.0   , 0.0     ],
                        [ 6, F, 5/Integer(2), 28329.5133    ,   0.0       , 0.0   , 0.0     ],
                        [ 6, G, 7/Integer(2), 28352.4444    ,   0.0       , 0.0   , 0.0     ],
                        [ 6, G, 9/Integer(2), 28352.4460    ,   0.0       , 0.0   , 0.0     ],
                        [10, P, 1/Integer(2), 28726.8123    ,  13.9       , 0.0   , 0.0     ],
                        [10, P, 3/Integer(2), 28753.6769    ,   2.485     , 0.0   , 0.0     ],

                        [ 9, D, 3/Integer(2), 28828.6820    ,   2.38      , 0.0   , 0.0     ],
                        [ 9, D, 5/Integer(2), 28835.79192   ,  -0.45      , 0.0   , 0.0     ],
                        [11, S, 1/Integer(2), 29131.73004   ,  39.4       , 0.0   , 0.0     ],
                        [ 7, F, 7/Integer(2), 29147.90818   ,   0.0       , 0.0   , 0.0     ],
                        [ 7, F, 5/Integer(2), 29147.98188   ,   0.0       , 0.0   , 0.0     ],
                        [ 7, G, 7/Integer(2), 29163.07206   ,   0.0       , 0.0   , 0.0     ],
                        [ 7, G, 9/Integer(2), 29163.0731    ,   0.0       , 0.0   , 0.0     ],
                        [11, P, 1/Integer(2), 29403.42310   ,   0.0       , 0.0   , 0.0     ],
                        [11, P, 3/Integer(2), 29420.824     ,   1.600     , 0.0   , 0.0     ],
                        [10, D, 3/Integer(2), 29468.2878    ,   1.54      , 0.0   , 0.0     ],
                        [10, D, 5/Integer(2), 29472.93995   ,  -0.35      , 0.0   , 0.0     ],
                        [12, S, 1/Integer(2), 29668.80336   ,  26.31      , 0.0   , 0.0     ],

                        [ 8, F, 7/Integer(2), 29678.68970   ,   0.0       , 0.0   , 0.0     ],
                        [ 8, F, 5/Integer(2), 29678.74280   ,   0.0       , 0.0   , 0.0     ],
                        [ 8, G, 7/Integer(2), 29689.13795   ,   0.0       , 0.0   , 0.0     ],
                        [ 8, G, 9/Integer(2), 29689.1388    ,   0.0       , 0.0   , 0.0     ],
                        [12, P, 1/Integer(2), 29852.43153   ,   0.0       , 0.0   , 0.0     ],
                        [12, P, 3/Integer(2), 29864.345     ,   1.10      , 0.0   , 0.0     ],
                        [11, D, 3/Integer(2), 29896.3399    ,   1.055     , 0.0   , 0.0     ],
                        [11, D, 5/Integer(2), 29899.54646   ,   0.24      , 0.0   , 0.0     ],
                        [13, S, 1/Integer(2), 30035.78836   ,  18.4       , 0.0   , 0.0     ],
                        [ 9, F, 7/Integer(2), 30042.27515   ,   0.0       , 0.0   , 0.0     ],
                        [ 9, F, 5/Integer(2), 30042.31405   ,   0.0       , 0.0   , 0.0     ],
                        [ 9, G, 7/Integer(2), 30049.75317   ,   0.0       , 0.0   , 0.0     ],
                        [ 9, G, 9/Integer(2), 30049.7545    ,   0.0       , 0.0   , 0.0     ],
                        [13, P, 1/Integer(2), 30165.66826   ,   0.0       , 0.0   , 0.0     ],
                        [13, P, 3/Integer(2), 30174.178     ,   0.77      , 0.0   , 0.0     ],

                        [12, D, 3/Integer(2), 30196.7963    ,   0.758     , 0.0   , 0.0     ],
                        [12, D, 5/Integer(2), 30199.09821   ,   0.19      , 0.0   , 0.0     ],
                        [14, S, 1/Integer(2), 30297.64510   ,   13.4      , 0.0   , 0.0     ],
                        [14, F, 7/Integer(2), 30302.13624   ,   0.0       , 0.0   , 0.0     ],
                        [10, F, 5/Integer(2), 30302.16537   ,   0.0       , 0.0   , 0.0     ],
                        [10, G, 7/Integer(2), 30307.66076   ,   0.0       , 0.0   , 0.0     ],
                        [10, G, 9/Integer(2), 30307.6617    ,   0.0       , 0.0   , 0.0     ],
                        [10, P, 1/Integer(2), 30392.8718    ,   0.0       , 0.0   , 0.0     ],
                        [14, P, 3/Integer(2), 30399.163     ,   0.0       , 0.0   , 0.0     ],

                        [13, D, 3/Integer(2), 30415.7533    ,   0.556     , 0.0   , 0.0     ],
                        [13, D, 5/Integer(2), 30417.46075   ,   0.14      , 0.0   , 0.0     ],
                        [15, S, 1/Integer(2), 30491.02346   ,  10.1       , 0.0   , 0.0     ],
                        [11, F, 7/Integer(2), 30494.26583   ,   0.0       , 0.0   , 0.0     ],
                        [11, F, 5/Integer(2), 30494.28809   ,   0.0       , 0.0   , 0.0     ],
                        [11, G, 9/Integer(2), 30498.4556    ,   0.0       , 0.0   , 0.0     ],
                        [11, G, 7/Integer(2), 30498.45695   ,   0.0       , 0.0   , 0.0     ],
                        [15, P, 1/Integer(2), 30562.90893   ,   0.0       , 0.0   , 0.0     ],
                        [15, P, 3/Integer(2), 30567.688     ,   0.0       , 0.0   , 0.0     ],
                        [14, D, 3/Integer(2), 30580.2267    ,   0.425     , 0.0   , 0.0     ],
                        [14, D, 5/Integer(2), 30581.52758   ,   0.0       , 0.0   , 0.0     ],
                        [16, S, 1/Integer(2), 30637.88276   ,   7.73      , 0.0   , 0.0     ],

                        [12, F, 7/Integer(2), 30640.30287   ,   0.0       , 0.0   , 0.0     ],
                        [12, F, 5/Integer(2), 30640.32028   ,   0.0       , 0.0   , 0.0     ],
                        [12, G, 7/Integer(2), 30643.55484   ,   0.0       , 0.0   , 0.0     ],
                        [16, P, 1/Integer(2), 30693.47416   ,   0.0       , 0.0   , 0.0     ],
                        [16, P, 3/Integer(2), 30697.191     ,   0.0       , 0.0   , 0.0     ],
                        [15, D, 3/Integer(2), 30706.9003    ,   0.325     , 0.0   , 0.0     ],
                        [15, D, 5/Integer(2), 30707.91378   ,   0.0       , 0.0   , 0.0     ],
                        [17, S, 1/Integer(2), 30752.03412   ,   6.06      , 0.0   , 0.0     ],

                        [13, F, 7/Integer(2), 30753.89018   ,   0.0       , 0.0   , 0.0     ],
                        [13, F, 5/Integer(2), 30753.90406   ,   0.0       , 0.0   , 0.0     ],
                        [13, G, 7/Integer(2), 30756.46241   ,   0.0       , 0.0   , 0.0     ],
                        [17, P, 1/Integer(2), 30795.90702   ,   0.0       , 0.0   , 0.0     ],
                        [17, P, 3/Integer(2), 30798.852     ,   0.0       , 0.0   , 0.0     ],
                        [16, D, 3/Integer(2), 30806.5283    ,   0.255     , 0.0   , 0.0     ],
                        [16, D, 5/Integer(2), 30807.33297   ,   0.0       , 0.0   , 0.0     ],
                        [18, S, 1/Integer(2), 30842.51775   ,   0.0       , 0.0   , 0.0     ],
                        [14, F, 7/Integer(2), 30843.97365   ,   0.0       , 0.0   , 0.0     ],
                        [14, F, 5/Integer(2), 30843.98488   ,   0.0       , 0.0   , 0.0     ],
                        [14, G, 7/Integer(2), 30846.04223   ,   0.0       , 0.0   , 0.0     ],

                        [18, P, 1/Integer(2), 30877.74761   ,   0.0       , 0.0   , 0.0     ],
                        [18, P, 3/Integer(2), 30880.1228    ,   0.0       , 0.0   , 0.0     ],
                        [17, D, 3/Integer(2), 30886.2959    ,   0.190     , 0.0   , 0.0     ],
                        [17, D, 5/Integer(2), 30886.94513   ,   0.0       , 0.0   , 0.0     ],
                        [19, S, 1/Integer(2), 30915.45262   ,   0.0       , 0.0   , 0.0     ],
                        [15, F, 7/Integer(2), 30916.61661   ,   0.0       , 0.0   , 0.0     ],
                        [15, F, 5/Integer(2), 30916.62583   ,   0.0       , 0.0   , 0.0     ],
                        [15, G, 7/Integer(2), 30918.30448   ,   0.0       , 0.0   , 0.0     ],
                        [19, P, 1/Integer(2), 30944.16859   ,   0.0       , 0.0   , 0.0     ],
                        [19, P, 3/Integer(2), 30946.113     ,   0.0       , 0.0   , 0.0     ],
                        [18, D, 5/Integer(2), 30951.1511    ,   0.160     , 0.0   , 0.0     ],
                        [18, D, 3/Integer(2), 30951.68259   ,   0.0       , 0.0   , 0.0     ],
                        [20, S, 1/Integer(2), 30975.10034   ,   0.0       , 0.0   , 0.0     ],

                        [16, F, 7/Integer(2), 30976.04620   ,   0.0       , 0.0   , 0.0     ],
                        [16, F, 5/Integer(2), 30976.05385   ,   0.0       , 0.0   , 0.0     ],
                        [16, G, 7/Integer(2), 30977.44090   ,   0.0       , 0.0   , 0.0     ],
                        [20, P, 1/Integer(2), 30998.79      ,   0.0       , 0.0   , 0.0     ],
                        [20, P, 3/Integer(2), 31000.40      ,   0.0       , 0.0   , 0.0     ],
                        [19, D, 3/Integer(2), 31004.5900    ,   0.0       , 0.0   , 0.0     ],
                        [19, D, 5/Integer(2), 31005.03231   ,   0.0       , 0.0   , 0.0     ],
                        [21, S, 1/Integer(2), 31024.50355   ,   0.0       , 0.0   , 0.0     ],

                        [17, F, 7/Integer(2), 31025.28272   ,   0.0       , 0.0   , 0.0     ],
                        [17, F, 5/Integer(2), 31025.28907   ,   0.0       , 0.0   , 0.0     ],
                        [17, G, 7/Integer(2), 31026.44832   ,   0.0       , 0.0   , 0.0     ],
                        [21, P, 1/Integer(2), 31044.31315   ,   0.0       , 0.0   , 0.0     ],
                        [21, P, 3/Integer(2), 31045.664     ,   0.0       , 0.0   , 0.0     ],
                        [20, D, 3/Integer(2), 31049.1456    ,   0.0       , 0.0   , 0.0     ],
                        [20, D, 5/Integer(2), 31049.51701   ,   0.0       , 0.0   , 0.0     ],
                        [22, S, 1/Integer(2), 31065.88056   ,   0.0       , 0.0   , 0.0     ],

                        [18, F, 7/Integer(2), 31066.53042   ,   0.0       , 0.0   , 0.0     ],
                        [18, F, 5/Integer(2), 31066.53582   ,   0.0       , 0.0   , 0.0     ],
                        [18, G, 7/Integer(2), 31067.51435   ,   0.0       , 0.0   , 0.0     ],
                        [22, P, 1/Integer(2), 31082.5979    ,   0.0       , 0.0   , 0.0     ],# not experimental (isoelectronic fitting).
                        [22, P, 3/Integer(2), 31083.77      ,   0.0       , 0.0   , 0.0     ],
                        [21, D, 3/Integer(2), 31086.6824    ,   0.0       , 0.0   , 0.0     ],
                        [21, D, 5/Integer(2), 31086.99717   ,   0.0       , 0.0   , 0.0     ],
                        [23, S, 1/Integer(2), 31100.88052   ,   2.4       , 0.0   , 0.0     ],
                        [19, F, 7/Integer(2), 31101.42859   ,   0.0       , 0.0   , 0.0     ],
                        [19, F, 5/Integer(2), 31101.43321   ,   0.0       , 0.0   , 0.0     ],
                        [19, G, 7/Integer(2), 31102.26518   ,   0.0       , 0.0   , 0.0     ],# not experimental (isoelectronic fitting).
                        [23, P, 1/Integer(2), 31115.11733   ,   0.56      , 0.0   , 0.0     ],
                        [23, P, 3/Integer(2), 31116.0904    ,   0.0       , 0.0   , 0.0     ],
                        [22, D, 3/Integer(2), 31118.60105   ,   0.0       , 0.0   , 0.0     ],
                        [22, D, 5/Integer(2), 31118.86983   ,   0.0       , 0.0   , 0.0     ],
                        [24, S, 1/Integer(2), 31130.74987   ,   0.0       , 0.0   , 0.0     ],

                        [20, F, 7/Integer(2), 31131.21615   ,   0.0       , 0.0   , 0.0     ],
                        [20, F, 5/Integer(2), 31131.22012   ,   0.0       , 0.0   , 0.0     ],
                        [20, G, 7/Integer(2), 31131.93570   ,   0.0       , 0.0   , 0.0     ],
                        [24, P, 1/Integer(2), 31142.9734    ,   0.0       , 0.0   , 0.0     ],# not experimental (isoelectronic fitting).
                        [24, P, 3/Integer(2), 31143.84      ,   0.0       , 0.0   , 0.0     ],
                        [23, D, 3/Integer(2), 31145.9690    ,   0.0       , 0.0   , 0.0     ],
                        [23, D, 5/Integer(2), 31146.20007   ,   0.0       , 0.0   , 0.0     ],
                        [25, S, 1/Integer(2), 31156.44439   ,   1.4       , 0.0   , 0.0     ],

                        [21, F, 7/Integer(2), 31156.8447    ,   0.0       , 0.0   , 0.0     ],# not experimental (isoelectronic fitting).
                        [21, F, 5/Integer(2), 31156.8482    ,   0.0       , 0.0   , 0.0     ],# not experimental (isoelectronic fitting).
                        [21, G, 7/Integer(2), 31157.46595   ,   0.0       , 0.0   , 0.0     ],# not experimental (isoelectronic fitting).
                        [25, P, 1/Integer(2), 31167.01727   ,   0.40      , 0.0   , 0.0     ],
                        [25, P, 3/Integer(2), 31167.74257   ,   0.0       , 0.0   , 0.0     ],
                        [24, D, 3/Integer(2), 31169.6144    ,   0.0       , 0.0   , 0.0     ],
                        [24, D, 5/Integer(2), 31169.81187   ,   0.0       , 0.0   , 0.0     ],
                        [22, F, 7/Integer(2), 31179.05392   ,   0.0       , 0.0   , 0.0     ],
                        [22, F, 5/Integer(2), 31179.05692   ,   0.0       , 0.0   , 0.0     ],
                        [22, G, 7/Integer(2), 31179.59567   ,   0.0       , 0.0   , 0.0     ],
                        [25, D, 3/Integer(2), 31190.17569   ,   0.0       , 0.0   , 0.0     ],
                        [25, D, 5/Integer(2), 31190.35063   ,   0.0       , 0.0   , 0.0     ],
                        [23, F, 7/Integer(2), 31198.4257    ,   0.0       , 0.0   , 0.0     ],# not experimental (isoelectronic fitting).
                        [23, F, 5/Integer(2), 31198.4283    ,   0.0       , 0.0   , 0.0     ],# not experimental (isoelectronic fitting).
                        [23, G, 7/Integer(2), 31198.89936   ,   0.0       , 0.0   , 0.0     ],# not experimental (isoelectronic fitting).
                        [24, F, 7/Integer(2), 31215.42383   ,   0.0       , 0.0   , 0.0     ],
                        [24, F, 5/Integer(2), 31215.42620   ,   0.0       , 0.0   , 0.0     ],
                        [24, G, 7/Integer(2), 31215.84176   ,   0.0       , 0.0   , 0.0     ],
                        [25, F, 7/Integer(2), 31230.4208    ,   0.0       , 0.0   , 0.0     ],# not experimental (isoelectronic fitting).
                        [25, F, 5/Integer(2), 31230.4229    ,   0.0       , 0.0   , 0.0     ],# not experimental (isoelectronic fitting).
                        [25, G, 7/Integer(2), 31230.79013   ,   0.0       , 0.0   , 0.0     ],# not experimental (isoelectronic fitting).
                        ]

            else:
                s = "The isotope "+str(isotope)+str(element)
                s += " is not in the database."
                raise ValueError, s
            # We rewrite the table in Hz
            nivfin=[ nivfin[ii][:3] + [nivfin[ii][3]*c*100 ] +[nivfin[ii][4]*1e6 ]
                                    + [nivfin[ii][5]*1e6   ] +[nivfin[ii][6]*1e6 ] for ii in range(len(nivfin)) ]

        else:
            s = "The element "+str(element)+" is not in the database."
            raise ValueError, s

        self.i = i

        # We check whether the quantum numbers given are in the database.
        # We find the energy of the state up to fine structure.
        in_database = False
        for ii in range(len(nivfin)):
            # print [n,l,j],nivfin[ii][:3]
            if [n, l, j] == nivfin[ii][:3]:
                nufin = nivfin[ii][3]
                A = nivfin[ii][4]; B = nivfin[ii][5]; C = nivfin[ii][6]
                self.Ahfs = A; self.Bhfs = B; self.Chfs = C
                in_database = True
                break
        if not in_database:
            s = "The values of n,l,k: "+str(n)+", "+str(l)+", "+str(j)
            s += " are not in the database for "+str(isotope)+str(element)+"."
            raise ValueError, s

        # We check the value of f.
        fmin = int(abs(j-i)); nf = int(2*min(j, i)+1)
        fperm = [fmin]+[fmin + ii for ii in range(1, nf)]
        if f is not None:
            if f in fperm:
                self.quantum_numbers.append(f)
                mperm = [ii for ii in range(-f, f+1)]
                if m is not None:
                    if m in mperm:
                        self.quantum_numbers.append(m)
                    else:
                        raise ValueError, 'm = '+str(m)+' is not allowed.'
                else:
                    self.mperm = mperm
            else:
                raise ValueError, 'f = '+str(f)+' is not allowed.'
        else:
            self.fperm = fperm

        # We establish the hyperfine correction to the energies.
        if f is not None:
            K = f*(f+1)-i*(i+1)-j*(j+1)
            deltanu_hfin = A*K/Integer(2)
            if j != 1/Integer(2):
                num = (3*K*(K+1)/Integer(2) - 2*i*(i+1)*j*(j+1))
                den = (4*i*(2*i-1)*j*(2*j-1))
                deltanu_hfin += B*num/den

                num = (5*K**2*(K/4+1) + K*(i*(i+1) + j*(j+1) + 3 - 3*i*(i+1)*j*(j+1)) - 5*i*(i+1)*j*(j+1))
                den = (i*(i-1)*(2*i-1)*j*(j-1)*(2*j-1))

                deltanu_hfin = deltanu_hfin + C*num/den

            deltanu_hfin = float(deltanu_hfin)
            self.nu = nufin + deltanu_hfin
            self.omega = 2*pi*self.nu

        else:
            hyperfine_structure = []
            for f in fperm:
                K = f*(f+1)-i*(i+1)-j*(j+1)
                deltanu_hfin = A*K/Integer(2)
                if j != Integer(1)/Integer(2):
                    deltanu_hfin = deltanu_hfin + B*(3*K*(K+1)/Integer(2) - 2*i*(i+1)*j*(j+1))/(4*i*(2*i-1)*j*(2*j-1))

                    num = (5*K**2*(K/4+1) + K*(i*(i+1) + j*(j+1) + 3 - 3*i*(i+1)*j*(j+1)) - 5*i*(i+1)*j*(j+1))
                    den = (i*(i-1)*(2*i-1)*j*(j-1)*(2*j-1))

                    deltanu_hfin = deltanu_hfin + C*num/den

                hyperfine_structure.append((deltanu_hfin, f))
            hyperfine_structure.sort()

            self.hyperfine_structure = hyperfine_structure
            self.nu = nufin
            self.omega = 2*pi*self.nu

    def __repr__(self):
        r"""The representation routine for states.

        >>> State("Rb",85,5,0,1/Integer(2)).__repr__()
        '85Rb 5S_1/2'

        >>> State("Rb",85,5,0,1/Integer(2),2).__repr__()
        '85Rb 5S_1/2^2'

        >>> State("Rb",85,5,0,1/Integer(2),2,2).__repr__()
        '85Rb 5S_1/2^2,2'

        """
        if self.l == 0: l = 'S'
        elif self.l == 1: l = 'P'
        elif self.l == 2: l = 'D'
        elif self.l == 3: l = 'F'
        elif self.l == 4: l = 'G'
        elif self.l == 5: l = 'H'
        elif self.l == 6: l = 'I'
        else: l = "(L="+str(self.l)+")"

        if self.f is None:
            s = str(self.isotope)+self.element+' '+str(self.n)+l+'_'
            s += str(self.j)
        else:
            s = str(self.isotope)+self.element+' '+str(self.n)+l+'_'
            s += str(self.j)+'^'+str(self.f)

        if self.m is not None:
            s += ','+str(self.m)
        return s

    def _latex_(self):
        r"""The LaTeX routine for states.

        >>> State("Rb",85,5,0,1/Integer(2))._latex_()
        '^{85}\\mathrm{Rb}\\ 5S_{1/2}'

        >>> State("Rb",85,5,0,1/Integer(2),2)._latex_()
        '^{85}\\mathrm{Rb}\\ 5S_{1/2}^{2}'

        >>> State("Rb",85,5,0,1/Integer(2),2,2)._latex_()
        '^{85}\\mathrm{Rb}\\ 5S_{1/2}^{2,2}'

        """
        if self.l == 0: l = 'S'
        elif self.l == 1: l = 'P'
        elif self.l == 2: l = 'D'
        elif self.l == 3: l = 'F'
        else: l = str(self.l)

        if self.f is None:
            s = '^{'+str(self.isotope)+'}\\mathrm{'+self.element+'}\\ '
            s += str(self.n)+l+'_{'+str(self.j)+'}'
        else:
            s = '^{'+str(self.isotope)+'}\\mathrm{'+self.element+'}\\ '
            s += str(self.n)+l+'_{'+str(self.j)+'}^{'+str(self.f)+'}'

        if self.m is not None:
            s = s[:-1] + ','+str(self.m)+'}'

        return s

    def __eq__(self, other):
        r"""Two states are equal if all of their identifiers are the same.

        >>> g1=State("Rb",87,5,0,1/Integer(2),1,1)
        >>> g2=State("Rb",87,5,0,1/Integer(2),2,2)
        >>> g1 == g2
        False
        """
        return self.quantum_numbers == other.quantum_numbers


class Transition(object):
    r"""This class describes a transition between different atomic states.

    For instance, the transition between the hyperfine ground states of cesium
    used in the definition of the second.

    >>> g1 = State("Cs", 133, 6, 0, 1/Integer(2),3)
    >>> g2 = State("Cs", 133, 6, 0, 1/Integer(2),4)
    >>> clock =Transition(g2,g1)
    >>> clock
    133Cs 6S_1/2^4 --/--> 133Cs 6S_1/2^3

    """

    def __init__(self, e1, e2, verbose=1):
        r'''This class describes a transition between atomic states. For
        instance, the transition between the hyperfine ground states of cesium
        used in the definition of the second.

        >>> g1 = State("Cs", 133, 6, 0, 1/Integer(2),3)
        >>> g2 = State("Cs", 133, 6, 0, 1/Integer(2),4)
        >>> clock =Transition(g2,g1)
        >>> clock
        133Cs 6S_1/2^4 --/--> 133Cs 6S_1/2^3

        Useful properties of transitions are whether they are electric-dipole
        allowed:
        >>> clock.allowed
        False

        The absolute frequency of the transition (in Hz):
        >>> print clock.nu
        9192631770.0

        The angular frequency of the transition (in rad/s):
        >>> print clock.omega
        57759008871.6

        The wavelength of the transition (in m):
        >>> print clock.wavelength
        0.0326122557175

        The Einstein A coefficient at the fine structure level (in Hz):
        >>> print clock.einsteinA
        0.0

        The states that form the transition:
        >>> clock.e1, clock.e2
        (133Cs 6S_1/2^4, 133Cs 6S_1/2^3)

        A transition between two different atoms will raise an error:
        >>> g2Rb85 = State("Rb", 85, 5, 0, 1/Integer(2), 3)
        >>> Transition(g2Rb85, g1)
        Traceback (most recent call last):
        ...
        NotImplementedError: Transition between different elements.

        A transition between two different isotopes will raise an error:
        >>> g1Rb87 = State("Rb", 87, 5, 0, 1/Integer(2), 1)
        >>> Transition(g2Rb85, g1Rb87)
        Traceback (most recent call last):
        ...
        NotImplementedError: Transition between different isotopes.

        '''
        if e1.element != e2.element:
            raise NotImplementedError, 'Transition between different elements.'
        if e1.isotope != e2.isotope:
            raise NotImplementedError, 'Transition between different isotopes.'

        self.e1 = e1; self.e2 = e2

        # We will now determine whether the transition is allowed by electric
        # dipole transition rules.
        allowed = True
        # The states must have oposite parity.
        # In particular, transitions between identical fine states are
        # forbidden.
        if e1 == e2: allowed = False

        # allowed=True
        # l1=e1.l; l2=e2.l
        # j1=e1.j; j2=e2.j

        l1 = e1.l; l2 = e2.l
        j1 = e1.j; j2 = e2.j

        # DeltaL=+-1
        if not (abs(l2-l1) in [-1, 1]): allowed = False
        # DeltaJ=0,+-1
        if not (abs(j2-j1) in [-1, 0, 1]): allowed = False

        if e1.f is not None and e2.f is not None:
            f1 = e1.f; f2 = e2.f
            # DeltaF=0,+-1
            if not (abs(f2-f1) in [-1, 0, 1]): allowed = False

            if f1 == 0 and f2 == 0: allowed = False

        if e1.m is not None and e2.m is not None:
            m1 = e1.m; m2 = e2.m
            if abs(m1) not in range(0, f1+1): allowed = False
            if abs(m2) not in range(0, f2+1): allowed = False
            if m1 == 0 and m2 == 0:
                if not (abs(f2-f1) in [-1, 1]): allowed = False
            if not (abs(m2-m1) in [-1, 0, 1]): allowed = False

        self.allowed = allowed
        ####################################################################
        self.nu = e1.nu-e2.nu
        if self.nu == 0: self.wavelength = float('inf')
        else: self.wavelength = 299792458.0 / self.nu
        self.omega = self.nu*2*pi

        # We find the Einstein A and B coefficients up to the fine structure
        # according to literature (in Hz).
        n1 = e1.n; l1 = e1.l; j1 = e1.j
        n2 = e2.n; l2 = e2.l; j2 = e2.j
        ord1 = [n1, l1, j1, n2, l2, j2]
        ord2 = [n2, l2, j2, n1, l1, j1]

        ####################################################################
        # This is the database of Einstein A coefficients.
        # The quantum numbers of the lower states followed by those of excited
        # states followed by Einstein A coefficients in 2*pi/s (angular
        # frequency).
        element = e1.element
        if element == "Rb":
            pairs = [[5, 1, 1/Integer(2), 5, 0, Integer(1)/Integer(2), 3.61031827750539e7],
                     [5, 1, 3/Integer(2), 5, 0, Integer(1)/Integer(2), 3.81075188880442e7],
                     [4, 2, 5/Integer(2), 5, 1, Integer(3)/Integer(2), 1.06814150222053e7],
                     [4, 2, 3/Integer(2), 5, 1, Integer(1)/Integer(2), 9.42477796076938e6],
                     [4, 2, 3/Integer(2), 5, 1, Integer(3)/Integer(2), 1.88495559215388e6],
                     [6, 0, 1/Integer(2), 5, 1, Integer(1)/Integer(2), 1.09704415463356e7], # The branching ratios here are not well referenced.
                     [6, 0, 1/Integer(2), 5, 1, Integer(3)/Integer(2), 1.09704415463356e7], # The branching ratios here are not well referenced.
                     [6, 1, 3/Integer(2), 5, 0, Integer(1)/Integer(2), 1.87867240684670e6],
                     [6, 1, 3/Integer(2), 4, 2, Integer(5)/Integer(2), 179699.099785336],
                     [6, 1, 3/Integer(2), 4, 2, Integer(3)/Integer(2), 1.61729189806803e6],
                     [6, 1, 3/Integer(2), 6, 0, Integer(1)/Integer(2), 4.49247749463340e6],
                     [5, 2, 5/Integer(2), 5, 1, Integer(3)/Integer(2), 3.10264947105589e6],
                     [5, 2, 5/Integer(2), 6, 1, Integer(3)/Integer(2), 1.09012008442504e6]]
        elif element == "Cs":
            # Reference [3]
            ###################################################################
            # [1] Steck
            # [2] Georgiades, Polzik, Kimble - 1994 - Two-photon spectroscopy
            # the 6S 1 2 -
            # [3] Measurement of the lifetime of the atomic cesium 5(2)D(5/2)
            # state with diode-laser excitation.
            gam = 2*pi*3.2e6
            # If the branching ratios between
            # 6D_5/2 -> 6P_3/2 and 6D_5/2 -> 7P_3/2
            # are 0.74 and 0.26 as they are for
            # 5D_5/2 -> 5P_3/2 and 5D_5/2 -> 6P_3/2 in rubidium, then
            gam1 = gam*0.74
            gam2 = gam*0.26
            ###################################################################
            pairs = [[6, 1, Integer(1)/Integer(2), 6, 0, Integer(1)/Integer(2), 2*pi*4.575e6],#[1]
                     [6, 1, Integer(3)/Integer(2), 6, 0, Integer(1)/Integer(2), 2*pi*5.234e6],#[1]

                     [6, 2, Integer(5)/Integer(2), 7, 1, Integer(3)/Integer(2), gam2],#[2]
                     [6, 2, Integer(5)/Integer(2), 6, 1, Integer(3)/Integer(2), gam1]]#[2]

        self.einsteinA = 0.0
        # if ord1 == ord2: self.einsteinA=0.0
        if self.allowed: self.einsteinA = None
        for pair in pairs:
            if pair[:-1] == ord1 or pair[:-1] == ord2:
                self.einsteinA = pair[-1]

        # print self.allowed, self.einsteinA
        if self.allowed and self.einsteinA is None and verbose > 0:
            s = "Warning: the transition" + str(self)
            s += "is electric dipole-allowed, but the Einstein A is not in "
            s += "the database."
            print s
            self.einsteinA = "Unknown, but different from zero!"

    def __repr__(self):
        r"""The representation routine for transitions.

        >>> g1 = State("Cs", 133, 6, 0, 1/Integer(2),3)
        >>> g2 = State("Cs", 133, 6, 0, 1/Integer(2),4)
        >>> Transition(g2,g1).__repr__()
        '133Cs 6S_1/2^4 --/--> 133Cs 6S_1/2^3'

        """

        if self.allowed:
            return self.e1.__repr__()+' -----> '+self.e2.__repr__()
        elif not self.allowed:
            return self.e1.__repr__()+' --/--> '+self.e2.__repr__()
        else:
            return self.e1.__repr__()+' --?--> '+self.e2.__repr__()

    def _latex_(self):
        r"""The representation routine for transitions.

        >>> g1 = State("Cs", 133, 6, 0, 1/Integer(2),3)
        >>> g2 = State("Cs", 133, 6, 0, 1/Integer(2),4)
        >>> Transition(g2,g1)._latex_()
        '^{133}\\mathrm{Cs}\\ 6S_{1/2}^{4}\\ \\nrightarrow \\ ^{133}\\mathrm{Cs}\\ 6S_{1/2}^{3}'

        """
        if self.allowed:
            return self.e1._latex_()+'\\ \\rightarrow \\ '+self.e2._latex_()
        elif not self.allowed:
            return self.e1._latex_()+'\\ \\nrightarrow \\ '+self.e2._latex_()
        else:
            return self.e1._latex_()+'\\ \\rightarrow^? \\ '+self.e2._latex_()

        return self.e1._latex_()+'\\ \\nleftrightarrow \\ '+self.e2._latex_()

    def __eq__(self, other):
        r"""The equals routine for transitions.

        Two transitions are equal if their constituent states are equal:
        >>> g1 = State("Cs", 133, 6, 0, 1/Integer(2),3)
        >>> g2 = State("Cs", 133, 6, 0, 1/Integer(2),4)
        >>> t1=Transition(g2,g1)._latex_()
        >>> G1 = State("Rb",  85, 5, 0, 1/Integer(2),2)
        >>> G2 = State("Rb",  85, 5, 0, 1/Integer(2),3)
        >>> t2=Transition(G2,G1)._latex_()
        >>> t1 == t2
        False

        """
        return self.e1 == other.e1 and self.e2 == other.e2


def split_fine_to_hyperfine(state):
    if type(state) == list:
        mag = []
        for s in state:
            mag += split_fine_to_hyperfine(s)
        return mag

    if len(state.quantum_numbers) != 4:
        s = str(state)+' no es un estado fino.'
        raise ValueError, s

    return [State(state.element, state.isotope,
            state.n, state.l, state.j, f) for f in state.fperm]

def split_hyperfine_to_magnetic(state):
    if type(state) == list:
        mag = []
        for s in state:
            mag += split_hyperfine_to_magnetic(s)
        return mag

    if len(state.quantum_numbers) != 5:
        raise ValueError, str(state)+' is not a hyperfine state.'

    return [State(state.element, state.isotope,
            state.n, state.l, state.j, state.f, m) for m in state.mperm]

def split_fine_to_magnetic(state):
    if type(state) == list:
        mag = []
        for s in state:
            mag += split_fine_to_magnetic(s)
        return mag

    if len(state.quantum_numbers) != 4:
        raise ValueError, str(state)+' is not a fine state.'

    hip = split_fine_to_hyperfine(state)

    mag = []
    for ll in [split_hyperfine_to_magnetic(h) for h in hip]:
        mag += ll
    return mag

def order_by_energy(states):
    iso = states[0].isotope
    for ss in states:
        if ss.isotope != iso:
            raise ValueError, 'We have a wrong isotope in this list: '+str(ss)
    aux = [(ee.nu, ee) for ee in states]
    aux.sort()

    return [i[1] for i in aux]


def make_list_of_states(states, structure=None, verbose=1):
    if structure is None:
        return order_by_energy(states)
    elif structure == 'hyperfine':
        l1 = order_by_energy(states)
        l2 = split_fine_to_hyperfine(l1)
        l3 = order_by_energy(l2)
        if l2 != l3:
            if verbose > 0:
                s = 'Warning: the ordering of the hyperfine states has been'
                s += ' changed to ensure they are ordered by ascending energy.'
                print s
        return l3
    elif structure == 'magnetic':
        l1 = order_by_energy(states)
        l2 = split_fine_to_hyperfine(l1)
        l3 = order_by_energy(l2)
        if l2 != l3:
            if verbose > 0:
                s = 'Warning: the ordering of the hyperfine states has been'
                s += ' changed to ensure they are ordered by ascending energy.'
                print
        l4 = split_hyperfine_to_magnetic(l3)
        return l4


def calculate_omega_matrix(states, Omega=1):
    """Calculate the matrix of transition frequencies.

    This function recieves a list of states and returns the corresponding
    omega_ij matrix, rescaled to units of Omega. These are understood to be
    absolute frequencies (as opposed angular frequencies).
    """
    N = len(states)
    omega = [[2*Pi*(states[i].nu-states[j].nu)/Omega
              for j in range(N)] for i in range(N)]

    return omega


def get_einstein_A_matrix(fine_states,Omega=1):
    einsteinA = [[0.0 for jj in range(len(fine_states))]
                 for ii in range(len(fine_states))]

    for ii in range(len(fine_states)):
        i = fine_states[ii]
        for jj in range(ii):
            j = fine_states[jj]
            t = Transition(i, j)
            Aij = t.einsteinA

            if Aij is None:
                s = 'The Einstein coeficcient A_ij between '+str(i)
                s += ' and '+str(j)+' is unknown, '
                s += 'please add it to the database in the code of the class '
                s += '"Transition" to proceed.'
                raise NotImplementedError, s
            einsteinA[ii][jj] = Aij/Omega
            einsteinA[jj][ii] = -Aij/Omega

    return einsteinA


def calculate_gamma_matrix(magnetic_states, Omega=1):
    r"""Calculate the matrix of decay between states.

    This function calculates the matrix $\gamma_{ij}$ of decay rates between
    states |i> and |j> (in the units specified by the Omega argument).

    >>> g=State("Rb",87,5,0,1/Integer(2))
    >>> e=State("Rb",87,5,1,3/Integer(2))
    >>> magnetic_states=make_list_of_states([g,e],"magnetic")

    To return the rates in rad/s:
    >>> print calculate_gamma_matrix(magnetic_states)
    [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -12702506.296014734, -15878132.870018415, -15878132.870018415, -0.0, -19053759.4440221, -9526879.72201105, -3175626.5740036834, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -12702506.296014734, -15878132.870018415, -0.0, -15878132.870018415, -0.0, -9526879.72201105, -12702506.296014734, -9526879.72201105, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -12702506.296014734, -0.0, -15878132.870018415, -15878132.870018415, -0.0, -0.0, -3175626.5740036834, -9526879.72201105, -19053759.4440221, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -3810751.8888044204, -0.0, -0.0, -12702506.296014734, -6351253.148007367, -0.0, -0.0, -0.0, -38107518.8880442, -12702506.296014734, -2540501.2592029474, -0.0, -0.0, -0.0, -0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -1905375.9444022102, -1905375.9444022102, -0.0, -6351253.148007367, -3175626.5740036834, -9526879.72201105, -0.0, -0.0, -0.0, -25405012.592029467, -20324010.07362358, -7621503.77760884, -0.0, -0.0, -0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -635125.3148007367, -2540501.259202947, -635125.3148007367, -0.0, -9526879.72201105, -0.0, -9526879.72201105, -0.0, -0.0, -0.0, -15243007.55521768, -22864511.33282652, -15243007.55521768, -0.0, -0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -1905375.9444022102, -1905375.9444022102, -0.0, -0.0, -9526879.72201105, -3175626.5740036834, -6351253.148007367, -0.0, -0.0, -0.0, -7621503.77760884, -20324010.07362358, -25405012.592029467, -0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -3810751.8888044204, -0.0, -0.0, -0.0, -6351253.148007367, -12702506.296014734, -0.0, -0.0, -0.0, -0.0, -2540501.2592029474, -12702506.296014734, -38107518.8880442], [12702506.296014734, 12702506.296014734, 12702506.296014734, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [15878132.870018415, 15878132.870018415, 0.0, 3810751.8888044204, 1905375.9444022102, 635125.3148007367, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [15878132.870018415, 0.0, 15878132.870018415, 0.0, 1905375.9444022102, 2540501.259202947, 1905375.9444022102, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 15878132.870018415, 15878132.870018415, 0.0, 0.0, 635125.3148007367, 1905375.9444022102, 3810751.8888044204, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [19053759.4440221, 0.0, 0.0, 12702506.296014734, 6351253.148007367, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [9526879.72201105, 9526879.72201105, 0.0, 6351253.148007367, 3175626.5740036834, 9526879.72201105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [3175626.5740036834, 12702506.296014734, 3175626.5740036834, 0.0, 9526879.72201105, 0.0, 9526879.72201105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 9526879.72201105, 9526879.72201105, 0.0, 0.0, 9526879.72201105, 3175626.5740036834, 6351253.148007367, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 19053759.4440221, 0.0, 0.0, 0.0, 6351253.148007367, 12702506.296014734, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 38107518.8880442, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 12702506.296014734, 25405012.592029467, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 2540501.2592029474, 20324010.07362358, 15243007.55521768, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 7621503.77760884, 22864511.33282652, 7621503.77760884, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 15243007.55521768, 20324010.07362358, 2540501.2592029474, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 25405012.592029467, 12702506.296014734, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 38107518.8880442, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]

    To return the rates in 10^6 rad /s:
    >>> gamma = calculate_gamma_matrix(magnetic_states,Omega=1e6)
    >>> gamma
    [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -12.702506296014734, -15.878132870018417, -15.878132870018417, -0.0, -19.053759444022102, -9.526879722011051, -3.1756265740036835, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -12.702506296014734, -15.878132870018417, -0.0, -15.878132870018417, -0.0, -9.526879722011051, -12.702506296014734, -9.526879722011051, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -12.702506296014734, -0.0, -15.878132870018417, -15.878132870018417, -0.0, -0.0, -3.1756265740036835, -9.526879722011051, -19.053759444022102, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -3.8107518888044205, -0.0, -0.0, -12.702506296014734, -6.351253148007367, -0.0, -0.0, -0.0, -38.107518888044204, -12.702506296014734, -2.5405012592029474, -0.0, -0.0, -0.0, -0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -1.9053759444022103, -1.9053759444022103, -0.0, -6.351253148007367, -3.1756265740036835, -9.526879722011051, -0.0, -0.0, -0.0, -25.40501259202947, -20.32401007362358, -7.62150377760884, -0.0, -0.0, -0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.6351253148007368, -2.540501259202947, -0.6351253148007368, -0.0, -9.526879722011051, -0.0, -9.526879722011051, -0.0, -0.0, -0.0, -15.24300755521768, -22.86451133282652, -15.24300755521768, -0.0, -0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -1.9053759444022103, -1.9053759444022103, -0.0, -0.0, -9.526879722011051, -3.1756265740036835, -6.351253148007367, -0.0, -0.0, -0.0, -7.62150377760884, -20.32401007362358, -25.40501259202947, -0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -3.8107518888044205, -0.0, -0.0, -0.0, -6.351253148007367, -12.702506296014734, -0.0, -0.0, -0.0, -0.0, -2.5405012592029474, -12.702506296014734, -38.107518888044204], [12.702506296014734, 12.702506296014734, 12.702506296014734, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [15.878132870018417, 15.878132870018417, 0.0, 3.8107518888044205, 1.9053759444022103, 0.6351253148007368, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [15.878132870018417, 0.0, 15.878132870018417, 0.0, 1.9053759444022103, 2.540501259202947, 1.9053759444022103, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 15.878132870018417, 15.878132870018417, 0.0, 0.0, 0.6351253148007368, 1.9053759444022103, 3.8107518888044205, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [19.053759444022102, 0.0, 0.0, 12.702506296014734, 6.351253148007367, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [9.526879722011051, 9.526879722011051, 0.0, 6.351253148007367, 3.1756265740036835, 9.526879722011051, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [3.1756265740036835, 12.702506296014734, 3.1756265740036835, 0.0, 9.526879722011051, 0.0, 9.526879722011051, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 9.526879722011051, 9.526879722011051, 0.0, 0.0, 9.526879722011051, 3.1756265740036835, 6.351253148007367, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 19.053759444022102, 0.0, 0.0, 0.0, 6.351253148007367, 12.702506296014734, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 38.107518888044204, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 12.702506296014734, 25.40501259202947, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 2.5405012592029474, 20.32401007362358, 15.24300755521768, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 7.62150377760884, 22.86451133282652, 7.62150377760884, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 15.24300755521768, 20.32401007362358, 2.5405012592029474, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 25.40501259202947, 12.702506296014734, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 38.107518888044204, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]

    Let us test if all D2 lines decay at the expected rate (6.065 MHz):
    >>> Gamma =[ sum([ gamma[i][j] for j in range(i)])/2/pi for i in range(len(magnetic_states))][8:]
    >>> for Gammai in Gamma: print Gammai
    6.065
    6.065
    6.065
    6.065
    6.065
    6.065
    6.065
    6.065
    6.065
    6.065
    6.065
    6.065
    6.065
    6.065
    6.065
    6.065

    """
    Ne = len(magnetic_states)

    II = magnetic_states[0].i

    gamma = [[0.0 for j in range(Ne)] for i in range(Ne)]
    for i in range(Ne):
        for j in range(i):
            ei = magnetic_states[i]
            ej = magnetic_states[j]
            einsteinAij = Transition(ei, ej).einsteinA

            if einsteinAij != 0:
                ji = ei.j; jj = ej.j
                fi = ei.f; fj = ej.f
                mi = ei.m; mj = ej.m

                gammaij = (2.0*ji+1)
                gammaij *= (2.0*fi+1)
                gammaij *= (2.0*fj+1)

                gammaij *= float(wigner_6j(ji, fi, II, fj, jj, 1)**2)
                gammaij *= sum([float(wigner_3j(fj, 1, fi, -mj, q, mi)**2)
                                for q in [-1, 0, 1]])

                gammaij *= einsteinAij/Omega
                gammaij = float(gammaij)

                gamma[i][j] = gammaij
                gamma[j][i] = -gammaij

    return gamma


def calculate_reduced_matrix_elements(fine_states):
    '''Calculate the reduced matrix elements for a list of fine states.

    This function calculates the reduced matrix elments
                <N,L,J||T^1(r)||N',L',J'>
    given a list of fine states.

    >>> g=State("Rb",87,5,0,1/Integer(2))
    >>> e1=State("Rb",87,5,1,1/Integer(2))
    >>> e2=State("Rb",87,5,1,3/Integer(2))
    >>> red=calculate_reduced_matrix_elements([g,e1,e2])
    >>> print red[1][0]
    4.23143658816481
    >>> print red[2][0]
    8.45396724128841

    '''
    # We calculate the reduced matrix elements starting from the list of
    # fine_states. The factor composed of physical quantities.
    factor = sqrt(3*c**3*me**2*e**2/(16*Pi*epsilon0*hbar**3))
    # We read the database to obtain the Einstein A coefficients in Hz.
    einsteinA = get_einstein_A_matrix(fine_states)
    # We read the database to obtain the transition frequencies in Hz.
    omega_fine = calculate_omega_matrix(fine_states)

    reduced_matrix_elements = [[0.0 for jj in range(len(fine_states))]
                               for ii in range(len(fine_states))]

    for ii in range(len(fine_states)):
        i = fine_states[ii]
        for jj in range(ii):
            j = fine_states[jj]
            t = Transition(i, j)
            einsteinAij = einsteinA[ii][jj]
            omega0 = omega_fine[ii][jj]

            # The formula is valid only for i =/= j so that omega0 =/= 0.
            # Because fine states are asumed to be ordered by their energies we
            # can asume that i decays in j.
            Ji = i.j; Jj = j.j

            rij = (2.0*Ji+1)/sqrt(2.0*Jj+1)*sqrt(einsteinAij/omega0**3)
            rij = factor*rij

            reduced_matrix_elements[ii][jj] = rij
            # We add the matrix elements on the other side of the diagonal.
            reduced_matrix_elements[jj][ii] = rij*(-1)**(Ji-Jj)

    return reduced_matrix_elements


def calculate_r_matrices(fine_states, reduced_matrix_elements):
    magnetic_states = make_list_of_states(fine_states, 'magnetic', verbose=0)
    index_list_fine, index_list_hyperfine = calculate_boundaries(fine_states, magnetic_states)

    Ne = len(magnetic_states)

    r = [[[0.0 for j in range(Ne)] for i in range(Ne)] for p in range(3)]

    II = fine_states[0].i

    for p in [-1, 0, 1]:
        for i in range(Ne):
            ei = magnetic_states[i]
            ii = fine_index(i, index_list_fine)

            for j in range(Ne):
                ej = magnetic_states[j]
                jj = fine_index(j, index_list_fine)

                reduced_matrix_elementij = reduced_matrix_elements[ii][jj]
                if reduced_matrix_elementij != 0:

                    ji = ei.j; jj = ej.j
                    fi = ei.f; fj = ej.f
                    mi = ei.m; mj = ej.m

                    rpij = (-1)**(fi-mi)
                    rpij *= wigner_3j(fi, 1, fj, -mi, p, mj)

                    rpij *= (-1)**(fj+ji+1+II)
                    rpij *= sqrt(2*fj+1)
                    rpij *= sqrt(2*fi+1)
                    rpij *= wigner_6j(ji, jj, 1, fj, fi, II)

                    rpij *= reduced_matrix_elementij

                    r[p+1][i][j] = float(rpij)
    return r


def calculate_matrices(states, Omega=1):
    r"""Calculate the matrices omega_ij, gamma_ij, r_pij.

    This function calculates the matrices omega_ij, gamma_ij and r_pij given a
    list of atomic states. The states can be arbitrarily in their fine,
    hyperfine or magnetic detail.
    """
    # We check that all states belong to the same element and the same isotope.
    iso = states[0].isotope
    element = states[0].element
    for state in states[1:]:
        if state.element != element:
            raise ValueError, 'All states must belong to the same element.'
        if state.isotope != iso:
            raise ValueError, 'All states must belong to the same isotope.'

    # We find the fine states involved in the problem.
    fine_states = find_fine_states(states)

    # We find the full magnetic states. The matrices will be first calculated
    # for the complete problem and later reduced to include only the states of
    # interest.
    full_magnetic_states = make_list_of_states(fine_states, 'magnetic',
                                               verbose=0)

    # We calculate the indices corresponding to each sub matrix of fine and
    # hyperfine levels.

    # We calculate the frequency differences between states.
    omega_full = calculate_omega_matrix(full_magnetic_states, Omega)
    # We calculate the matrix gamma.

    gamma_full = calculate_gamma_matrix(full_magnetic_states, Omega)
    # We calculate the reduced matrix elements
    reduced_matrix_elements = calculate_reduced_matrix_elements(fine_states)

    # We calculate the matrices r_-1, r_0, r_1
    r_full = calculate_r_matrices(fine_states, reduced_matrix_elements)

    # Reduction to be implemented
    omega = omega_full
    r = r_full
    gamma = gamma_full
    return omega, gamma, r


def calculate_boundaries(fine_states, full_magnetic_states):
    r"""Calculate the boundary indices within a list of magnetic states.

    This function calculates the boundary indices of each fine state
    and each hyperfine state within a list of magnetic states. The output
    is a list of tuples (a,b) with a the starting index of a state and
    b it's ending index.

    >>> g=State("Rb", 87, 5, 0, 1/Integer(2))
    >>> full_magnetic_states=make_list_of_states([g],"magnetic")
    >>> calculate_boundaries([g], full_magnetic_states)
    ([(0, 8)], [(0, 3), (3, 8)])

    """
    N_magnetic = len(full_magnetic_states)

    # We calculate the boundaries of the various detail levels
    # First we will make a list of indices of that will tell where each fine
    # level begins and ends.
    fq = full_magnetic_states[0].quantum_numbers[:4]
    index_list_fine = []; start_fine = 0
    # And another list of indices that will tell where each hyperfine level
    # begins and ends.
    hq = full_magnetic_states[0].quantum_numbers[:5]
    index_list_hyperfine = []; start_hyperfine = 0

    for i in range(N_magnetic):
        magnetic = full_magnetic_states[i]
        if magnetic.quantum_numbers[:4] != fq:
            index_list_fine += [(start_fine, i)]
            start_fine = i
            fq = magnetic.quantum_numbers[:4]

        if magnetic.quantum_numbers[:5] != hq:
            index_list_hyperfine += [(start_hyperfine, i)]
            start_hyperfine = i
            hq = magnetic.quantum_numbers[:5]

        if i == N_magnetic-1:
            index_list_fine += [(start_fine, i+1)]
            index_list_hyperfine += [(start_hyperfine, i+1)]
    return index_list_fine, index_list_hyperfine

def fine_index(magnetic_index,index_list_fine):
	""""""
	N_fine=len(index_list_fine)
	for i in range(N_fine):
		if index_list_fine[i][0] <= magnetic_index < index_list_fine[i][1]:
			return i


def quaver(isotope, p, J, F, M, Jp, Fp, Mp, numeric=False, verbose=False):
    if isotope == 85:
        II = Integer(5)/Integer(2)
    elif isotope == 87:
        II = Integer(3)/Integer(2)

    qu = (-1)**(F-M+J+II+Fp+1)
    qu *= sqrt((2*F+1)*(2*Fp+1))
    if not sage_included:
        II = float(II)
        J = float(J)
        Jp = float(Jp)
    qu *= wigner_3j(F, 1, Fp, -M, p, Mp)*wigner_6j(II, J, F, 1, Fp, Jp)

    if numeric: return float(qu)
    else: return qu


def find_fine_states(magnetic_states):
    fine_states = []
    for state in magnetic_states:
        fq = state.quantum_numbers[:4]
        fine_state = State(state.element, fq[0], fq[1], fq[2], fq[3])
        if fine_state not in fine_states: fine_states += [fine_state]
    return fine_states


def exclude_states(omega, gamma, r, Lij, states, excluded_states):
    """Exclude states from matrices.

    This function takes the matrices and excludes the states listed in
    excluded_states.
    """
    Ne = len(omega)
    excluded_indices = [i for i in range(Ne) if states[i] in excluded_states]

    omega_new = []; gamma_new = []; r_new = [[], [], []]; Lij_new = []
    for i in range(Ne):
        row_om = []; row_ga = []; row_L = []
        for j in range(Ne):
            if j not in excluded_indices:
                row_om += [omega[i][j]]
                row_ga += [gamma[i][j]]
                row_L += [Lij[i][j]]
        if i not in excluded_indices:
            omega_new += [row_om]
            gamma_new += [row_ga]
            Lij_new += [row_L]

    for p in range(3):
        for i in range(Ne):
            row_r = []
            for j in range(Ne):
                if j not in excluded_indices:
                    row_r += [r[p][i][j]]
            if i not in excluded_indices:
                r_new[p] += [row_r]

    states_new = [states[i] for i in range(Ne) if i not in excluded_indices]

    return omega_new, gamma_new, r_new, Lij_new, states_new

def reduce_magnetic_to_hyperfine(omega,gamma,r,Lij,magnetic_states,hyperfine_states,isotropic_r=False):
	#We find the fine states involved in the problem.
	fine_states=find_fine_states(magnetic_states)

	#We calculate the indices corresponding to each sub matrix of fine and hyperfine levels.
	index_list_fine,index_list_hyperfine=calculate_boundaries(fine_states,magnetic_states)

	#We determine which indices will be reduced
	index_list_to_reduce0=[]
	for a,b in index_list_hyperfine:
		qn = magnetic_states[a].quantum_numbers[:5]
		for hs in hyperfine_states:
			if hs.quantum_numbers == qn:
				index_list_to_reduce0+=[(a,b)]

	Ne_magnetic=len(magnetic_states)
	index_list_to_reduce=[]
	for i in range(Ne_magnetic):
		band=True
		for a,b in index_list_to_reduce0:
			if a <= i < b:
				if (a,b) not in index_list_to_reduce:
					index_list_to_reduce+=[(a,b)]
				band=False
				break
		if band:
			index_list_to_reduce+=[(i,i+1)]

	#We calculate the x,y,z matrices in order to reduce them.
	x=[[   (r[0][i][j]-r[2][i][j])/sqrt(2.0) for j in range(Ne_magnetic)] for i in range(Ne_magnetic)]
	y=[[1j*(r[0][i][j]+r[2][i][j])/sqrt(2.0) for j in range(Ne_magnetic)] for i in range(Ne_magnetic)]
	z=[[r[1][i][j]                           for j in range(Ne_magnetic)] for i in range(Ne_magnetic)]

	#We make the reduction
	Ne_reduced = len(index_list_to_reduce)
	omega_red=[]; gamma_red=[]; xyz=[[],[],[]]; Lij_red=[]; states_red=[]
	for ii in range(Ne_reduced):
		aii,bii=index_list_to_reduce[ii]
		omega_row=[]; gamma_row=[]; x_row,y_row,z_row=[],[],[]; Lij_row=[]
		for jj in range(Ne_reduced):
			ajj,bjj=index_list_to_reduce[jj]

			giijj=0; xiijj,yiijj,ziijj=0,0,0; Liijj=[]
			for i in range(aii,bii):
				for j in range(ajj,bjj):
					#We can simply take whichever submatrix element for omega.
					oiijj=omega[i][j]
					#We sum the elements of gamma.
					giijj+=gamma[i][j]
					#We make the quadrature sum the elements of x,y,z.
					xiijj+=x[i][j].conjugate()*x[i][j]
					yiijj+=y[i][j].conjugate()*y[i][j]
					ziijj+=z[i][j].conjugate()*z[i][j]
					#For Lij we include whichever l exists in the submatrix in the reduced Lij.
					for l in Lij[i][j]:
						if l not in Liijj: Liijj+=[l]

			omega_row+=[oiijj]
			gamma_row+=[giijj]
			x_row+=[sqrt(xiijj.real)]
			y_row+=[sqrt(yiijj.real)]
			z_row+=[sqrt(ziijj.real)]
			Lij_row+=[Liijj]
		#For the states we check wether there is a reduction, and add the apropiate states.
		element_n=magnetic_states[aii].element
        qn=magnetic_states[aii].quantum_numbers[:5]
        hs=State(element_n,qn[0], qn[1], qn[2], qn[3], qn[4])
        if bii-aii==1:
			if split_hyperfine_to_magnetic([hs])==[magnetic_states[aii]]:
				if hs in hyperfine_states:
					states_red+=[hs]
				else:
					states_red+=[magnetic_states[aii]]
			else:
				states_red+=[magnetic_states[aii]]
        else:
			qn=magnetic_states[aii].quantum_numbers[:5]
			element_n=magnetic_states[aii].element
			states_red+=[State(element_n,qn[0], qn[1], qn[2], qn[3], qn[4])]
        omega_red+=[omega_row]
        gamma_red+=[gamma_row]
        xyz[0]+=[x_row]; xyz[1]+=[y_row]; xyz[2]+=[z_row]
        Lij_red+=[Lij_row]

	#We calculate the reduced r from xyz
	if isotropic_r:
		#We can impose the isotropy of r
		xyz[0]=xyz[2][:]; xyz[1]=xyz[2][:]
	rm1_red=[[ ( xyz[0][i][j]-1j*xyz[1][i][j])/sqrt(2.0) for j in range(Ne_reduced)] for i in range(Ne_reduced)]
	r0_red = xyz[2]
	rp1_red=[[ (-xyz[0][i][j]-1j*xyz[1][i][j])/sqrt(2.0) for j in range(Ne_reduced)] for i in range(Ne_reduced)]
	r_red=[rm1_red, r0_red, rp1_red]

	return omega_red,gamma_red,r_red,Lij_red,states_red

def calculate_reduced_matrix_elements_0(fine_states):
	'''This function calculates the reduced matrix elments <N,L,J||T^1(r)||N',L',J'> given a list of fine states.'''
	#We calculate the reduced matrix elements starting from the list of fine_states

	#The factor composed of physical quantities.
	factor=sqrt(3*c**3*me**2*e**2/(16*Pi*epsilon0*hbar**3))
	#We read the database to obtain the Einstein A coefficients in Hz.
	einsteinA=get_einstein_A_matrix(fine_states)
	#We read the database to obtain the transition frequencies in Hz.
	omega_fine=calculate_omega_matrix(fine_states)

	reduced_matrix_elements=[[0.0 for jj in range(len(fine_states))] for ii in range(len(fine_states))]

	for ii in range(len(fine_states)):
		i=fine_states[ii]
		for jj in range(ii):
			j=fine_states[jj]
			t=Transition(i,j)
			einsteinAij=einsteinA[ii][jj]
			omega0=omega_fine[ii][jj]

			#The formula is valid only for i =/= j so that omega0 =/= 0.
			#Because fine states are asumed to be ordered by their energies we can asume that
			# i decays in j.
			Ji=i.j; Jj=j.j

			rij=sqrt((2.0*Ji+1)/(2*Jj+1))*sqrt(einsteinAij/omega0**3)
			rij=factor*rij

			reduced_matrix_elements[ii][jj]=rij
			#We add the matrix elements on the other side of the diagonal.
			reduced_matrix_elements[jj][ii]=rij*(-1)**(Jj-Ji)

	return reduced_matrix_elements

def calculate_reduced_matrix_elements_steck(fine_states):
	'''This function calculates the reduced matrix elments <N,L,J||T^1(r)||N',L',J'> given a list of fine states.'''
	#We calculate the reduced matrix elements starting from the list of fine_states

	#The factor composed of physical quantities.
	factor=sqrt(3*c**3*me**2*e**2/(16*Pi*epsilon0*hbar**3))
	#We read the database to obtain the Einstein A coefficients in Hz.
	einsteinA=get_einstein_A_matrix(fine_states)
	#We read the database to obtain the transition frequencies in Hz.
	omega_fine=calculate_omega_matrix(fine_states)

	reduced_matrix_elements=[[0.0 for jj in range(len(fine_states))] for ii in range(len(fine_states))]

	for ii in range(len(fine_states)):
		i=fine_states[ii]
		for jj in range(ii):
			j=fine_states[jj]
			t=Transition(i,j)
			einsteinAij=einsteinA[ii][jj]
			omega0=omega_fine[ii][jj]

			#The formula is valid only for i =/= j so that omega0 =/= 0.
			#Because fine states are asumed to be ordered by their energies we can asume that
			# i decays in j.
			Ji=i.j; Jj=j.j

			rij=sqrt((2.0*Ji+1)/(2*Jj+1))*sqrt(einsteinAij/omega0**3)
			rij=factor*rij

			reduced_matrix_elements[ii][jj]=rij
			#We add the matrix elements on the other side of the diagonal.
			reduced_matrix_elements[jj][ii]=rij*(-1)**(Jj-Ji) * sqrt(2.0*Ji+1)/sqrt(2.0*Jj+1)

	return reduced_matrix_elements

def calculate_r_matrices_0(fine_states, reduced_matrix_elements):

	full_magnetic_states=make_list_of_states(fine_states,'magnetic',verbose=0)
	index_list_fine, index_list_hyperfine = calculate_boundaries(fine_states, full_magnetic_states)

	N_magnetic=len(full_magnetic_states)
	r=[]
	for p in [-1,0,1]:
		mat=[]
		for i in range(N_magnetic):
			row=[]
			ISO,N,L,J,F,M=full_magnetic_states[i].quantum_numbers
			ii=fine_index(i, index_list_fine)
			for j in range(N_magnetic):
				ISOp,Np,Lp,Jp,Fp,Mp=full_magnetic_states[j].quantum_numbers
				jj=fine_index(j, index_list_fine)

				red=reduced_matrix_elements[ii][jj]
				qu=quaver(ISO,p, J, F, M, Jp, Fp, Mp,numeric=True)
				row+=[red*qu]
			mat+=[row]
		r+=[mat]
	return r

def calculate_r_matrices_steck(fine_states, reduced_matrix_elements):

	magnetic_states=make_list_of_states(fine_states,'magnetic',verbose=0)
	index_list_fine, index_list_hyperfine = calculate_boundaries(fine_states, magnetic_states)

	Ne=len(magnetic_states)

	r=[[[0.0 for j in range(Ne)] for i in range(Ne)] for p in range(3)]

	if fine_states[0].isotope==85:
		II=Integer(5)/Integer(2)
	elif fine_states[0].isotope==87:
		II=Integer(3)/Integer(2)

	for p in [-1,0,1]:
		for i in range(Ne):
			ei=magnetic_states[i]
			ii=fine_index(i, index_list_fine)

			for j in range(Ne):
				ej=magnetic_states[j]
				jj=fine_index(j, index_list_fine)

				reduced_matrix_elementij=reduced_matrix_elements[ii][jj]
				if reduced_matrix_elementij != 0:

					ji=ei.j; jj=ej.j
					fi=ei.f; fj=ej.f
					mi=ei.m; mj=ej.m

					rpij =(-1)**(fj-1-mi)
					rpij*=sqrt(2*fi+1)
					rpij*=wigner_3j(fj,1,fi,mj,p,-mi)

					rpij*=(-1)**(fj+ji+1+II)
					rpij*=sqrt(2*fj+1)
					rpij*=sqrt(2*ji +1)
					rpij*=wigner_6j(ji,jj,1,fj,fi,II)

					rpij*=reduced_matrix_elementij

					r[p+1][i][j]=float(rpij)
	return r


r"""
    >>> print vapour_pressure(25.0 + 273.15,"Rb")
    5.31769896107e-05
    >>> print vapour_pressure(39.3 + 273.15,"Rb")
    0.000244249795696
    >>> print vapour_pressure(90.0 + 273.15,"Rb")
    0.0155963687128
    >>> print vapour_pressure(25.0 + 273.15,"Cs")
    0.000201461144963
    >>> print vapour_pressure(28.5 + 273.15,"Cs")
    0.000297898928349
    >>> print vapour_pressure(90.0 + 273.15,"Cs")
    0.0421014384667

    The element must be in the database.

    >>> print vapour_pressure(90.0 + 273.15,"Ca")
    Traceback (most recent call last):
    ...
    ValueError: Ca is not an element in the database for this function.

    References:
    [1] Daniel A. Steck, “Cesium D Line Data,” available online at
    http://steck.us/alkalidata (revision 2.1.4, 23 December 2010).
    [2] Daniel A. Steck, “Rubidium 85 D Line Data,” available online at
    http://steck.us/alkalidata (revision 2.1.5, 19 September 2012).
    [3] Daniel A. Steck, “Rubidium 87 D Line Data,” available online at
    http://steck.us/alkalidata (revision 2.1.5, 19 September 2012).
"""


def vapour_pressure(Temperature, element):
    r"""Return the vapour pressure of rubidium or cesium in Pascals.

    This function receives as input the temperature in Kelvins and the
    name of the element.

    >>> print vapour_pressure(25.0 + 273.15,"Rb")
    5.31769896107e-05
    >>> print vapour_pressure(39.3 + 273.15,"Rb")
    0.000244249795696
    >>> print vapour_pressure(90.0 + 273.15,"Rb")
    0.0155963687128
    >>> print vapour_pressure(25.0 + 273.15,"Cs")
    0.000201461144963
    >>> print vapour_pressure(28.5 + 273.15,"Cs")
    0.000297898928349
    >>> print vapour_pressure(90.0 + 273.15,"Cs")
    0.0421014384667

    The element must be in the database.

    >>> print vapour_pressure(90.0 + 273.15,"Ca")
    Traceback (most recent call last):
    ...
    ValueError: Ca is not an element in the database for this function.

    References:
    [1] Daniel A. Steck, "Cesium D Line Data," available online at
        http://steck.us/alkalidata (revision 2.1.4, 23 December 2010).
    [2] Daniel A. Steck, "Rubidium 85 D Line Data," available online at
        http://steck.us/alkalidata (revision 2.1.5, 19 September 2012).
    [3] Daniel A. Steck, "Rubidium 87 D Line Data," available online at
        http://steck.us/alkalidata (revision 2.1.5, 19 September 2012).
    """
    if element == "Rb":
        Tmelt = 39.30+273.15  # K.
        if Temperature < Tmelt:
            P = 10**(2.881+4.857-4215.0/Temperature)  # Torr.
        else:
            P = 10**(2.881+4.312-4040.0/Temperature)  # Torr.
    elif element == "Cs":
        Tmelt = 28.5 + 273.15  # K.
        if Temperature < Tmelt:
            P = 10**(2.881+4.711-3999.0/Temperature)  # Torr.
        else:
            P = 10**(2.881+4.165-3830.0/Temperature)  # Torr.
    else:
        s = str(element)
        s += " is not an element in the database for this function."
        raise ValueError(s)

    P = P * 101325.0/760.0  # Pascals.
    return P


def vapour_number_density(Temperature, element):
    r"""Return the number of atoms in a rubidium or cesium vapour in m^-3.

    It receives as input the temperature in Kelvins and the
    name of the element.

    >>> print vapour_number_density(90.0 + 273.15,"Cs")
    8.39706962725e+18

    """
    return vapour_pressure(Temperature, element)/k_B/Temperature


def vapour_density(Temperature,element,isotope=None):
    r"""This function returns the density in a rubidium or cesium
    vapour in kg/m^-3. It receives as input the temperature in Kelvins, the
    name of the element, and optionally the isotope. If no isotope is
    specified, the density of a vapour with the natural abundances will
    be returned.

    >>> print vapour_density(90.0 + 273.15,"Cs",133)
    1.85318869181e-06

    If no isotope is specified, the natural abundances are used to calculate
    the density.

    >>> print vapour_density(25.0 + 273.15,"Rb")
    1.83339788085e-09

    """

    atom=Atom(element,isotope)
    if atom.isotope==None:
        rho=0.0
        for iso in atom.isotopes:
            atom=Atom(element,iso)
            rho+=vapour_number_density(Temperature,element)*atom.mass*atom.abundance
        return rho
    else:
        return vapour_number_density(Temperature,element)*atom.mass


def speed_likely(Temperature,element,isotope):
    r"""This function calculates the most likely speed (in meters per second)
    of an atom in a vapour assuming a Maxwell-Boltzmann velocity distribution.
    This is simply

    sqrt(2*k_B*T/m)

    where k_B is Boltzmann's constant, T is the temperature (in Kelvins) and
    m is the mass of the atom (in kilograms).

    >>> print speed_likely(25+273.15,"Rb",85)
    241.638108688

    >>> print speed_likely(25+273.15,"Cs",133)
    193.142579342

    """
    atom = Atom(element, isotope)
    return sqrt(2*Temperature*k_B/atom.mass)

def speed_average(Temperature,element,isotope):
    r"""This function calculates the average speed (in meters per second)
    of an atom in a vapour assuming a Maxwell-Boltzmann velocity distribution.
    This is simply

    sqrt(8*k_B*T/m/pi)

    where k_B is Boltzmann's constant, T is the temperature (in Kelvins) and
    m is the mass of the atom (in kilograms).

    >>> print speed_average(25+273.15,"Rb",85)
    272.65940782

    >>> print speed_average(25+273.15,"Cs",133)
    217.938062809

    """
    atom = Atom(element, isotope)
    return sqrt(8*k_B*Temperature/atom.mass/pi)


def speed_rms(Temperature,element,isotope):
    r"""This function calculates the average speed (in meters per second)
    of an atom in a vapour assuming a Maxwell-Boltzmann velocity distribution.
    This is simply

    sqrt(8*k_B*T/m/pi)

    where k_B is Boltzmann's constant, T is the temperature (in Kelvins) and
    m is the mass of the atom (in kilograms).

    >>> print speed_rms(25+273.15,"Rb",85)
    295.945034349

    >>> print speed_rms(25+273.15,"Cs",133)
    236.550383496

    """
    atom = Atom(element, isotope)
    return sqrt(3*Temperature*k_B/atom.mass)


def collision_rate(Temperature, element, isotope):
    r"""This function recieves the temperature of an atomic vapour (in Kelvin),
    the element, and the isotope of the atoms, and returns the angular
    frequency rate of collisions (in rad/s) in a vapour assuming a
    Maxwell-Boltzmann velocity distribution, and taking the cross section
    of the collision to be

        sigma=pi*(2*r)**2

    where r is the atomic radius. colission rate returned is

        gamma_col=2*pi* ( sigma * v * n )

    where v is the average velocity of the distribution, and n is the
    number density of the vapour.

    A few examples (in Hz):
    >>> print collision_rate(25 + 273.15, "Cs", 133)/2/pi
    9.0607260277

    For cesium collisions become important for temperatures above 120 Celsius.
    >>> print collision_rate(120 + 273.15, "Cs", 133)/2/pi
    10519.235289

    """
    atom=Atom(element,isotope)

    sigma=pi*(2*atom.radius)**2
    v=speed_average(Temperature,element,isotope)
    n=vapour_number_density(Temperature,element)
    return 2*pi*sigma*v*n


# [1] Wavelengths, Transition Probabilities, and Energy Levels for the
#     Spectra of Cesium (Cs I–Cs LV),
#     J. E. Sansonetti,
#     J. Phys. Chem. Ref. Data 38, 761–923 (2009)
#     DOI:10.1063/1.3132702
#
# [2] Measurement of the 6DJ hyperfine structure of cesium using resonant
#     two-photon sub-Doppler spectroscopy
#     Kortyna
#
# [3] Cesium D Line Data
#     Daniel Adam Steck
#
# [4] http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=html&isotype=some
#
# [5] http://www.nndc.bnl.gov/nudat2/chartNuc.jsp
#
# [6] Rubidium 85 D Line Data
#     Daniel Adam Steck
# [7]
# @article{atomic_radii-Slater,
# author = {J. C. Slater},
# title = {Atomic Radii in Crystals},
# journal = {The Journal of Chemical Physics},
# volume = {41},
# number = {10},
# pages = {3199-3204},
# year = {1964},
# doi = {10.1063/1.1725697},
#
# URL = {
#        http://dx.doi.org/10.1063/1.1725697
#
# },
# eprint = {
#        http://dx.doi.org/10.1063/1.1725697
# }
# }
#
# [8] Rubidium 87 D Line Data
#     Daniel Adam Steck


# Not reviewed yet:
#
# [21] Arimondo E, Inguscio M and Violino P 1977
# Rev. Mod. Phys.49 31
#
# [22] Gustafsson J, Rojas D and Axner O 1997
# Spectrochim. Acta B 52 1937
#
# [23] Rapol UD, Krishna A and Natarajan V 2003
# Eur. Phys. J. D 23 185
#
# [24] Banerjee A, Das D and Natarajan V 2004
# Europhys. Lett. 65 172
