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

"""The basic configuration of FAST."""

from fast import __file__
# Whether to use parallelization through OpenMP.
parallel = True
parallel = False
# Whether to use NETCDF binary files for data communication.
use_netcdf = True
use_netcdf = False

# An integer between 0 and 2 to control which tests are ran.
run_long_tests = 0

# The install directory for FAST:
fast_path = __file__[:-len("__init__.pyc")]
