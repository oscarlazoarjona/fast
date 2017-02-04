# -*- coding: utf-8 -*-

#************************************************************************
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
#************************************************************************

########################################################################
#We build the config.py file

#We get the current directory (from which python setup.py is running).
import os
install_dir=os.path.dirname(os.path.abspath("setup.py"))
install_dir=install_dir[:-5]

config_text="""#Whether to use parallelization through OpenMP.
parallel=True; parallel=False
#Whether to use NETCDF binary files for data communication.
use_netcdf=True; use_netcdf=False

#An integer between 0 and 2 to control which tests are ran.
run_long_tests=0

#The install directory for FAST:
"""
config_text+='fast_path="'+install_dir+"/fast"+'"\n'

config=file("config.py","w")
config.write(config_text)
config.close()


########################################################################
#We add an extra line to the .bashrc file.

extra_line='''export PYTHONPATH="${PYTHONPATH}:'''+install_dir+'"'


bashrc_file=file(os.path.expanduser('~')+'/.bashrc',"r")
bashrc_text=bashrc_file.read()
bashrc_file.close()

if extra_line not in bashrc_text:
	bashrc_text+="\n"+extra_line+"\n"
	
	bashrc_file=file(os.path.expanduser('~')+'/.bashrc',"w")
	bashrc_file.write(bashrc_text)
	bashrc_file.close()
	
	print "added line to ~/.bashrc:"
	print extra_line

	print "reloading ~/.bashrc"
	os.system("bash")

try:
	import fast
	print "Install finished."
	print "Try"
	print "$ python tests.py"
	print "to verify everything runs smoothly."
except:
	raise RuntimeError,'Fast could not be found, something went wrong with the installation.'
	
