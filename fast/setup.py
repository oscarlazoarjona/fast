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

import os

bashrc_file=file(os.path.expanduser('~')+'/.bashrc',"r")
bashrc_text=bashrc_file.read()
bashrc_file.close()


#We get the current directory (from which python setup.py)
install_dir= os.path.dirname(os.path.abspath("setup.py"))
install_dir=install_dir[:-5]

#We find out if config.py already has the variable fast_path setup
#and add it if not.
f=file("config.py","r")
config=f.read()
f.close()

if "fast_path=" not in config:
    config+="\n"
    config+='fast_path="'+install_dir+"/fast"+'"\n'
    f=file("config.py","w")
    f.write(config)
    f.close()

extra_line='''export PYTHONPATH="${PYTHONPATH}:'''+install_dir+'"'

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
	from fast.all import *
	print "Install finished."
	print "Try"
	print "$ python tests.py"
	print "to verify everything runs smoothly."
except:
	raise RuntimeError,'Fast could not be found, something went wrong with the installation.'
	
