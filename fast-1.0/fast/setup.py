# -*- coding: utf-8 -*-
#Oscar Gerardo Lazo Arjona

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
	
