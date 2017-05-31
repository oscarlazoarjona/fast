This is FAST, an acronym for FAST Atomic Spectroscopy from Theory.

Installation instructions:

This software requieres gfortran, numpy, matplotlib, lapack, sympy and scipy to run.

All of this can be installed with the following commands in Ubuntu:

1.- Install the fortran dependencies:
 $ sudo apt-get install gfortran liblapack-dev

2.- Install Anaconda, that is Continuum's python distribution:
 $ wget https://repo.continuum.io/archive/Anaconda2-4.3.0-Linux-x86_64.sh
 $ chmod u+x Anaconda2-4.3.0-Linux-x86_64.sh
 $ ./Anaconda2-4.3.0-Linux-x86_64.sh

The installer will ask you to accept a few things: the terms and conditions,
the default location for the install, and adding anaconda to your .bashrc
file. Agree to all of them. All python and pip commands in this readme
should be ran using the Anaconda distribution, and not the default Python
distribution that might be available in your operative system.

3.- Open a new terminal in order to reload you .bashrc file and continue
from the new terminal.

4.- Install tabulate:
 $ pip install tabulate

5.- Install FAST itself
../fast$ python setup.py

########################################################################
                              Using FAST
########################################################################
Once this completed tests can be run with
../fast$ python examples/tests.py

or the notebooks can be browsed with

 $ jupyter notebook

Enjoy!
