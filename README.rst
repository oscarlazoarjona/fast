FAST
====
This is FAST, an acronym for FAST Atomic Spectroscopy from Theory.

This software derives and solves optical Bloch equations efficiently. This can
be done for an arbitrary number of atomic states and radiation fields. The
equations are used primarily to produce theoretical spectra.

A library of properties of alkali atoms (currently rubidium and cesium) is
included to produce physically accurate equations from first principles or
from tables of measured quantities. For instance, the saturation of the
Rubidium D2 line

.. image:: https://raw.githubusercontent.com/oscarlazoarjona/fast/master/examples/folder_03___Rb87_one_photon/suite_levels.png
 
 is calculated from first principles

.. image:: https://raw.githubusercontent.com/oscarlazoarjona/fast/master/examples/folder_03___Rb87_one_photon/suite_3power.png

Symbolic derivations of the equations are also implemented.

Installing
----------
This software requieres gfortran, numpy, matplotlib, lapack, sympy and scipy to
run. This is done in the following steps:

**1.-** Install the Python dependencies by installing Anaconda, that is Continuum's
Python distribution from https://www.continuum.io/downloads

**2.-** Install the non-Python dependencies.
In Ubuntu, the fortran dependencies (and git) can be installed with:
::

    $ sudo apt-get install gfortran liblapack-dev git

Other operative systems are not supported at the moment, but FAST should become
OS independent in the following releases (with Fortran being optional).

**3.-** Install FAST.
To get the latest stable version of FAST use the command:
::

    $ pip install fast-atomic

To get the latest unstable version (this requires git):
::

    $ pip install git+git://github.com/oscarlazoarjona/fast

To upgrade to the latest stable version:
::

    $ pip install fast-atomic --upgrade

To upgrade to the latest unstable version (this requires git):
::

    $ pip install git+git://github.com/oscarlazoarjona/fast --upgrade

To uninstall:
::

    $ pip uninstall fast-atomic

Using FAST
----------

Once this completed example jupyter notebooks can be downloaded from
https://github.com/oscarlazoarjona/fast-notebooks

And they can be ran using
::

    $ jupyter notebook

Documentation
-------------
The current documentation can be found in

https://oscarlazoarjona.github.io/fast

Enjoy!
