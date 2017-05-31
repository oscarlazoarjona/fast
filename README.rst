FAST
====
This is FAST, an acronym for FAST Atomic Spectroscopy from Theory.

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

Enjoy!
