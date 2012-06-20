Installation
############


Disclaimer
==========

HiPart is developed and tested in modern Linux environments. The
installation instructions below are given for a Linux system only. If you want
to use HiPart on other operating systems such as Windows or OSX, you should
have a minimal computer geek status to get it working. We are always interested
in hearing from your installation adventures.

At some point you may think: "Hmm, this program is slow." You may be right.
The source code focuses on readability and extensibility, not on speed. It
is mostly written in Python, with several parts in C and Fortran. I mainly
care about the outcome, and less about the time it takes to obtain the outcome.


MolMod dependency
=================

`MolMod <http://molmod.github.com/molmod/>`_ is a Python library used by most
Python programs developed at the CMM. It must be installed before HiPart can
be used or installed. Installation and download instructions can be found in the
`molmod documentation <http://molmod.github.com/molmod/tutorial/install.html>`_.
The instructions below only work if the MolMod package is installed.


External dependencies
=====================

Some software packages should be installed before HiPart can be installed or
used. It is recommended to use the software package management of your Linux
distribution to install these dependencies.

The following software must be installed for HiPart:

* Python 2.5, 2.6 or 2.7 (including the header files): http://www.python.org/doc/
* Numpy 1.0 or later: http://numpy.scipy.org/
* A Fortran and a C compiler supported by the F2PY module in Numpy, e.g.
  gfortran and gcc: http://gcc.gnu.org/
* Git: http://git-scm.com/

In the tutorial we also use the following:

* wget: http://www.gnu.org/software/wget/

Most Linux distributions can install this software with just a single command:

* Ubuntu 12.4::

    sudo apt-get install python python-dev python-numpy gfortran gcc git-core wget

* Debian 5. One first has to become root because the sudo program is not
  configured by default. ::

    su -
    apt-get install python python-dev python-numpy gfortran gcc git-core wget
    exit

* Fedora 17. ::

    sudo yum install python-devel numpy numpy-f2py gcc-gfortran gcc git wget

* Suse 11.2::

    sudo zypper install python-devel python-numpy gcc gcc-fortran git wget

  There seems to be something odd going on with the default Python configuration
  on Suse installations. You have to edit the file
  ``/usr/lib64/python2.4/distutils/distutils.cfg`` or
  ``/usr/lib32/python2.4/distutils/distutils.cfg``, depending on the CPU
  architecture, to comment out the line ``prefix=/usr/local`` with a ``#``
  symbol. Otherwise it is impossible to install Python packages in the home
  directory, as we will do below.


Installing the latest version of HiPart
===========================================

The following series of commands will download the latest version of HiPart,
and will then install it into your home directory. ::

    cd ~/build/
    git clone git://github.com/molmod/hipart.git
    (cd hipart; ./setup.py install --home=~)

You are now ready to start using HiPart!


A few quick checks
==================

It may be interesting to double check your installation before proceeding,
unless you `feel lucky`. The HiPart files are installed in the
following directories:

* Scripts: ``~/bin``
* Modules: ``~/lib/python`` or ``~/lib64/python``

There should be at least some files present in these directories.

All HiPart scripts start with ``hi-``, so when you enter ``hi-`` on the command
line and then hit the ``TAB`` key twice, you should see the following programs
in your terminal::

    hi-atomdb.py                hi-multipoles.py
    hi-bond-orders.py           hi-net-overlap.py
    hi-charges.py               hi-overlap-matrices-orb.py
    hi-dipoles.py               hi-overlap-matrices.py
    hi-esp-test.py              hi-spin-charges.py

The Python modules should be accessible from any Python session. This can be
checked by starting Python interactively and loading the modules manually. There
should be no errors when importing the modules::

    $ python
    Python 2.6.5 (r265:79063, Apr 16 2010, 13:57:41)
    [GCC 4.4.3] on linux2
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import hipart
    >>> import molmod
    >>> quit()


Upgrading to the latest version of HiPart and MolMod
====================================================

In case you want to upgrade HiPart to the latests development version after a
previous install, then execute the following commands (in the same directory)::

    cd ~/build/
    (cd molmod; git pull; rm -r ~/lib*/python/molmod*; ./setup.py install --home=~)
    (cd hipart; git pull; rm ~/bin/hi-*.py; rm -r ~/lib*/python/hipart*; ./setup.py install --home=~)


Testing your installation
=========================

For the development and testing one needs to install these additional packages:

* Nosetests >= 0.11: http://somethingaboutorange.com/mrl/projects/nose/0.11.2/
* Sympy >= 0.7: http://www.sympy.org/
* Sphinx >= 1.0: http://sphinx.pocoo.org/
* Scipy: http://www.scipy.org/

New Linux distributions can install this software with just a single terminal command:

* Ubuntu 12.4::

    sudo apt-get install python-nose python-sphinx python-scipy python-sympy

* Debian 5 does not have Python 2.6. Hipart does work on Debian 5, but some of
  the development tools will not work and some tests do not run. ::

    su -
    apt-get install python-nose python-sphinx python-scipy python-sympy
    exit

* Fedora 17::

    sudo yum install python-nose sphinx scipy sympy

* Suse 11.2. One needs to add a repository, but a recent Sympy is already present::

    sudo zypper ar http://download.opensuse.org/repositories/devel:/languages:/python/openSUSE_11.2/devel:languages:python.repo
    sudo zypper install python-sympy python-scipy python-nose python-sphinx

Sympy-0.6.7 can be installed as follows if your Linux distribution does not have recent version::

    wget 'http://sympy.googlecode.com/files/sympy-0.6.7.tar.gz'
    tar -xzf sympy-0.6.7.tar.gz
    cd sympy-0.6.7
    ./setup.py install --home=~

Once these dependecies are installed, go to the directory where the HiPart
source code was downloaded and execute the following commands::

    cd ~/build/hipart
    ./cleanfiles.sh
    ./setup.py build_ext -i
    nosetests -v

This will run a series of tests to check the validity of the outcomes generated
by HiPart. If some tests fail, post the output of the tests on the `mailing list
<http://molmod.ugent.be/code/wiki/HiPart/MailingList>`_.
