Installation
============


Preparing your mind
-------------------

HiPart is developed and tested in modern Linux environments. The
installation and usage will therefore be relatively easy on Linux. If you want
to use HiPart on other operating systems such as Windows or OSX, you should
have a minimal computer geek status to get it working, and this document may
provide at best some clues. We are always interested in hearing from your
installation adventures.

Also note that the current distribution of HiPart is not stable. You are
supposed to simply check out the latests development version and try it out.
There may be obvious and obscure bugs, although we try hard to get the code as
reliable as possible. Note that this is Open Source software and that no warranty
of any kind is implied or expressed. If the usage of this program ruins your
academic career, you probably deserve it. Always know what you and your
computers are doing.

At some point you may think: "Hmm, this program is slow." You may be right.
The source code focusses on readability and extensibility, not on speed. It
is mostly written in Python, with several parts in C and Fortran. I mainly
care about the outcome, and less about the time it takes to obtain the outcome.


Preparing your Linux system
---------------------------

Some software packages should be installed before HiPart can be installed or
used. It is recommended to use the software package management of your Linux
distribution to install these dependencies.

The following software must be installed for HiPart:

 * Python 2.5, 2.6 or 2.7 (including the header files): http://www.python.org/doc/
 * Numpy 1.0 or later: http://numpy.scipy.org/
 * A Fortran and a C compiler supported by the F2PY module in Numpy, e.g.
   gfortran and gcc: http://gcc.gnu.org/
 * Setuptools for Python: http://pypi.python.org/pypi/setuptools
 * Git: http://git-scm.com/

On Ubuntu, the following command will take care of all the installation
work::

    sudo apt-get install python python-dev python-numpy gfortran gcc python-setuptools git-core


Installing the bleeding edge version of HiPart
----------------------------------------------


The following series of commands will download the latest versions of the
MolMod package (required) and HiPart, and will then install them into your
home directory. Make sure you execute these commands in some sort of temporary
directory. ::

    git clone git://molmod.ugent.be/git/molmod.git
    git clone git://molmod.ugent.be/git/hipart.git
    (cd molmod; ./setup.py install --home=~)
    (cd hipart; ./setup.py install --home=~)

In order to activate the Python modules and the executable scripts installed
in your home directory, the following lines need to be added to your login
scripts:

* Bash users: add the following two lines to your ``~/.bashrc`` file::
    
    export PYTHONPATH=$HOME/lib/python:$PYTHONPATH
    export PATH=$HOME/bin:$PATH

* TC Shell users: add the lines to your ``~/.tcshrc`` file::

    setenv PYTHONPATH $HOME/lib/python:$PYTHONPATH
    setenv PATH $HOME/bin:$PATH

If you don't know which shell you are using, you are probably using Bash. Note
that some of these lines may already be present. Now log out and log in again.
You are ready to start using HiPart!


Upgrading to the latest development release
-------------------------------------------

In case you want to upgrade HiPart to the latests development version after a
previous install, then execute the following commands (in the same directory)::
    
    (cd molmod; git pull; rm ~/lib/python/molmod*; ./setup.py install --home=~)
    (cd hipart; git pull; rm ~/bin/hi-*.py; rm ~/lib/python/HiPart*; ./setup.py install --home=~)


Testing your installation
-------------------------

For the development and testing one needs to install two additonal packages:

 * Nosetests: http://somethingaboutorange.com/mrl/projects/nose/0.11.2/
 * Sympy: http://www.sympy.org/

On Ubuntu, the following command will take care of the installation::

    sudo apt-get install python-nose python-sympy

Once these are installed go to the directory where the HiPart source code was
downloaded and execute the following commands::

    cd hipart
    ./setup.py nosetests

This will run a series of tests to check the validity of the outcomes generated
by HiPart. If some tests fail, post the output of the tests on the mailing list.
