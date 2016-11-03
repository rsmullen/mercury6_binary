Mercury6
===============================
A version to handle binary stars
---------------------------------

We've created a modified version of the original Mercury code from John Chambers (http://www.arm.ac.uk/~jec/home.html) that can work for either single stars or binary stars.  This uses the exact same input formats as the original. A full explanation of all inputs and outputs is given in ``mercury.man.`` 

If you want a nice python wrapper for initializing Mercury simulations, see Adam Sutherland's code quicksilver at https://github.com/adamsutherland/quicksilver.

Notable contents of this repository
---------------------------

*    ``mercury6_ras.for:`` The main program file that has been modified from the original ``mercury6_2.for``.  It requires ``mercury.inc`` and ``swift.inc`` to compile.  These are the notable changes made within the code.
 
     +   ``mfo_user_centralradius`` allows the user to set the prescription for the central binary radius.
     +   When the user uses a central binary, the central star is considered as a big body instead of as a central object.
     +   The calculation of the Hill radius uses the object's location instead of its semi-major axis.  Also in this routine, we use a modified Jacobi routine that resorts planets on every call. 
     +   We apply the De Souza Torres & Anderson (2008) bug fix.
     +   Other changes can be found by searching ``RAS`` in the file.

*    ``close6_ras.for:`` The code to create close encounter files.  It requires ``mercury.inc`` and ``swift.inc`` to compile.
*    ``element6.for:`` The code to create output element files.  It requires ``mercury.inc`` and ``swift.inc`` to compile and has not been modified from the original version.
*    ``mercury.inc:``  This is the file that controls the binary.  At the bottom, I have added three options for the binary.

     +   ``isbinary:`` If you want to have a central binary, set this to ``.TRUE.;`` If you want to run Mercury like the original version, set this to ``.FALSE.``.
     +   ``cenname:`` The name you want the central object (the primary) to have.
     +   ``allowclose:`` A flag to allow or forbid close encounters/collisions between the central stars.  Use ``.FALSE.`` to forbid and ``.TRUE.`` to allow.  If you have a very close binary that is stable for the length of the integration, set this to false to speed up the program.

*     ``Kepler47/:``  This directory contains example input files to run a circumbinary example, Kepler 47.  The ephemeris was taken from Kratter and Shannon (2014).  To run, compile the code with the ``mercury.inc`` contained in the folder.
*     ``SolarSystem/:`` This directory contains example input files to run a single star example (the solar system), as taken from the original Chambers tar file. To run, compile the code with the ``mercury.inc`` contained in the folder.
*     ``Original/:``  This contains the unaltered code, just in case.


How to compile and run
----------------------

Use your favorite FORTRAN compiler, such as ``gfortran`` or ``f77``, to create an executable.  For instance, on Linux or Mac, try::

   gfortran -o mercury6 mercury6_ras.for
   gfortran -o close6   close6_ras.for
   gfortran -o element6 element6.for

There will likely be warnings due to the code being written in FORTRAN77, but it should compile.  Copy or link the executable wherever you want (wherever your input files are) to run your code using ``./mercury6``.

Tricks and Caveats
------------------

Unfortunately, the code needs to be recompiled any time parameters in the ``mercury.inc`` file get changed.

The binary stars are the central body in the ``param.in`` file and the first body in ``big.in``.

The coordinates must be in central body.  We've found it most reliable to draw our planets in Jacobi coordinates and then convert into central body after.  Similarly, we've found it easiest to not rely on the built-in orbital element converter or the output Jacobi coordinate conversion routine, so we prefer to output in central body coordinates and convert after the fact. 

We don't output the change in mass of the central body.  If that information is important, it can be reconstructed by adding the mass of bodies that collided.


Disclaimers
------------

* This is designed for a central binary.  However, it *should* work for a binary with s-type planets, although the radius calculations will have to be tinkered with (alternatively, but untested, try setting your stars as the central object and the first object in the big.in file, regardless of true order?). Use at your own risk.
* The changes have only been tested with the RADAU integrator.  Use other integrators at your own risk.
* A routine (``mco_h2jras``) uses a bubble sort algorithm.  This can slow things down if you have a lot of massive bodies.
* I've fixed all the errors I've found.  If you find a bug, let me know so we can try to fix it.  
