Mercury6
===============================
A version to handle binary stars
---------------------------------

We've created a modified version of the original Mercury code from John Chambers (http://www.arm.ac.uk/~jec/home.html) that can work for either single stars or binary stars.  This uses the exact same input formats as the original. A full explanation of all inputs and outputs is given in ``mercury.man``. 

Notable contents of this repository
---------------------------

* ``mercury6_ras.for``: The main program file that has been modified from the original ``mercury6_2.for``.  It requires ``mercury.inc`` and ``swift.inc`` to compile.
* ``close6_ras.for``: The code to create close encounter files.  It requires ``mercury.inc`` and ``swift.inc`` to compile.
* ``element6.for``: The code to create output element files.  It requires ``mercury.inc`` and ``swift.inc`` to compile and has not been modified from the original version.
* ``mercury.inc``:  This is the file that controls the binary.  At the bottom, I have added three options for the binary

    * ``isbinary``: If you want to have a central binary, set this to ``.TRUE.``.  If you want to run Mercury like the original veriosn, set this to ``.FALSE.``


How to compile and run
----------------------
Use your favorite fortran compiler, such as ``gfortran`` or ``f77``, to create an executable.  For instance, on Linux or Mac, try::

   gfortran -o mercury6 mercury6_ras.for
   gfortran -o close6   close6_ras.for
   gfortran -o element6 element6.for

There will likely be warnings due to the code being written in Fortran77, but it should compile

Disclaimers
------------

* This is designed for a central binary.  However, it *should* work for a binary with s-type planets, although the radius calculations will have to be tinkered with (alternatively, try setting your stars as the central object and the first object in the big.in file, regardless of true order?). Use at your own risk.
* The changes have only been tested with the RADAU integrator.  Use other integrators at your own risk.
* A routine (``mco_h2jras``) uses a bubble sort algorithm.  This can slow things down if you have a lot of massive bodies.
* I've fixed all the errors I've found.  If you find a bug, let me know so we can try to fix it.
