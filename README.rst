Mercury6
===============================
A version to handle binary stars
---------------------------------

Contents of this repository
---------------------------
mercury


How to compile and run
----------------------
Use your favorite fortran compiler to create an executable.  For instance, on Linux or Mac, try

.. code:: python

   import rebound
   sim = rebound.Simulation()
   sim.add(m=1.0)
   sim.add(m=1.0e-3, a=1.0)
   sim.integrate(1000.)
   sim.status()

.. code::
   gfortran -o mercury6 mercury6_ras.for
   gfortran -o close6 close6_ras.for
   gfortran -o element6 element6.for


Disclaimers
------------
*This is designed for a central binary.  However, it *should* work for a 
binary with s-type planets, although the radius calculations will have to
be tinkered with (alternatively, try setting your stars as the central 
object and the first object in the big.in file, regardless of true 
order?). Use at your own risk.
*The changes have only been tested with the RADAU integrator.  Use other 
integrators at your own risk.
*A routine (mco_h2jras) uses a bubble sort algorithm.  This can slow 
things down if you have a lot of massive bodies.
*I've fixed all the errors I've found.  If you find a bug, let me know
so we can try to fix it.
