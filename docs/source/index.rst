.. goss documentation master file, created by
   sphinx-quickstart on Fri Feb 18 11:47:26 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

goss
====


``goss`` is a C++ library for solver ordinary differential equations. ``pygoss`` is a python library containing bindings to the ``goss`` library.
The general idea is that you define your ODE in a gotran_ and hand the ode over to ``goss``.

.. _gotran: https://github.com/ComputationalPhysiology/gotran

Content
-------

.. toctree::
   :maxdepth: 1

   install
   CONTRIBUTING
   presentation
   demo

Source code
-----------
The source code can be found on `GitHub <https://github.com/ComputationalPhysiology/goss>`__

License
-------
``goss`` is licensed under the GNU LGPL, version 3 or (at your option) any later version.
``goss`` is Copyright (2011-2022) by the authors and Simula Research Laboratory.

Authors
-------
``goss`` is developed at `Simula Research Laboratory <https://www.simula.no>`__. The core developers are

- Henrik Finsberg (henriknf@simula.no)
- Cécile Daversin-Catty (cecile@simula.no)
- Johan Hake


API documentation
------------------
.. toctree::
   :maxdepth: 1

   modules


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
