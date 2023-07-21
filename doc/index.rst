.. gillespy documentation master file, created by
   sphinx-quickstart on Mon Sep 27 16:17:39 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

gillespy: stochastic simulation algorithms
==========================================

The `stochastic simulation algorithm (SSA) <https://en.wikipedia.org/wiki/Gillespie_algorithm>` is a computer-oriented procedure utilizing the pseudorandomness to numerically simulate Markov processes. It's a common and popular tool for investigating reactive chemical mixtures such as the biological cell, but can be quite useful in other settings where Markovian and similar processes are the subject of investigation. This package provides pure-Python, object-oriented implementations of existing SSAs, and a convenient JSON interchange format for specifying the process that is to be simulated.


Installation
------------

This package is available from pyPI::

    python3 -m pip install gillespy

or from Anaconda::

    conda install gillespy

See the toctree at the bottom of this page for links to reference documentation.


Sources and contributing
------------------------

Package soruces are `available on GitHub <https://www.github.com/wefatherley/gillespy>`, which includes a `CONTRIBUTING` file specifying ways to contribute.


Reference and citation
----------------------

When utilizing this package, please consider citing 


Copyright and license (MIT)
---------------------------

Copyright (c) 2021 Will Fatherley

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


.. toctree::
   :maxdepth: 2
   
   ssa_details
   usage_and_cookbook
   json_specification
   api_reference