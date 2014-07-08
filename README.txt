README file for the Quantum Information Toolkit.
Version 0.10.0 (beta)
Released 2014-07-??


Introduction
============

Quantum Information Toolkit (QIT) is a free, open source
MATLAB 7.6 toolkit for various quantum information and computing
-related purposes, distributed under GPL.
There is also a Python version available, with equivalent functionality.
The latest version can be downloaded from the project website,

  http://sourceforge.net/projects/qit/

The toolkit is installed by simply unzipping it, or downloading it
directly from the Git repository. To initialize the toolkit in MATLAB, run
the init.m script in the toolkit root directory. To get an overview of
the features and capabilities of the toolkit, run examples/tour.m.


Octave users
============

The object-oriented features used in the various classes, as well as
the package directories require MATLAB 7.6 to work. Unfortunately they
do not work at all in Octave 3.6. Apparently there are plans to add
MATLAB classdef compatibility to Octave eventually. In the meantime,
users looking for a free implementation are encouraged to try the
Python version of QIT.


License
=======

QIT is released under the GNU General Public License version 3.
This basically means that you can freely use, share and modify it as
you wish, as long as you give proper credit to the authors and do not
change the terms of the license. See LICENSE.txt for the details.


Design notes
============

The main design goals for this toolkit are ease of use and
comprehensiveness. It is primarily meant to be used as a tool for
hypothesis testing, small simulations, and learning, not for
computationally demanding simulations. Hence optimal efficiency of the
algorithms used is not a number one priority.
However, if you think an algorithm could be improved without
compromising accuracy or maintainability, please let the authors know
or become a contributor yourself!


Bibliography
============

Some of the m-files have literature references relevant to the
algorithms or concepts used. Each reference is on its own line starting
with the characters "%!". One can compile a list of all the references
in the toolkit using the shell command "grep '%!' */*.m".
All the references are also listed in refs.bib in Bibtex format.


Contributing
============

QIT is an open source project and your contributions are welcome.
To keep the code readable and maintainable, we ask you to follow these
coding guidelines:

* Use the file template.m as a template for new functions
* Instead of using multiple similar functions, use a single function
   performing multiple related tasks, see e.g. @state/state.m
* Fully document the functions and classes (purpose, calling syntax,
   output, approximations used, assumptions made...)
* Add relevant literature references using the %! syntax
* Use the error() function on invalid input
* Use code indention to improve readability
* Use variables sparingly, give them descriptive (but short) names
* Use comments to explain the logic of your code
* When you add new functions also add testing scripts for validating
  your code. If you modify existing code, make sure you didn't break
  anything by checking that the testing scripts still run flawlessly.


Authors
=======

Ville Bergholm                2008-2014
Jacob D. Biamonte             2008-2009
James D. Whitfield            2009-2010
