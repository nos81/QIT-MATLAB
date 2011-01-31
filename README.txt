README file for the Quantum Information Toolkit.
Version 0.9.9 (beta)
Released 2011-01-31



=== Introduction ===

Quantum Information Toolkit (QIT) is collection of free, open source
MATLAB 7.6 scripts, functions and classes for various quantum
information and computing -related purposes, distributed under GPL.
The latest version can be downloaded from the project website,

  http://sourceforge.net/projects/qit/

The toolkit is installed by simply unzipping it, or downloading it
directly from the SVN server. To initialize the toolkit in MATLAB, run
the init.m script in the toolkit root directory. To get an overview of
the features and capabilities of the toolkit, run examples/tour.m.



=== Octave 3.0 users ===

The object-oriented features used in the @state class, as well as the
namespace directories require MATLAB 7.6 to work. Unfortunately they
do not work at all in Octave 3.0. (It would be possible to fix this
by re-implementing @state and namespaces using traditional syntax, but
this would subtract from the usability.) However, the toolkit has some
workarounds which make it at least partially usable in Octave.

After installing the toolkit, run the shell script ./octavize.
This will perform the necessary modifications. You can restore the
toolkit to its original state by running the shell script ./restore.
To initialize the toolkit, run octave_init.m instead of init.m.



=== License ===

QIT is released under the GNU General Public License version 3.
This basically means that you can freely use, share and modify it as
you wish, as long as you give proper credit to the authors and do not
change the terms of the license. See LICENSE.txt for the details.



=== Design notes ===

The main design goals for this toolkit are ease of use and
comprehensiveness. It is primarily meant to be used as a tool for
hypothesis testing, small simulations, and learning, not for
computationally demanding simulations. Hence optimal efficiency of the
algorithms used is not a number one priority.
However, if you think an algorithm could be improved without
compromising accuracy or maintainability, please let the authors know
or become a contributor yourself!



=== Bibliography ===

Some of the m-files have literature references relevant to the
algorithms or concepts used. Each reference is on its own line starting
with the characters "%!". One can compile a list of all the references
in the toolkit using the shell command "grep '%!' */*.m".



=== Contributing ===

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



=== Authors ===

Ville Bergholm                2008-2011
Jacob D. Biamonte             2008-2009
James D. Whitfield	      2009-2010
