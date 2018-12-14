How to contribute
=================

How to use git
--------------
If you are new to git, maybe you should read some documentation first. E.g. the
`Manual <https://git-scm.com/docs/user-manual.html>`_, a
`Tutorial <http://rogerdudler.github.com/git-guide/>`_, a 
`CheatSheet <https://services.github.com/on-demand/downloads/github-git-cheat-sheet.pdf>`_.

For your convenience here is the workflow you should use if you want to
contribute:

* Log into GitHub with your account and fork ComPWA
* Get copy of repository: ``git clone git@github.com:YOURACCOUNT/ComPWA.git``
* Add the main repository as a remote: 
  ``git remote add upstream git@github.com:ComPWA/ComPWA.git``

.. note::
   You can name the repository with any name you wish. ``upstream`` is a common
   label for the main repository.
   
   Note the remote from which you cloned the repository initially is by default
   named ``origin`` (here: your fork). A local ``master`` branch is
   automatically checked out from the origin after the clone as well. You can
   list all branches with ``git branch -a``.

The following steps you repeat until your contribution is finished, and can be
added to the main repository.

.. tip::
   Remember to commit frequently to not have single commits with too many
   changes! Rebasing can help you later on to group and alter commit messages,
   so don't worry.

* ... make changes ...
* Check changes: ``git status`` and/or ``git diff``
* Stage updated files for commit:  ``git add -u``
  or add new files ``git add <list of files>``
* Commit changes: ``git commit`` (opens up editor for commit message)
* Enter a meaningful commit message. First line is a overall summary.
  Then a more detailed description if neccessary.
* Synchronize with the changes from the main repository/upstream
  
  * Fetch new changes ``git fetch upstream``
  * Reapply your current branch commits to the head of upstream master branch
    ``git rebase -i upstream/master``
  * At this point conflicts of your and the changes from the main repository
    may occur. If no conflicts appeared then you are finished and can continue
    coding or push your work onto you fork.
    Otherwise repeat these steps until finshed (you can abort the whole rebase
    process via ``git rebase --abort``):
    
    * Review the conflicts (:ref:`VS Code <contribute-vscode>` is great for
      this)
    * Mark them as resolved ``git add <filename>``
    * Continue the rebase ``git rebase --continue``
* Push your changes to your fork: ``git push origin <branchname>``
  This step synchronizes your local and fork branch, but is not required after
  every commit. But it is certainly neccessary once you are ready to merge your
  code into ``upstream``.

  .. tip::
     It can be useful to push your local branch to your fork under a different
     name using ``git push origin <local-branchname>:<remote-branchname>``

Once you think your contribution is finished and can be merged into the main
repository

* Make sure your changes are rebased and pushed into your fork
* Log into GitHub with your account and create a pull request (merge your fork
  branch into the main repository master branch)
* While open, commits pushed to your fork branch will automatically update the
  pull request

.. _contribute-report-issues:

Reporting Issues
----------------
Use the `ComPWA github issues page <https://github.com/ComPWA/ComPWA/issues>`_
to

* report problems/issues 
* file a feature request
* request modifications to existing "unpleasant" code

Please don't hesitate to report any issues, but try make sure not to post
duplicates

We are also very glad if you want to take it into your own hands and contribute
to ComPWA! 

Continuous Integration
----------------------
The master branch is automatically build using TravisCI. Probably it is 
interesting to check out the `log file <https://travis-ci.org/ComPWA/ComPWA>`_
and the projects TravisCI configuration file 
`travisCI.yml <https://github.com/ComPWA/ComPWA/blob/master/.travis.yml>`_.



Code Quality & Conventions
--------------------------

A highly recommended read for learning how to write good code:
**Clean Code, by Robert C. Martin**

Try and follow his advice, and keep in mind the boy scout rule::

  "Leave behind the code cleaner, then you found it"

C++ Code
^^^^^^^^

Specifically for C++ the `C++ Core Guidelines <https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines>`_
are extremly valuable. So please take them to heart!

LLVM Conventions 
""""""""""""""""
We use the `LLVM Coding Conventions <http://llvm.org/docs/CodingStandards.html>`_
for the ComPWA source code. Some more conventions are defined below. For some
IDE plugins for code formatting see :ref:`Plugins <plugins-label>`.

Additional Conventions
""""""""""""""""""""""
The following coding conventions should be used when adding code to the
framework. This increases readability and makes working with code of different
teams easier.

Classes, Functions, Variables
"""""""""""""""""""""""""""""
* Use camel casing for class and function names (camel case: if a new word
  starts within the name, this word starts again with a capital letter)
* No underscores in class or function names!
* Names of classes should begin with a capital letter (example: MyClass)
* Names of functions should begin with lower case letter (example:
  myFunction())
* Member variables should start with a capital letter!

Use meaningful types and Names
""""""""""""""""""""""""""""""
* Try to come up with names for classes or functions, that describe it well
* Try to make the name as short as possible, but avoid short forms like 
  ``getAmpMaxVal()``
* Try to use meaningful types! (Example: To save indices corresponding to a
  container, use: ``std::vector<unsigned int> IndexList;`` NOT: 
  ``std::vector<int> IndexList;``)

Const correctness
"""""""""""""""""
Try to follow const correctness. So member functions that do not alter the
class instance state should have the const keywords at the end. And try to use
const references instead of copies (except base types) when you can. Example

.. code-block:: c++

   std::vector<unsigned int> MyClass::findEvenNumbers(
       const std::vector<unsigned int>& number_list) const {
     ...  
   }

Forward declarations
""""""""""""""""""""
Try to forward declare as much as possible

Pointers and references
"""""""""""""""""""""""
Use ``int *pi; int &ri; `` instead of 
`` int* pi; int& ri;``.

Spaces
""""""
Use space in the following manner:

.. code-block:: c++

   i = x + 1;
   a = method(a, b);
   if (true) {
     //do something
   }

Comparison
""""""""""
When comparing a variable with a constant always use the constant as left hand
side. E.g. ``float *pf; if (NULL == pf);``

Python Code
^^^^^^^^^^^

We use pep8. Available automatic source formatters are `flake8` and `autopep8`.

Documentation
-------------

Generally try to code in such a way that it is self explanatory and its
documentation is not neccessary. Of course this ideal case is not achieved in
reality, but avoid useless comments such as ``getValue() # get's value``. Also
try to comment only parts, which really need an explanation. Because keeping 
the documentation in sync with the code is crucial, and is a lot of work.

The documentation is built with sphinx using the "read the docs" theme. For the
python code/modules ``sphinx-apidoc`` is used. The comment style is following
the pep8 conventions.

`Doxygen <http://www.doxygen.org>`_ (`manual <http://www.stack.nl/~dimitri/doxygen/>`_)
is used for documentation of the c++ code/modules. We run doxygen to produce
xml output which is translated to rst files via breathe.
We use the comment style as suggested by the
`LLVM Coding Conventions <http://llvm.org/docs/CodingStandards.html>`_.

See `here <http://www.stack.nl/~dimitri/doxygen/formulas.html>`_ in order to learn
how to use latex equations in your comments. Further tutorials on the usage of
doxygen can be found
`here <http://www.stack.nl/~dimitri/doxygen/docblocks.html#docexamples>`_ and 
`here <http://justcheckingonall.wordpress.com/2008/07/20/simple-doxygen-guide>`_.


.. _plugins-label:

Plugins
-------

Eclipse
^^^^^^^
To switch the default formatter of eclipse to a LLVM-style one, first install
the `marketplace <http://www.eclipse.org/mpc/>`_ via the Eclipse update::

   -> Help
   -> Install new Software
   -> All available sties
   -> type "marketplace" in the search box
   -> install  

Then install the CppStyle plugin with the 
`marketplace <https://marketplace.eclipse.org/content/cppstyle#group-details>`_::

   -> Help
   -> Marketplace
   -> type "CppStyle" in the search box
   -> install  

Afterwards, go to::

   -> Window
   -> Preferences
   -> C++
   -> Code Style
   -> Formatter
   -> Switch "Code Formatter" from "[built in]" to "CppStyle (clang-format)"  

When you let format your code by Eclipse it is now based on clang-format with
the standard Google style.

XCode
^^^^^
Since Xcode 8.0 third party plugins are pretty much restricted. Nevertheless,
you can try `XcodeClangFormat <https://github.com/mapbox/XcodeClangFormat>`_.

.. _contribute-vscode:

Visual Studio Code (VS Code)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
`VS Code <https://code.visualstudio.com/>`_ is a performant, feature rich, and
beautiful source code editor. Although for c++ is not as feature rich as
Eclipse or XCode, it is still quite useful (i.e. more lightweight and good git
integrations, ...). However for python development it is highly recommendable!

It can also be extended with various 
`extensions <https://code.visualstudio.com/docs/editor/extension-gallery>`_. To
bring up the extensions view, either press ``Ctrl+Shift+X`` or click on the 
Extensions icon on the lefthand sidebar.

Recommendable plugins are ``C/C++``, ``CMake``, ``Python``, 
``reStructuredText``, ...