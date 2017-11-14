== LLVM Conventions ==
We use the [http://llvm.org/docs/CodingStandards.html LLVM Coding Conventions ] for the ComPWA source code. Some more conventions are defined below. For some IDE plugins for code formatting see [[#Plugins| Plugins]].

== Additional Conventions ==
The following coding conventions should be used when adding code to the framework. This increases readability and makes working with code of different teams easier.

=== Classes === 
*  Names of classes should begin with a capital letter
*  If a new word starts within the name, this word starts again with a capital letter
*  No underscores in class names
*  Member variables should start with an underscore

=== Forward declarations ===
Try to avoid all forward declarations in public headers! In case it is necessary comment on it.

=== Pointers and References ===
Use <code>int *pi; int &ri; </code> instead of <code> int* pi; int& ri;</code>.
=== Spaces ===
Use space in the following manner:<source>
    i = x + 1 ; <br>
    a = method(a, b);
    if (true) {
       //do something
    }
</source>
=== Comparison ===
When comparing a variable with a constant always use the constant as left hand site. E.g. <code>float *pf; if ( NULL == pf );</code>
=== Documentation ===
[http://www.doxygen.org Doxygen] ([http://www.stack.nl/~dimitri/doxygen/ manual]) is used for documentation. We use the comment style as suggested by the [http://llvm.org/docs/CodingStandards.html LLVM Coding Conventions]. See [http://www.stack.nl/~dimitri/doxygen/formulas.html here] in order to learn how to use latex equations in your comments. Further tutorials on the usage of doxygen can be found
[http://www.stack.nl/~dimitri/doxygen/docblocks.html#docexamples here] and 
[http://justcheckingonall.wordpress.com/2008/07/20/simple-doxygen-guide/ here].

== Plugins ==
=== Eclipse ===
To switch the default formatter of eclipse to a LLVM-style one, first install the [http://www.eclipse.org/mpc/ marketplace ] via the Eclipse update:  
* Help -> Install new Software -> All available sties -> type "marketplace" in the search box -> install  

Then install the CppStyle plugin with the [https://marketplace.eclipse.org/content/cppstyle#group-details marketplace]:  
* Help -> Marketplace -> type "CppStyle" in the search box -> install  

Afterwards, go to  
* Window -> Preferences -> C++ -> Code Style -> Formatter  
* Switch "Code Formatter" from "[built in]" to "CppStyle (clang-format)"  

When you let format your code by Eclipse it is now based on clang-format with the standard Google style.
=== XCode ===
Since Xcode 8.0 third party plugins are pretty much restricted. Nevertheless, you can try [https://github.com/mapbox/XcodeClangFormat XcodeClangFormat].
