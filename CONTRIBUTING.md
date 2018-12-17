Read our [documentation](https://compwa.readthedocs.io/en/latest/contribute.html) on how to contribute!
For your convenience, you can find the coding conventions below.

# Coding Conventions

## C++ Code

### LLVM Conventions
We use the [LLVM Coding Conventions](http://llvm.org/docs/CodingStandards.html)
for the ComPWA source code. Some more conventions are defined below. For some
IDE plugins for code formatting see [Plugins](https://compwa.readthedocs.io/en/latest/contribute.html#plugins).

### Additional Conventions
The following coding conventions should be used when adding code to the
framework. This increases readability and makes working with code of different
teams easier.

### Classes, Functions, Variables
* Use camel casing for class and function names (camel case: if a new word
  starts within the name, this word starts again with a capital letter)
* No underscores in class or function names!
* Names of classes should begin with a capital letter (example: MyClass)
* Names of functions should begin with lower case letter (example:
  myFunction())
* Member variables should start with a capital letter!

### Use meaningful types and Names
* Try to come up with names for classes or functions, that describe it well
* Try to make the name as short as possible, but avoid short forms like 
  ``getAmpMaxVal()``
* Try to use meaningful types! (Example: To save indices corresponding to a
  container, use: ``std::vector<unsigned int> IndexList;`` NOT: 
  ``std::vector<int> IndexList;``)

### Const correctness
Try to follow const correctness. So member functions that do not alter the
class instance state should have the const keywords at the end. And try to use
const references instead of copies (except base types) when you can. Example

``` c++

   std::vector<unsigned int> MyClass::findEvenNumbers(
       const std::vector<unsigned int>& number_list) const {
     ...  
   }
```
### Forward declarations
Try to forward declare as much as possible

### Pointers and references
Use ``int *pi; int &ri; `` instead of 
`` int* pi; int& ri;``.

### Spaces
Use space in the following manner:

``` c++

   i = x + 1;
   a = method(a, b);
   if (true) {
     //do something
   }
```
### Comparison
When comparing a variable with a constant always use the constant as left hand
side. E.g. ``float *pf; if (NULL == pf);``

## Python Code

We use pep8. Available automatic source formatters are `flake8` and `autopep8`.
