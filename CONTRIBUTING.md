Read our documentation on
[how to contribute](https://compwa.github.io/contribute.html)!
For your convenience, you can find the coding conventions below.

# Coding Conventions

## LLVM Conventions

Our coding conventions for the ComPWA source code are based on the
[LLVM Coding Conventions](http://llvm.org/docs/CodingStandards.html) and are
specified in the [`.clang-format` file](./.clang-format). We highly recommend
you to use automatic source formatter plugins in your code editor or IDE
(preferably using [ClangFormat](https://clang.llvm.org/docs/ClangFormat.html)).
Some more conventions are defined below.

## Conventions Worth Mentioning

The following coding conventions cannot be enforced by `clang-format`, but
**should** also be used!

### Classes, Functions, Variables

- Use camel casing for class and function names (camel case: if a new word
  starts within the name, this word starts again with a capital letter)
- No underscores in class or function names!
- Names of classes should begin with a capital letter (example: MyClass)
- Names of functions should begin with lower case letter (example: myFunction())
- Member variables should start with a capital letter!

### Use meaningful types and Names

- Try to come up with names for classes or functions, that describe it well
- Try to make the name as short as possible, but avoid short forms like
  `getAmpMaxVal()` and rather use `getAmplitudeMaximum()`
- Try to use meaningful types! (Example: To save indices corresponding to a
  container, use: `std::vector<unsigned int> IndexList;` NOT:
  `std::vector<int> IndexList;`)

### Const correctness

Try to follow const correctness. So member functions that do not alter the class
instance state should have the const keywords at the end. Same goes for
references. Either use copies or const references instead of plain references.

Example:

```c++

   std::vector<unsigned int> MyClass::findEvenNumbers(
       const std::vector<unsigned int>& number_list) const {
     ...
   }
```

### Forward declarations

Try to forward declare as much as possible

# Git conventions

- In the master branch, it should be possible to compile and test the framework
  **in each commit**. In your own topic branches, it is recommended to commit
  frequently (WIP keyword), but
  [squash those commits](https://git-scm.com/book/en/v2/Git-Tools-Rewriting-History)
  to compilable commits upon submitting a merge request.
- Please use [conventional commit messages](https://www.conventionalcommits.org/):
  start the commit subject line with a semantic keyword (see e.g.
  [Angular](https://github.com/angular/angular/blob/master/CONTRIBUTING.md#type)
  or [these examples](https://seesparkbox.com/foundry/semantic_commit_messages),
  followed by [a column](https://git-scm.com/docs/git-interpret-trailers), then
  the message. The subject line should be in imperative mood—just imagine the
  commit to give a command to the code framework. So for instance:
  `feat: add coverage report tools` or `fix: remove`. The message should be in
  present tense, but you can add whatever you want there (like hyperlinks for
  references).
- Try to keep test coverage high. The coverage of is evaluated for each merge
  request (see an example
  [here](https://github.com/ComPWA/ComPWA/pull/288#issuecomment-581402267)).
