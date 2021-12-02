# Yeast_CPP
The purpose of this readme is to notify the user about the usage and dependencies of the yeast model.

## Necessary external librares
The program uses Boost, Eigen and MathGL. Most of these are too large for me to upload here along with the source code of the yeast model, so I will provide the links for download here:

- Eigen 3.3.9: https://eigen.tuxfamily.org/index.php?title=Main_Page

- Boost 1.75.0: https://www.boost.org/users/download/

- MathGl 2.4.4: https://sourceforge.net/projects/mathgl/

## Compiler setup
Before the model can be run, these libraries have to be present on the computer and they must also be linked in the compiler of choice.
During development I have been using CodeBlocks, in which linking libraries is fairly easy. Please note the following compiler settings as well:

- GNU C++17 C++ language standard is necessary

- O3 optimization is recommended for speed

In CodeBlocks, all of these can be set up in the Project -> Build options (project level setup) or the Settings -> Compiler (setup for the whole of CodeBlocks) menus.
Linking of external libraries can be done here as well.
