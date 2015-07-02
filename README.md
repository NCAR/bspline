This is the source distribution for the NCAR/EOL BSpline C++ library. 

## Background

The BSpline package provides an implementation of a Cubic B-Spline method
devised by Vic Ooyama, and brought to our attention by James
Franklin. Franklin employed the B-Spline for general purpose filtering in
his dropsonde quality control processing program, known as `editsonde`.

NCAR used Franklin's code, and a paper published by Ooyama, to build this
generic B-Spline class library.

## Legal

See the [COPYRIGHT](./COPYRIGHT) file in the source distribution.

## Requirements

BSpline builds on Windows, Linux and Mac. Here are the build environments
for each operating system:

 * Windows: Microsoft Visual Studio 2005
 * Linux: scons and g++
 * Mac: scons and g++ (both provided by Xcode)

BSpline uses the [Doxygen](www.doxygen.org) documentation system. Doxygen
is not required to build BSpline, but if available it can be used to create
html documentation.

## Organization

 * BSpline/
   * Contains the source code to create a bspline library, and a header file to access the BSpline from your code.
 * Tests/
   * Contains several subdirectories with code for evaluating
      the code. The C++ directory is the only one that is currently compiled
      and known to work. The R directory contains code for plotting
      results using the R statistical langauge. The Fortran directory contains
      code for some unknown application.
 * Design/
   * Contains notes and and some original Ooyama code.

## Building

#### Linux/Mac

[Scons] (http://www.scons.org) is used to build BSpline.

BSpline is built by running scons in the top directory:

```
scons
```

This will create `BSpline/libbspline.a` and `Tests/C++/bspline`. The library 
will contain float and double instantiations of the BSpline
templates.

#### Windows

A Visual Studio solution (`bspline.sln`) is provided for building the
same products.

## Documentation

The BSpline code is documented using the Doxygen embedded comments.  Run
`doxygen` in the top directory, where it will find the Doxygen
configuration file:

```
doxygen
```
   
This will produce html formatted documentation in the `doc/` directory.
Open `doc/index.html` with your favorite browser.

The generated documentation is also committed to the GitHub Pages for this
project, so it can be viewed online:

 * http://ncareol.github.io/bspline/

Here are the steps to update the documentation in GitHub pages from the
source tree:

In a clone of the bspline project, clone it again into a `doc` directory:

```
cd bspline
git clone git@github.com:ncareol/bspline.git doc
```

Checkout the gh-pages branch in the `doc` subdirectory.

```
cd doc
git checkout origin/gh-pages -b gh-pages
```

Generate the doxygen output into the doc directory:

```
cd ..
scons doc
```

Commit the changes to the gh-pages branch:

```
cd doc
git commit -a -m'generated updated documentation'
```
