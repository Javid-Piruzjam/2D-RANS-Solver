This is a generic makefile for C, C++ and Fortran projects.

Features:
* Generic: No adoption necessary for changing source files
* Configurable location of src and object files
* Automatic intelligent dependency generation (C, C++)
* Separation of build configuration and Makefile
* Support for multiple build configurations
* Multiple builds in the same directory

Problems:
There are certain scenarios with renaming files which can cause errors.
Still a fix is not trivial and I wanted to keep things simple. After
renaming of files a make clean && make  fixes the problem.

Usage:
The actual build configurations are in make include files named include_<TAG>.mk.
The tag to build is configured in the Makefile. Here also the location of the
target name, the build dir where all files are generated to and the src dir can
be configured. 

Standard configuration is silent execution of the commands. This can changed by
adopting the variable Q or calling make with make Q= .
This makefile can be also adopted to build C++ or Fortran programs.

More information on the science :-) of automatic dependency file generation can be found  here:
http://www.cs.berkeley.edu/~smcpeak/autodepend/autodepend.html
http://make.paulandlesley.org/autodep.html


