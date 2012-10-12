// The PyGOSS extension module
%module(package="goss.cpp", directors="1") cpp

%{
#include "goss/goss.h"

// NumPy includes
#define PY_ARRAY_UNIQUE_SYMBOL PyGOSS
#include <numpy/arrayobject.h>

%}

%init%{
import_array();
%}

%include "goss/swig/shared_ptr_classes.i"
%include "goss/swig/typemaps.i"

// Global exceptions
%include <exception.i>
%include "goss/swig/exceptions.i"

// Do not expand default arguments in C++ by generating an extra
// function in the SWIG layer. This reduces code bloat.
%feature("compactdefaultargs");

// STL SWIG string class
%include <std_string.i>

// Include information about swig version
%include "goss/swig/version.i"

// Enable automatic docstring generation
// FIXME: Consider generate from C++
%feature("autodoc", "1");

// Include the interface with pre and post modifications
%include "goss/swig/pre.i"
%include "goss/swig/goss_interface.i"
%include "goss/swig/post.i"
