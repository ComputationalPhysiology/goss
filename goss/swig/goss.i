/* -*- C -*- */
// Copyright (C) 2012 Johan Hake
//
// This file is part of GOSS.
//
// GOSS is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GOSS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GOSS. If not, see <http://www.gnu.org/licenses/>.

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
