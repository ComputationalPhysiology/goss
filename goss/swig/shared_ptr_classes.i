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
//
// First added:  2012-10-11
// Last changed: 2012-10-11

//=============================================================================
// SWIG directives for the shared_ptr stored classes in PyGOSS
//
// Objects of these classes can then be passed to c++ functions
// demanding a boost::shared_ptr<type>
//=============================================================================

//-----------------------------------------------------------------------------
// Un-comment these lines to use std::tr1, only works with swig version >=1.3.37
//-----------------------------------------------------------------------------
//#define SWIG_SHARED_PTR_NAMESPACE std
//#define SWIG_SHARED_PTR_SUBNAMESPACE tr1

//-----------------------------------------------------------------------------
// Include macros for shared_ptr support
//-----------------------------------------------------------------------------
%include <boost_shared_ptr.i>

//-----------------------------------------------------------------------------
// define to make SWIG_SHARED_PTR_QNAMESPACE available in inlined C++ code
//-----------------------------------------------------------------------------
%{
#define SWIG_SHARED_PTR_QNAMESPACE boost
%}

