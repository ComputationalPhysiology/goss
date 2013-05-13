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

//%{
//  // Type conversion function for state and values in ODE::eval
//  SWIGINTERNINLINE double* _state_array_typemap(ODE* self, PyObject* input)
//  {
//    // Check type
//    if (!PyArray_Check(input))
//      SWIG_exception(SWIG_TypeError, "Numpy array expected");
//
//  // Get PyArrayObject
//  PyArrayObject *xa = reinterpret_cast<PyArrayObject*>($input);
//
//  // Check data type
//  if (!(PyArray_ISCONTIGUOUS(xa) && PyArray_TYPE(xa) == NPY_DOUBLE))
//    SWIG_exception(SWIG_TypeError, "Contigous numpy array of doubles expected."
//           " Make sure the numpy array is contiguous, and uses dtype=np.float_.");
//
//  // Check size of passed array
//  if ( PyArray_SIZE(xa) != self->num_states() )
//    SWIG_exception(SWIG_ValueError, "Expected a numpy array of the same size "
//		   "as number of states.");
//  
//  $1 = (double *)PyArray_DATA(xa);
//
//  }
//
//%}
//
//
//%extend goss::ODE
//{
//  void eval(PyObject* states, double t, PyObject* values)
//  {
//    
//  }
//}

%extend goss::Progress {

void _add(int incr) {
    for (int j=0;j<incr; ++j)
        (*self)++;
}

void _set(double value) {
    *self = value;
}

%pythoncode
%{
def __iadd__(self, other):
    if isinstance(other, int):
        self._add(other)
    elif isinstance(other, float):
        self._set(other)
    return self

def update(self, other):
    "Update the progress with given number"
    if isinstance(other, float):
        self._set(other)
%}

}

//-----------------------------------------------------------------------------
// Use traceback in debug message
// Reimplement info
//-----------------------------------------------------------------------------
%pythoncode %{
def debug(message):
    import traceback
    file, line, func, txt = traceback.extract_stack(None, 2)[0]
    __debug(file, line, func, message)

%}
