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

%typemap(in) (double* component_field_states, goss::uint num_components)
{

  // Check type
  if (!PyArray_Check($input))
    SWIG_exception(SWIG_TypeError, "Numpy array expected");

  // Get PyArrayObject
  PyArrayObject *xa = reinterpret_cast<PyArrayObject*>($input);

  // Check data type
  if (!(PyArray_ISCONTIGUOUS(xa) && PyArray_TYPE(xa) == NPY_DOUBLE))
    SWIG_exception(SWIG_TypeError, "Contigous numpy array of doubles expected."
           " Make sure the numpy array is contiguous, and uses dtype=np.float_.");

  // Check size of passed array
  goss::uint size = PyArray_SIZE(xa);
  $2 = size / arg1->num_nodes();
  
  if (size % arg1->num_nodes() != 0)
    SWIG_exception(SWIG_ValueError, "Expected a numpy array of the same size "
		   "as an integer number times the number of nodes.");
  
  if ($2>=arg1->ode()->num_field_states())
    SWIG_exception(SWIG_ValueError, "Expected a numpy array of the same size "
		   "or less than the number of field states times the number "
                   "of nodes..");

  $1 = (double *)PyArray_DATA(xa);
}

%typemap(in) (const double* component_field_states, goss::uint num_components)
{

  // Check type
  if (!PyArray_Check($input))
    SWIG_exception(SWIG_TypeError, "Numpy array expected");

  // Get PyArrayObject
  PyArrayObject *xa = reinterpret_cast<PyArrayObject*>($input);

  // Check data type
  if (!(PyArray_ISCONTIGUOUS(xa) && PyArray_TYPE(xa) == NPY_DOUBLE))
    SWIG_exception(SWIG_TypeError, "Contigous numpy array of doubles expected."
           " Make sure the numpy array is contiguous, and uses dtype=np.float_.");

  // Check size of passed array
  goss::uint size = PyArray_SIZE(xa);
  $2 = size / arg1->num_nodes();
  
  if (size % arg1->num_nodes() != 0)
    SWIG_exception(SWIG_ValueError, "Expected a numpy array of the same size "
		   "as an integer number times the number of nodes.");
  
  if ($2>=arg1->ode()->num_field_states())
    SWIG_exception(SWIG_ValueError, "Expected a numpy array of the same size "
		   "or less than the number of field states times the number "
                   "of nodes..");

  $1 = (double *)PyArray_DATA(xa);
}

%typemap(in) (const goss::uint* components)
{

  // Check type
  if (!PyArray_Check($input))
    SWIG_exception(SWIG_TypeError, "Numpy array expected");

  // Get PyArrayObject
  PyArrayObject *xa = reinterpret_cast<PyArrayObject*>($input);

  // Check data type
  if (!(PyArray_ISCONTIGUOUS(xa) && PyArray_TYPE(xa) == NPY_UINT32))
    SWIG_exception(SWIG_TypeError, "Contigous numpy array of uint32 expected."
           " Make sure the numpy array is contiguous, and uses dtype=np.uint32.");

  // Get size
  goss::uint size = PyArray_SIZE(xa);

  // Check size of passed array
  if (size <= arg1->ode()->num_field_states())
    SWIG_exception(SWIG_ValueError, "Expected a numpy array of the same size "
		   "or less than the number of field states.");
  
  // Check size compared to the previously calculated size from the double* array
  if (size!=arg3)
    SWIG_exception(SWIG_ValueError, "Expected a numpy array of the same size "
		   "as len(component_field_states)/num_field_states.");

  $1 = (goss::uint *)PyArray_DATA(xa);
  
  for (goss::uint i = 0; i<size; i++)
    if ($1[i]>=arg1->ode()->num_field_states())
      SWIG_exception(SWIG_IndexError, "Expected the indices in the component array to "\
                     "be less than the number of field states.");
}

%typemap(in) double* system_field_states
{
  // Check type
  if (!PyArray_Check($input))
    SWIG_exception(SWIG_TypeError, "Numpy array expected");

  // Get PyArrayObject
  PyArrayObject *xa = reinterpret_cast<PyArrayObject*>($input);

  // Check data type
  if (!(PyArray_ISCONTIGUOUS(xa) && PyArray_TYPE(xa) == NPY_DOUBLE))
    SWIG_exception(SWIG_TypeError, "Contigous numpy array of doubles expected."
           " Make sure the numpy array is contiguous, and uses dtype=np.float_.");

  // Check size of passed array
  if ( PyArray_SIZE(xa) != arg1->num_nodes()*arg1->ode()->num_field_states() )
    SWIG_exception(SWIG_ValueError, "Expected a numpy array of the same size "
		   "as number of field states times the number of nodes.");
  
  $1 = (double *)PyArray_DATA(xa);
}

%typemap(in) const double* system_field_states
{
  // Check type
  if (!PyArray_Check($input))
    SWIG_exception(SWIG_TypeError, "Numpy array expected");

  // Get PyArrayObject
  PyArrayObject *xa = reinterpret_cast<PyArrayObject*>($input);

  // Check data type
  if (!(PyArray_ISCONTIGUOUS(xa) && PyArray_TYPE(xa) == NPY_DOUBLE))
    SWIG_exception(SWIG_TypeError, "Contigous numpy array of doubles expected."
           " Make sure the numpy array is contiguous, and uses dtype=np.float_.");

  // Check size of passed array
  if ( PyArray_SIZE(xa) != arg1->num_nodes()*arg1->ode()->num_field_states() )
    SWIG_exception(SWIG_ValueError, "Expected a numpy array of the same size "
		   "as number of field states times the number of nodes.");
  
  $1 = (double *)PyArray_DATA(xa);
}

%typemap(in) const double* system_field_params
{
  // Check type
  if (!PyArray_Check($input))
    SWIG_exception(SWIG_TypeError, "Numpy array expected");

  // Get PyArrayObject
  PyArrayObject *xa = reinterpret_cast<PyArrayObject*>($input);

  // Check data type
  if (!(PyArray_ISCONTIGUOUS(xa) && PyArray_TYPE(xa) == NPY_DOUBLE))
    SWIG_exception(SWIG_TypeError, "Contigous numpy array of doubles expected."
           " Make sure the numpy array is contiguous, and uses dtype=np.float_.");

  // Check size of passed array
  if ( PyArray_SIZE(xa) != arg1->num_nodes()*arg1->ode()->num_field_parameters())
    SWIG_exception(SWIG_ValueError, "Expected a numpy array of the same size "
		   "as number of field parameters times the number of nodes.");
  
  $1 = (double *)PyArray_DATA(xa);
}

%typemap(in) double* values
{
  // Check type
  if (!PyArray_Check($input))
    SWIG_exception(SWIG_TypeError, "Numpy array expected");

  // Get PyArrayObject
  PyArrayObject *xa = reinterpret_cast<PyArrayObject*>($input);

  // Check data type
  if (!(PyArray_ISCONTIGUOUS(xa) && PyArray_TYPE(xa) == NPY_DOUBLE))
    SWIG_exception(SWIG_TypeError, "Contigous numpy array of doubles expected."
           " Make sure the numpy array is contiguous, and uses dtype=np.float_.");

  // Check size of passed array
  if ( PyArray_SIZE(xa) != arg1->num_states() )
    SWIG_exception(SWIG_ValueError, "Expected a numpy array of the same size "
		   "as number of states.");
  
  $1 = (double *)PyArray_DATA(xa);
}

%typemap(in) const double* states
{
  // Check type
  if (!PyArray_Check($input))
    SWIG_exception(SWIG_TypeError, "Numpy array expected");

  // Get PyArrayObject
  PyArrayObject *xa = reinterpret_cast<PyArrayObject*>($input);

  // Check data type
  if (!(PyArray_ISCONTIGUOUS(xa) && PyArray_TYPE(xa) == NPY_DOUBLE))
    SWIG_exception(SWIG_TypeError, "Contigous numpy array of doubles expected."
           " Make sure the numpy array is contiguous, and uses dtype=np.float_.");

  // Check size of passed array
  if ( PyArray_SIZE(xa) != arg1->num_states() )
    SWIG_exception(SWIG_ValueError, "Expected a numpy array of the same size "
		   "as number of states.");
  
  $1 = (double *)PyArray_DATA(xa);
}

%typemap(in) double* monitored
{
  // Check type
  if (!PyArray_Check($input))
    SWIG_exception(SWIG_TypeError, "Numpy array expected");

  // Get PyArrayObject
  PyArrayObject *xa = reinterpret_cast<PyArrayObject*>($input);

  // Check data type
  if (!(PyArray_ISCONTIGUOUS(xa) && PyArray_TYPE(xa) == NPY_DOUBLE))
    SWIG_exception(SWIG_TypeError, "Contigous numpy array of doubles expected."
           " Make sure the numpy array is contiguous, and uses dtype=np.float_.");

  // Check size of passed array
  if ( PyArray_SIZE(xa) != arg1->num_monitored() )
    SWIG_exception(SWIG_ValueError, "Expected a numpy array of the same size "
		   "as number of monitored intermediates.");
  
  $1 = (double *)PyArray_DATA(xa);
}

%typemap(in) double* field_params
{
  // Check type
  if (!PyArray_Check($input))
    SWIG_exception(SWIG_TypeError, "Numpy array expected");

  // Get PyArrayObject
  PyArrayObject *xa = reinterpret_cast<PyArrayObject*>($input);

  // Check data type
  if (!(PyArray_ISCONTIGUOUS(xa) && PyArray_TYPE(xa) == NPY_DOUBLE))
    SWIG_exception(SWIG_TypeError, "Contigous numpy array of doubles expected."
           " Make sure the numpy array is contiguous, and uses dtype=np.float_.");

  // Check size of passed array
  if ( PyArray_SIZE(xa) != arg1->num_field_parameters() )
    SWIG_exception(SWIG_ValueError, "Expected a numpy array of the same size "
		   "as number of field parameters.");
  
  $1 = (double *)PyArray_DATA(xa);
}

// Apply typemap to const
%apply double* field_params { const double* field_params }

// Apply states typemap to y (used in ODESolver::forward)
%apply const double* states { double* y }

// The typecheck
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY) double *
{
    $1 = PyArray_Check($input) ? 1 : 0;
}

%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY) const double *
{
    $1 = PyArray_Check($input) ? 1 : 0;
}

//-----------------------------------------------------------------------------
// Out typemap for const std::vector<std::string>
//-----------------------------------------------------------------------------
%typemap(out) const std::vector<std::string>&
{
  int size = $1->size();
  PyObject* ret = PyList_New(size);
  PyObject* tmp_Py_str = 0;
  for (int i=0; i < size; i++)
  {
    tmp_Py_str = PyString_FromString((*$1)[i].c_str());
    if (PyList_SetItem(ret, i, tmp_Py_str)<0)
    {
      PyErr_SetString(PyExc_ValueError,"something wrong happened when copying std::string to Python");
      return NULL;
    }
  }

  $result = ret;
}

//-----------------------------------------------------------------------------
// Out typemap for const std::vector<goss::uint>&
//-----------------------------------------------------------------------------
%typemap(out) const std::vector<goss::uint>&
{
  int size = $1->size();
  PyObject* ret = PyList_New(size);
  PyObject* tmp_Py_int = 0;
  for (int i=0; i < size; i++)
  {
    tmp_Py_int = PyInt_FromLong(static_cast<long>((*$1)[i]));
    if (PyList_SetItem(ret, i, tmp_Py_int) < 0)
    {
      PyErr_SetString(PyExc_ValueError, "something wrong happened when copying uint to Python");
      return NULL;
    }
  }

  $result = ret;
}



