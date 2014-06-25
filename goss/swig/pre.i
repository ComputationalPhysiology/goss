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

%ignore goss::DoubleVector::data;
%ignore goss::DoubleVector2D::data;

%ignore goss::ODESystemSolver::states;
%rename(_eval) goss::ODE::eval(const double*, double, double*);
%rename(eval_component) goss::ODE::eval(uint, const double*, double);

//-----------------------------------------------------------------------------
// Need to ignore these dues to SWIG confusion of overloaded functions
//-----------------------------------------------------------------------------
%ignore goss::Table::set(std::string,std::string,std::size_t);

//-----------------------------------------------------------------------------
// Ignore operators so SWIG stop complaining
//-----------------------------------------------------------------------------
%ignore goss::TableEntry::operator std::string;
%ignore goss::Progress::operator++;
%ignore goss::Progress::operator=;
%ignore goss::Table::operator=;
%ignore goss::TableEntry::operator=;

//-----------------------------------------------------------------------------
// Ignore GOSS C++ stream handling
//-----------------------------------------------------------------------------
%ignore goss::LogStream;
%ignore goss::cout;
%ignore goss::endl;

// ---------------------------------------------------------------------------
// Renames and ignores for Parameter
// For some obscure reason we need to rename Parameter
// ---------------------------------------------------------------------------
%rename (ParameterValue) goss::Parameter;
%rename (__nonzero__) goss::Parameter::operator bool;
%rename (__int__) goss::Parameter::operator int;
%rename (__float__) goss::Parameter::operator double;
%rename (__str__) goss::Parameter::operator std::string;
%rename (_assign) goss::Parameter::operator=;
%rename (_get_int_range) goss::Parameter::get_range(int& min_value, int& max_value) const;
%rename (_get_double_range) goss::Parameter::get_range(double& min_value, double& max_value) const;
%rename (_get_string_range) goss::Parameter::get_range(std::set<std::string>& range) const;
%ignore goss::Parameter::operator std::size_t;

// ---------------------------------------------------------------------------
// Renames and ignores for Parameters
// ---------------------------------------------------------------------------
%rename (_assign_bool) goss::Parameter::operator= (bool value);
%rename (_add) goss::Parameters::add;
%rename (_add_bool) goss::Parameters::add(std::string key, bool value);
%rename (_get_parameter_keys) goss::Parameters::get_parameter_keys;
%rename (_get_parameter_set_keys) goss::Parameters::get_parameter_set_keys;
%rename (_get_parameter_set) goss::Parameters::operator();
%rename (_get_parameter) goss::Parameters::operator[];
%rename (assign) goss::Parameters::operator=;
%ignore goss::Parameters::parse;
%ignore goss::Parameters::update;

// ---------------------------------------------------------------------------
// Typemaps (in) for std::set<std::string>
// ---------------------------------------------------------------------------
%typecheck(SWIG_TYPECHECK_STRING_ARRAY) std::set<std::string> {
    $1 = PySequence_Check($input) ? 1 : 0;
}

%typemap(in) std::set<std::string> (std::set<std::string> tmp) {
  int i;
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError,"expected a list of 'str'");
    return NULL;
  }
  int list_length = PyList_Size($input);
  if (!list_length > 0){
    PyErr_SetString(PyExc_ValueError,"expected a list with length > 0");
    return NULL;
  }
  for (i = 0; i < list_length; i++) {
    PyObject *o = PyList_GetItem($input,i);
    if (PyString_Check(o)) {
      tmp.insert(std::string(PyString_AsString(o)));
    } else {
      PyErr_SetString(PyExc_TypeError,"provide a list of strings");
      return NULL;
    }
  }
  $1 = tmp;
}

// ---------------------------------------------------------------------------
// Typemaps (argout) for std::vector<std::string>&
// ---------------------------------------------------------------------------
%typemap(in, numinputs=0) std::vector<std::string>& keys (std::vector<std::string> tmp_vec){
  $1 = &tmp_vec;
}

%typemap(argout) std::vector<std::string>& keys
{
  int size = $1->size();
  PyObject* ret = PyList_New(size);
  PyObject* tmp_Py_str = 0;
  for (int i=0; i < size; i++)
  {
    tmp_Py_str = PyString_FromString((*$1)[i].c_str());
    if (PyList_SetItem(ret,i,tmp_Py_str)<0)
    {
      PyErr_SetString(PyExc_ValueError,"something wrong happened when copying std::string to Python");
      return NULL;
    }
  }
  $result = ret;
}

// ---------------------------------------------------------------------------
// Typemaps (argout) for int &min_value, int &max_value
// ---------------------------------------------------------------------------
%typemap(in, numinputs=0) (int &min_value, int &max_value) (int min_temp, int max_temp){
  $1 = &min_temp; $2 = &max_temp;
}

%typemap(argout) (int &min_value, int &max_value)
{
  $result = Py_BuildValue("ii", *$1, *$2);
}

// ---------------------------------------------------------------------------
// Typemaps (argout) for real &min_value, real &max_value
// ---------------------------------------------------------------------------
%typemap(in, numinputs=0) (double &min_value, double &max_value) ( double min_temp, double max_temp){
  $1 = &min_temp; $2 = &max_temp;
}

%typemap(argout) (double &min_value, double &max_value)
{
  $result = Py_BuildValue("dd", *$1, *$2);
}

// ---------------------------------------------------------------------------
// Typemaps (argout) for std::set<std::string>&
// ---------------------------------------------------------------------------
%typemap(in, numinputs=0) std::set<std::string>& range (std::set<std::string> tmp_set){
  $1 = &tmp_set;
}

%typemap(argout) std::set<std::string>& range
{
  int size = $1->size();
  PyObject* ret = PyList_New(size);
  PyObject* tmp_Py_str = 0;
  std::set<std::string>::iterator it;
  int i = 0;
  for ( it=$1->begin() ; it != $1->end(); it++ )
  {
    tmp_Py_str = PyString_FromString(it->c_str());
    if (PyList_SetItem(ret, i, tmp_Py_str)<0)
    {
      PyErr_SetString(PyExc_ValueError,"something wrong happened when copying std::string to Python");
      return NULL;
    }
    i++;
  }
  $result = ret;
}
