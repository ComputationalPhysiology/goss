//#include <stdlib.h>
#include <iostream>
#include <stdexcept>
#include "ODE.h"


using namespace goss;

//-----------------------------------------------------------------------------
double ODE::eval(uint idx, const double* state, double t) 
{ 
  std::cout << "Warning: Calling base class ODE::eval component wise. "\
    "This is very slow." << std::endl;

  if (idx < 0 || idx > _system_size-1)
    throw std::runtime_error("Index out of range");

  double* values = new double[_system_size];
  eval(state, t, values);
  
  const double ret = values[idx];
  delete[] values;

  return ret;
}
//-----------------------------------------------------------------------------
std::string ODE::get_state_name(uint idx) const
{
  if (idx >=0 && idx < _system_size)
    return _state_descr[idx];
  else
    throw std::runtime_error("Index out of range");
}
//-----------------------------------------------------------------------------
std::string ODE::get_parameter_name(uint idx) const
{
  if (idx >=0 && idx < _parameter_size)
    return _parameter_descr[idx];
  else
    throw std::runtime_error("Index out of range");
}
//-----------------------------------------------------------------------------


