#include <iostream>
#include <stdexcept>
#include "ODE.h"


using namespace goss;

//-----------------------------------------------------------------------------
double ODE::eval(uint idx, const double* state, double t) 
{ 
  std::cout << "Warning: Calling base class ODE::eval component wise. "\
    "This is very slow." << std::endl;

  if (idx >= _num_states)
    throw std::runtime_error("Index out of range");

  double* values = new double[_num_states];
  eval(state, t, values);
  
  const double ret = values[idx];
  delete[] values;

  return ret;
}
//-----------------------------------------------------------------------------
