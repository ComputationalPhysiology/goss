#include <stdlib.h>
#include <iostream>
#include <stdexcept>
#include "ODE.h"


using namespace goss;

//-----------------------------------------------------------------------------
double ODE::eval(int idx, const double* state, double t) 
//-----------------------------------------------------------------------------
{ 
  std::cout << "Warning: Calling base class ODE::eval component wise. This is very slow." << std::endl;
  if (idx < 0 || idx > system_size-1)
    throw std::runtime_error("Index out of range");
  double* values = static_cast<double*>(malloc(sizeof(double)*system_size));
  eval(state, t, values);
  return values[idx];
}

//-----------------------------------------------------------------------------
double ODE::getParameter(int idx) const
//-----------------------------------------------------------------------------
{

}
//-----------------------------------------------------------------------------
void ODE::getParameters(goss::DoubleVector*) const
//-----------------------------------------------------------------------------
{

}

//-----------------------------------------------------------------------------
std::string ODE::getStateName(int idx) const
//-----------------------------------------------------------------------------
{
  if (idx >=0 && idx < system_size)
    return state_descr[idx];
  else
    throw std::runtime_error("Index out of range");
}

//-----------------------------------------------------------------------------
std::string ODE::getParameterName(int idx) const
//-----------------------------------------------------------------------------
{
  if (idx >=0 && idx < parameter_size)
    return parameter_descr[idx];
  else
    throw std::runtime_error("Index out of range");
}


