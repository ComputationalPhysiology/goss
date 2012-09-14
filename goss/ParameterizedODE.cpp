#include <iostream>
#include <stdexcept>
#include "ParameterizedODE.h"


using namespace goss;

//-----------------------------------------------------------------------------
std::string ParameterizedODE::get_state_name(uint idx) const
{
  if (idx < _system_size)
    return _state_descr[idx];
  else
    throw std::runtime_error("Index out of range");
}
//-----------------------------------------------------------------------------
std::string ParameterizedODE::get_parameter_name(uint idx) const
{
  if (idx < _parameter_size)
    return _parameter_descr[idx];
  else
    throw std::runtime_error("Index out of range");
}
//-----------------------------------------------------------------------------


