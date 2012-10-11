#ifndef DOUBLEVECTOR_H_IS_INCLUDED
#define DOUBLEVECTOR_H_IS_INCLUDED

#include <boost/scoped_array.hpp>

#include "types.h"

namespace goss {

  // A small wrapper around a double pointer
  struct DoubleVector
  {

    boost::scoped_array<double> data;
    uint n;

  };

  // A small wrapper around a double pointer
  struct DoubleVector2D 
  {

    boost::scoped_array<double> data;
    uint m;
    uint n;

  };

}

#endif
