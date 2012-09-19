#ifndef DOUBLEVECTOR_H_IS_INCLUDED
#define DOUBLEVECTOR_H_IS_INCLUDED

#include <boost/scoped_array.hpp>

namespace goss {

  // A small wrapper around a double pointer
  struct DoubleVector
  {

    boost::scoped_array<double> data;
    int n;

  };

  // A small wrapper around a double pointer
  struct DoubleVector2D 
  {

    boost::scoped_array<double> data;
    int m;
    int n;

  };

}

#endif
