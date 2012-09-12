#ifndef DOUBLEVECTORS_H_IS_INCLUDED
#define DOUBLEVECTORS_H_IS_INCLUDED

namespace pycc {

  struct DoubleVector {
    double *data;
    int n;
  };

  struct DoubleVector2D {
    double *data;
    int m;
    int n;
  };

}

#endif
