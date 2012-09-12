#ifndef DOUBLEVECTOR_H_IS_INCLUDED
#define DOUBLEVECTOR_H_IS_INCLUDED

namespace goss {

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
