#ifndef ImplicitEuler_h_IS_INCLUDED
#define ImplicitEuler_h_IS_INCLUDED

#include "ImplicitODESolver.h"
#include <vector>

namespace goss 
{

  // Implicit Euler
  class ImplicitEuler : public ImplicitODESolver
  {
  public:

    // Default constructor
    ImplicitEuler();

    // Constructor
    ImplicitEuler (ODE* ode, double ldt=-1.0);

    // Constructor
    ImplicitEuler(double ldt);

    // Constructor
    ~ImplicitEuler ();

    // Attach ODE
    virtual void attach(ODE* ode);

    // Step solver an interval of time forward
    void forward(double* y, double t, double interval);

    std::vector<int> newton_iter1;
    std::vector<int> newton_accepted1;
    std::vector<double> dt_v;

  protected:

    double* z1;
    bool justrefined;

  };

}
#endif
