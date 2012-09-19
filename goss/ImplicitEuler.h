#ifndef ImplicitEuler_h_IS_INCLUDED
#define ImplicitEuler_h_IS_INCLUDED

#include <boost/scoped_array.hpp>
#include <vector>

#include "ImplicitODESolver.h"

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

    // Copy constructor
    ImplicitEuler(const ImplicitEuler& solver);

    // Constructor
    ~ImplicitEuler ();

    // Attach ODE
    virtual void attach(ODE* ode);

    // Reset ODE
    virtual void reset();

    // Step solver an interval of time forward
    void forward(double* y, double t, double interval);

    std::vector<int> newton_iter1;
    std::vector<int> newton_accepted1;
    std::vector<double> dt_v;

  protected:

    boost::scoped_array<double> z1;
    bool justrefined;

  };

}
#endif
