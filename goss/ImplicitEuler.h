#ifndef ImplicitEuler_h_IS_INCLUDED
#define ImplicitEuler_h_IS_INCLUDED

#include <boost/scoped_array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
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
    ImplicitEuler (boost::shared_ptr<ODE> ode, double ldt=-1.0);

    // Constructor
    ImplicitEuler(double ldt);

    // Copy constructor
    ImplicitEuler(const ImplicitEuler& solver);

    // Return a copy of itself
    boost::shared_ptr<ODESolver> copy() const 
    { return boost::make_shared<ImplicitEuler>(*this); }

    // Constructor
    ~ImplicitEuler ();

    // Attach ODE
    virtual void attach(boost::shared_ptr<ODE> ode);

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
