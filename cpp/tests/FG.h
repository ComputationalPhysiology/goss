#ifndef FG_h_IS_INCLUDED
#define FG_h_IS_INCLUDED

#include <memory>
#include <goss/ODE.h>
#include <cmath>

namespace goss
{

  /*
     This is a test equation presented in a Paper by Fox and Goodwin

     L. Fox, E.R. Goodwun,
     "Some new methods for the numerical integration
     of ordinary differential equations",
     Proc. Cambridge Philos. Soc. 45 (1949) 373-388

     The ode is

     y'=-10y+6z
     z'=13.5y-10z

     The exact solution is

     y(t)=2/3*exp(exp(-t)+exp(-19t))
     z(t)=exp(exp(-t)-exp(-19t))

     This is a stiff ode, which is not easy to integrate using schemes whith
     poor damping properties. The existance of the exact solution makes it
     attractive as a test equation.
   */
  class FG : public ODE
  {
  public:

    FG () : ODE(2){}

    std::shared_ptr<ODE> copy() const
    {
      return std::make_shared<FG>(*this);
    }

    ~FG() {};

    void eval(const double* y, double t, double* f_vals)
    {
      const double y1 = y[0];
      const double z1 = y[1];
      f_vals[0] = -10*y1+6*z1;
      f_vals[1] = 13.5*y1-10*z1;
    }

    void get_ic(goss::DoubleVector *res) const
    {
      res->n = _num_states;
      res->data.reset(new double[_num_states]);
      exact(res->data.get(), 0.0);
    }

    void exact(double* y, double t) const
    {
      y[0] = 2.0/3.0*std::exp(1)*(std::exp(-t)+std::exp(-19*t));
      y[1] = std::exp(1)*(std::exp(-t)-std::exp(-19*t));
    }

  //protected:
  //
  //  virtual void add_descr()
  //  {
  //    state_descr[0] = "y(t) = 2/3*exp*(exp(-t) + exp(-19t))";
  //    state_descr[1] = "z(t) =     exp*(exp(-t) - exp(-19t))";
  //  }
  };
}

#endif
