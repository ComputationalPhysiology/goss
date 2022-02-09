#ifndef VDP_h_IS_INCLUDED
#define VDP_h_IS_INCLUDED

#include <memory>
#include <goss/ODE.h>
#include <stdexcept>

namespace goss
{

  class VDP : public ODE
  {
  public:
    double mu;

    VDP() : ODE(2), mu(10.0) {}

    virtual ~VDP() {}

    std::shared_ptr<ODE> copy() const
    {
      return std::make_shared<VDP>(*this);
    }

    void eval(const double* y, double t, double* f_vals)
    {
      const double y1 = y[0];
      const double y2 = y[1];
      f_vals[0] = y2;
      f_vals[1] = mu*mu*((1-y1*y1)*y2-y1);
    }

    virtual double eval(uint idx, const double* state, double t)
    {
      const double y1 = state[0];
      const double y2 = state[1];
      switch (idx)
      {
      case 0:
	return y2;
	break;
      case 1:
	return mu*mu*((1-y1*y1)*y2-y1);
	break;
      default:
	throw std::runtime_error("Index out of bounds");
	return 0.0;
      }
    }

    void get_ic(goss::DoubleVector *res) const
    {
      res->n = _num_states;
      res->data.reset(new double[_num_states]);
      res->data[0] = 2.0;
      res->data[1] = 0.0;
    }

    void linear_terms(uint* indices) const
    {
      indices[0] = 0;
      indices[1] = 1;
    }

    void linear_derivatives(const double* x, double t, double* y) const
    {
      const double x1 = x[0];
      y[1] = mu*mu*(1-x1*x1);
    }
  };
}

#endif
