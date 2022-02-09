#ifndef EulerRigidBody_h_IS_INCLUDED
#define EulerRigidBody_h_IS_INCLUDED

#include <memory>
#include <goss/ODE.h>

namespace goss
{

  class EulerRigidBody : public ODE
  {
  public:
    double I1, I2, I3, pi;
    EulerRigidBody() : ODE(3), I1(0.5), I2(2.0), I3(3.0), pi(2*acos(0.0))
    {}

    std::shared_ptr<ODE> copy() const
    {
      return std::make_shared<EulerRigidBody>(*this);
    }

    ~EulerRigidBody() {}

    void eval(const double* state, double t, double* f_vals)
    {
      const double y1 = state[0];
      const double y2 = state[1];
      const double y3 = state[2];
      f_vals[0] = (I2 - I3)/I1*y2*y3;
      f_vals[1] = (I3 - I1)/I2*y3*y1;
      f_vals[2] = (I1 - I2)/I3*y1*y2;

      if (t>=3*pi && t<=4*pi)
      {
	f_vals[2] += 0.25*pow(sin(t),2);
      }
    }

    void get_ic(goss::DoubleVector *res) const
    {
      res->n = _num_states;
      res->data.reset(new double[_num_states]);
      res->data[0] = 1.0;
      res->data[1] = 0.0;
      res->data[2] = 0.9;
    }
  };
}

#endif
