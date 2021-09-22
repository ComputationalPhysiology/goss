
#include <memory>
#include <goss/goss.h>
#include <cmath>

using namespace goss;

class Sin : public ODE
{
public:
double omega;

Sin () : ODE(2), omega(4*std::acos(0.0)) {}

std::shared_ptr<ODE> copy() const
{
    return std::make_shared<Sin>(*this);
}

~Sin() {}

void eval(const double* y, double t, double* f_vals)
{
    double y1 = y[0];
    double y2 = y[1];
    f_vals[0] =  omega*y2; 
    f_vals[1] = -omega*y1;
}

void get_ic(goss::DoubleVector *res) const
{
    res->n = _num_states;
    res->data.reset(new double[_num_states]);
    res->data[0] = 0.0; 
    res->data[1] = 1.0;//omega; 
}
};

int main(){

    std::shared_ptr<Sin> ode;
    ExplicitEuler solver(ode);
    DoubleVector x;

    double dt = 0.0001;
    double tstop = 10.0;
    

    const uint nstep = std::ceil(tstop/dt - 1.0E-12);
    
    double t = 0.0;
    
    for (uint i = 0; i < nstep; i++)
    {
      solver.forward(x.data.get(), t, dt);
      
      t += dt;
      info("t %.2f: y: %.2f, z: %.2f", t, x.data.get()[0], x.data.get()[1]);
      /*
      if (std::fmod(t, .1) < dt)
	info("t %.2f: y: %.2f, z: %.2f", t, x.data.get()[0], x.data.get()[1]);
      */
    }

    return 0;
}