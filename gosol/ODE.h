#ifndef ODE_h_IS_INCLUDED
#define ODE_h_IS_INCLUDED

#include <vector>
#include <string>
#include <iostream>

namespace gosol {

  struct DoubleVector; 
  class ODE 
  {
    public:
      ODE(int system_size, int parameter_size) : system_size(system_size), parameter_size(parameter_size) { 
        state_descr.resize(system_size); 
        parameter_descr.resize(parameter_size); 
      } 

      virtual ~ODE() {}

      inline int size() const { return system_size; }
      virtual void eval(const double* state, double t, double* f_vals) = 0;
      virtual double eval(int idx, const double* state, double t);

      virtual int numLinearTerms() const {return 0;}
      virtual int numParameters() const {return 0;}
      virtual void linearTerms(int* indices) const {}
      virtual void linearDerivatives(const double* x, double t, double* y) const {}

      virtual void evalIntermediates(const double* x, double t, double* y){ /*TODO: Implement*/}
      virtual double evalIntermediate(int i, const double* x, double t){ return 0.0; /*TODO: Implement*/}


      /* Method for extracting data from the ODE. Must be implemented in
       * subclasses. FIXME: Figure out a better system for this.
       */
      virtual void probe (double* y) { }

      virtual void setParameters(const double* values) {}
      virtual void setParameter(int idx, double value) {}
      virtual double getParameter(int idx) const;
      virtual void getParameters(gosol::DoubleVector*) const;
      virtual std::string getStateName(int idx) const;
      virtual std::string getParameterName(int idx) const;
      virtual void getIC(gosol::DoubleVector *res) const = 0;
      virtual void getDefaultInitialCondition(gosol::DoubleVector *res) const
      {
        std::cout << "Warning: getDefaultInitialCondition() is deprecated. Use getIC() instead." << std::endl;
        getIC(res);
      }

    protected: 
      virtual void setStateDescriptions() {};
      virtual void setParameterNames() {};
      const int system_size;
      const int parameter_size;
      std::vector<std::string> state_descr;
      std::vector<std::string> parameter_descr;
  };
}

#endif
