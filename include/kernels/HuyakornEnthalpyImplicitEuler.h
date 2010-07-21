
#include "ImplicitEuler.h"

#ifndef HUYAKORNENTHALPYIMPLICITEULER
#define HUYAKORNENTHALPYIMPLICITEULER

//Forward Declarations
class HuyakornEnthalpyImplicitEuler;

template<>
InputParameters validParams<HuyakornEnthalpyImplicitEuler>();

class HuyakornEnthalpyImplicitEuler : public ImplicitEuler
{
public:

  HuyakornEnthalpyImplicitEuler(std::string name, MooseSystem & moose_system, InputParameters parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

//  VariableValue & _temperature;
//  VariableValue & _temperature_old;
//  VariableValue & _rho;
//  VariableValue & _rho_old;
  
  MaterialProperty<Real> & _temperature;
  MaterialProperty<Real> & _temperature_old;
  MaterialProperty<Real> & _rho;
  MaterialProperty<Real> & _rho_old;

  MaterialProperty<Real> & _rho_r;
  MaterialProperty<Real> & _porosity;
  
};
#endif //HUYAKORNENTHALPYIMPLICITEULER
