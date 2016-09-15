/****************************************************************/
/*             DO NOT MODIFY OR REMOVE THIS HEADER              */
/*          FALCON - Fracturing And Liquid CONvection           */
/*                                                              */
/*       (c)     2012 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

/*----------------------------------------------------------------

Contributors: Yidong Xia (INL)

Descriptions: time derivative of mass balance equation
              in coupled pressure-temperature based
              T-H-M-C balance equations

----------------------------------------------------------------*/

#include "PTMassTimeDerivative_Full.h"

template<>
InputParameters validParams<PTMassTimeDerivative_Full>()
{
  InputParameters params = validParams<TimeDerivative>();
    
  params.addRequiredCoupledVar("coupled_temperature", " Use Coupled temperature to calculate the time derivative");
  return params;
}


/*******************************************************************************
Routine: PTMassTimeDerivative_Full -- constructor
Authors: Yidong Xia
*******************************************************************************/
PTMassTimeDerivative_Full::PTMassTimeDerivative_Full(const InputParameters & parameters):
  TimeDerivative(parameters),
  _has_coupled_temp(isCoupled("coupled_temperature")),
  _poro(getMaterialProperty<Real>("porosity")),
  _wrho(getMaterialProperty<Real>("density_water")),
  _drop(getMaterialProperty<Real>("partial_rho_over_partial_pres")),
  _drot(getMaterialProperty<Real>("partial_rho_over_partial_temp")),
  _dTdt(coupledDot("coupled_temperature")),
  _dTdtdT(coupledDotDu("coupled_temperature")),
  _temp_var(_has_coupled_temp ? coupled("coupled_temperature") : zero)
{}


/*******************************************************************************
Routine: computeQpResidual -- compute residual at quadrature point
Authors: Yidong Xia
*******************************************************************************/
Real
PTMassTimeDerivative_Full::
computeQpResidual()
{
  return _poro[_qp]*_drop[_qp]*TimeDerivative::computeQpResidual() + _test[_i][_qp]*_poro[_qp]*_drot[_qp]*_dTdt[_qp];
}


/*******************************************************************************
Routine: computeQpJacobian -- compute Jacobian at quadrature point
Authors: Yidong Xia
*******************************************************************************/
Real
PTMassTimeDerivative_Full::
computeQpJacobian()
{
  return _poro[_qp]*_drop[_qp]*TimeDerivative::computeQpJacobian();
}

Real
PTMassTimeDerivative_Full::
computeQpOffDiagJacobian(unsigned int jvar)
{
    Real r=0.0;
    
    if (jvar == _temp_var && _has_coupled_temp)
    
    r= _test[_i][_qp]*_poro[_qp]*_drot[_qp]*_phi[_j][_qp]*_dTdtdT[_qp];
    
    return r;
}




