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

#include "TimeDerivative.h"

#ifndef PTMASSTIMEDERIVATIVE_FULL_H
#define PTMASSTIMEDERIVATIVE_FULL_H

//Forward Declarations
class PTMassTimeDerivative_Full;

template<>
InputParameters validParams<PTMassTimeDerivative_Full>();

class PTMassTimeDerivative_Full : public TimeDerivative
{
  public:

    PTMassTimeDerivative_Full(const InputParameters & parameters);

  protected:

    virtual Real computeQpResidual();
    virtual Real computeQpJacobian();
    virtual Real computeQpOffDiagJacobian(unsigned int jvar);
    
    bool _has_coupled_temp;

    const MaterialProperty<Real> & _poro;
    const MaterialProperty<Real> & _wrho;
    const MaterialProperty<Real> & _drop;
    const MaterialProperty<Real> & _drot;
    
    const VariableValue & _dTdt;
    const VariableValue & _dTdtdT;
    
  private:
    unsigned int _temp_var;
};
#endif
