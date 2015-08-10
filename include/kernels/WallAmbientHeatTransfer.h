#include "Kernel.h"

#ifndef WALLAMBIENTHEATTRANSFER_H
#define WALLAMBIENTHEATTRANSFER_H

class WallAmbientHeatTransfer;

template<>
InputParameters validParams<WallAmbientHeatTransfer>();

class WallAmbientHeatTransfer : public Kernel
{
public:
  WallAmbientHeatTransfer(const InputParameters & parameters);
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  
private:
  const MaterialProperty<Real> & _wall_exterior_transfer_coeff;
  const MaterialProperty<Real> & _inner_dia;
  const MaterialProperty<Real> & _outer_dia;
  
  VariableValue & _ambient_temp;
  
};
#endif //WALLAMBIENTHEATTRANSFER_H