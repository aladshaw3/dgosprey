#include "WallAmbientHeatTransfer.h"

template<>
InputParameters validParams<WallAmbientHeatTransfer>()
{
  InputParameters params = validParams<Kernel>();
  params.addCoupledVar("coupled", "Name of the ambient temperature variable to be coupled.");
  return params;
}


WallAmbientHeatTransfer::WallAmbientHeatTransfer(const InputParameters & parameters)
:Kernel(parameters),
_wall_exterior_transfer_coeff(getMaterialProperty<Real>("wall_exterior_transfer_coeff")),
_inner_dia(getMaterialProperty<Real>("inner_dia")),
_outer_dia(getMaterialProperty<Real>("outer_dia")),
_ambient_temp(coupledValue("coupled"))
{
}

Real WallAmbientHeatTransfer::computeQpResidual()
{
  //Note: if _inner_dia >= _outer_dia, this will lead to an error
  Real _coef = (4.0 * _wall_exterior_transfer_coeff[_qp] * _outer_dia[_qp]) / ( (_outer_dia[_qp]*_outer_dia[_qp]) - (_inner_dia[_qp]*_inner_dia[_qp]) );
  return _test[_i][_qp] * _coef * (_u[_qp] - _ambient_temp[_qp]);
}

Real WallAmbientHeatTransfer::computeQpJacobian()
{
  Real _coef = (4.0 * _wall_exterior_transfer_coeff[_qp] * _outer_dia[_qp]) / ( (_outer_dia[_qp]*_outer_dia[_qp]) - (_inner_dia[_qp]*_inner_dia[_qp]) );
  
  return _test[_i][_qp] * _coef * _phi[_j][_qp];
}

