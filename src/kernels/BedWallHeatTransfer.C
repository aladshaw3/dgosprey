#include "BedWallHeatTransfer.h"

template<>
InputParameters validParams<BedWallHeatTransfer>()
{
  InputParameters params = validParams<Kernel>();
  params.addCoupledVar("coupled", "Name of the column temperature variable to be coupled.");
  return params;
}


BedWallHeatTransfer::BedWallHeatTransfer(const InputParameters & parameters)
:Kernel(parameters),
_bed_wall_transfer_coeff(getMaterialProperty<Real>("bed_wall_transfer_coeff")),
_inner_dia(getMaterialProperty<Real>("inner_dia")),
_outer_dia(getMaterialProperty<Real>("outer_dia")),
_column_temp(coupledValue("coupled"))
{
}

Real BedWallHeatTransfer::computeQpResidual()
{
  //Note: if _inner_dia >= _outer_dia, this will lead to an error
  Real _coef = (4.0 * _bed_wall_transfer_coeff[_qp] * _inner_dia[_qp]) / ( (_outer_dia[_qp]*_outer_dia[_qp]) - (_inner_dia[_qp]*_inner_dia[_qp]) );
  return -_test[_i][_qp] * _coef * (_column_temp[_qp] - _u[_qp]);
}

Real BedWallHeatTransfer::computeQpJacobian()
{
  Real _coef = (4.0 * _bed_wall_transfer_coeff[_qp] * _inner_dia[_qp]) / ( (_outer_dia[_qp]*_outer_dia[_qp]) - (_inner_dia[_qp]*_inner_dia[_qp]) );
  
  return _test[_i][_qp] * _coef * _phi[_j][_qp];
}

