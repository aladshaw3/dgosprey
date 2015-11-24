#include "WallHeatAccumulation.h"

template<>
InputParameters validParams<WallHeatAccumulation>()
{
  InputParameters params = validParams<TimeDerivative>();
  return params;
}


WallHeatAccumulation::WallHeatAccumulation(const InputParameters & parameters)
:TimeDerivative(parameters),
_wall_density(getMaterialProperty<Real>("wall_density")),
_wall_heat_capacity(getMaterialProperty<Real>("wall_heat_capacity"))
{
  
}

Real WallHeatAccumulation::computeQpResidual()
{
  return _wall_density[_qp] * _wall_heat_capacity[_qp] * TimeDerivative::computeQpResidual();
}

Real WallHeatAccumulation::computeQpJacobian()
{
  return _wall_density[_qp] * _wall_heat_capacity[_qp] * TimeDerivative::computeQpJacobian();
}
