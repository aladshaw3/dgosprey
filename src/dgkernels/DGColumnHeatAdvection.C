//
//  DGColumnHeatAdvection.C
//  DGOSPREY
//
//  Created by aladshaw3 on 1/16/15.
//
//

#include "DGColumnHeatAdvection.h"

template<>
InputParameters validParams<DGColumnHeatAdvection>()
{
	InputParameters params = validParams<DGAdvection>();
	return params;
}

DGColumnHeatAdvection::DGColumnHeatAdvection(const InputParameters & parameters) :
DGAdvection(parameters),
_vel(getMaterialProperty<Real>("velocity")),
_gas_density(getMaterialProperty<Real>("gas_density")),
_gas_heat_capacity(getMaterialProperty<Real>("gas_heat_capacity"))
{

}

Real DGColumnHeatAdvection::computeQpResidual(Moose::DGResidualType type)
{
	
	_velocity(0)=0.0;
	_velocity(1)=_vel[_qp]*_gas_heat_capacity[_qp]*_gas_density[_qp];
	_velocity(2)=0.0;
	
	return DGAdvection::computeQpResidual(type);
}

Real DGColumnHeatAdvection::computeQpJacobian(Moose::DGJacobianType type)
{
	
	_velocity(0)=0.0;
	_velocity(1)=_vel[_qp]*_gas_heat_capacity[_qp]*_gas_density[_qp];
	_velocity(2)=0.0;
	
	return DGAdvection::computeQpJacobian(type);
}
