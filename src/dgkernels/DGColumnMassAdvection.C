//
//  DGColumnMassAdvection.C
//  DGOSPREY
//
//  Created by aladshaw3 on 1/15/15.
//
//

#include "DGColumnMassAdvection.h"

template<>
InputParameters validParams<DGColumnMassAdvection>()
{
	InputParameters params = validParams<DGAdvection>();
	return params;
}

DGColumnMassAdvection::DGColumnMassAdvection(const InputParameters & parameters) :
DGAdvection(parameters),
_vel(getMaterialProperty<Real>("velocity"))
{

}

Real DGColumnMassAdvection::computeQpResidual(Moose::DGResidualType type)
{
	
	_velocity(0)=0.0;
	_velocity(1)=_vel[_qp];
	_velocity(2)=0.0;
	
	return DGAdvection::computeQpResidual(type);
}

Real DGColumnMassAdvection::computeQpJacobian(Moose::DGJacobianType type)
{
	
	_velocity(0)=0.0;
	_velocity(1)=_vel[_qp];
	_velocity(2)=0.0;
	
	return DGAdvection::computeQpJacobian(type);
}
