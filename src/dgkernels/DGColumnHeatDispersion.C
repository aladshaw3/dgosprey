//
//  DGColumnHeatDispersion.C
//  DGOSPREY
//
//  Created by aladshaw3 on 1/16/15.
//
//

#include "DGColumnHeatDispersion.h"

template<>
InputParameters validParams<DGColumnHeatDispersion>()
{
	InputParameters params = validParams<DGAnisotropicDiffusion>();
	return params;
}

DGColumnHeatDispersion::DGColumnHeatDispersion(const InputParameters & parameters) :
DGAnisotropicDiffusion(parameters),
_conductivity(getMaterialProperty<Real>("conductivity"))
{
	
}

Real DGColumnHeatDispersion::computeQpResidual(Moose::DGResidualType type)
{
	_Diffusion(0,0) =  _conductivity[_qp];
	_Diffusion(0,1) = std::pow(std::pow(_conductivity[_qp],2.0) + std::pow(_conductivity[_qp],2.0),0.5);
	_Diffusion(0,2) = 0.0;
	
	_Diffusion(1,0) = std::pow(std::pow(_conductivity[_qp],2.0) + std::pow(_conductivity[_qp],2.0),0.5);
	_Diffusion(1,1) = _conductivity[_qp];
	_Diffusion(1,2) = 0.0;
	
	_Diffusion(2,0) = 0.0;
	_Diffusion(2,1) = 0.0;
	_Diffusion(2,2) = 0.0;
	
	return DGAnisotropicDiffusion::computeQpResidual(type);
}

Real DGColumnHeatDispersion::computeQpJacobian(Moose::DGJacobianType type)
{
	
	_Diffusion(0,0) =  _conductivity[_qp];
	_Diffusion(0,1) = std::pow(std::pow(_conductivity[_qp],2.0) + std::pow(_conductivity[_qp],2.0),0.5);
	_Diffusion(0,2) = 0.0;
	
	_Diffusion(1,0) = std::pow(std::pow(_conductivity[_qp],2.0) + std::pow(_conductivity[_qp],2.0),0.5);
	_Diffusion(1,1) = _conductivity[_qp];
	_Diffusion(1,2) = 0.0;
	
	_Diffusion(2,0) = 0.0;
	_Diffusion(2,1) = 0.0;
	_Diffusion(2,2) = 0.0;
	
	return DGAnisotropicDiffusion::computeQpJacobian(type);
}
