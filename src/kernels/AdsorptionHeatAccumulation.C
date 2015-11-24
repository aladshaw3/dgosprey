//----------------------------------------
//  Created by Austin Ladshaw on 09/14/2015
//  Copyright (c) 2015
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

#include "AdsorptionHeatAccumulation.h"

template<>
InputParameters validParams<AdsorptionHeatAccumulation>()
{
	InputParameters params = validParams<Kernel>();
	params.addCoupledVar("solid_heats","Coupled variables for the solid heat of adsorption");
	return params;
}

AdsorptionHeatAccumulation::AdsorptionHeatAccumulation(const InputParameters & parameters) :
Kernel(parameters),
_porosity(getMaterialProperty<Real>("porosity")),
_pellet_density(getMaterialProperty<Real>("pellet_density"))
{
	unsigned int n = coupledComponents("solid_heats");
	_solid_heat.resize(n);
	_solid_heat_old.resize(n);
	
	for (unsigned int i = 0; i<_solid_heat.size(); ++i)
	{
		_solid_heat[i] = &coupledValue("solid_heats",i);
		_solid_heat_old[i] = &coupledValueOld("solid_heats",i);
	}
}

Real AdsorptionHeatAccumulation::computeQpResidual()
{
	double sum = 0.0;
	
	for (unsigned int i=0; i<_solid_heat.size(); i++)
	{
		sum = sum + (((*_solid_heat[i])[_qp]-(*_solid_heat_old[i])[_qp])/_dt);
	}
	return _pellet_density[_qp]*(1.0/1000.0)*sum*_test[_i][_qp];
}

Real AdsorptionHeatAccumulation::computeQpJacobian()
{
	return 0.0;
}

