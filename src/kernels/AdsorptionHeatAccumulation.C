/*!
 *  \file AdsorptionHeatAccumulation.h
 *	\brief Standard kernel for the heat of adsorption and its effect on the system temperature
 *  \author Austin Ladshaw
 *	\date 11/20/2015
 *	\copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science and was developed for use
 *			   by Idaho National Laboratory and Oak Ridge National Laboratory
 *			   engineers and scientists. Portions Copyright (c) 2015, all
 *             rights reserved.
 *
 *			   Austin Ladshaw does not claim any owership or copyright to the
 *			   MOOSE framework in which these kernels are constructed, only
 *			   the kernels themselves. The MOOSE framework copyright is held
 *			   by the Battelle Energy Alliance, LLC (c) 2010, all rights reserved.
 */

/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/


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

