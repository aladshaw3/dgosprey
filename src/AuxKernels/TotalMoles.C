/*!
 *  \file TotalMoles.h
 *	\brief Aux Kernel to calculate total moles of a gas in the bed.
 *	\details This file is responsible for calculating the total moles of a gas in the bed.
 *
 *  \author Austin Ladshaw
 *	\date 01/24/2017
 *	\copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science and was developed for use
 *			   by Idaho National Laboratory and Oak Ridge National Laboratory
 *			   engineers and scientists. Portions Copyright (c) 2016, all
 *             rights reserved.
 *
 *			   Austin Ladshaw does not claim any ownership or copyright to the
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

#include "TotalMoles.h"

template<>
InputParameters validParams<TotalMoles>()
{
	InputParameters params = validParams<AuxKernel>();
	params.addCoupledVar("solid","Solid concentration variable being coupled");
	params.addCoupledVar("gas", "Gas concentration variable being coupled");
	return params;
}

TotalMoles::TotalMoles(const InputParameters & parameters) :
AuxKernel(parameters),
_porosity(getMaterialProperty<Real>("porosity")),
_pellet_density(getMaterialProperty<Real>("pellet_density")),
_bed_length(getMaterialProperty<Real>("bed_length")),
_inner_dia(getMaterialProperty<Real>("inner_dia")),
_solid(coupledValue("solid")),
_gas(coupledValue("gas"))
{
	
}

Real
TotalMoles::computeValue()
{
	Real moles = 0.0;
	
	//double volume = M_PI * _inner_dia[_qp] * _inner_dia[_qp] / 4.0 * _bed_length[_qp];
	//moles = (volume/1000.0) * ((_porosity[_qp] * _gas[_qp]) + ((1.0-_porosity[_qp])*_pellet_density[_qp]*_solid[_qp]));
	moles = ((_porosity[_qp] * _gas[_qp]) + ((1.0-_porosity[_qp])*_pellet_density[_qp]*_solid[_qp]));
	
	return moles;
}
