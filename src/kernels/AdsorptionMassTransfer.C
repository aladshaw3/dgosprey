/*!
 *  \file AdsorptionMassTransfer.h
 *	\brief Standard kernel for the transfer of mass via adsorption
 *  \author Austin Ladshaw
 *	\date 01/29/2015
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


#include "AdsorptionMassTransfer.h"

template<>
InputParameters validParams<AdsorptionMassTransfer>()
{
	InputParameters params = validParams<Kernel>();
	params.addCoupledVar("solid_conc","Coupled variables for the solid adsorption concentrations");
	return params;
}

AdsorptionMassTransfer::AdsorptionMassTransfer(const InputParameters & parameters) :
Kernel(parameters),
_porosity(getMaterialProperty<Real>("porosity")),
_pellet_density(getMaterialProperty<Real>("pellet_density")),
_solid(coupledValue("solid_conc")),
_solid_old(coupledValueOld("solid_conc"))
{

}

Real AdsorptionMassTransfer::computeQpResidual()
{
	return -(1.0-_porosity[_qp])*_pellet_density[_qp]*((_solid[_qp]-_solid_old[_qp])/_dt)*_test[_i][_qp];
}

Real AdsorptionMassTransfer::computeQpJacobian()
{
	return 0.0;
}

