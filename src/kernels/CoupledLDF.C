/*!
 *  \file CoupledLDF.h
 *	\brief Advanced kernel for a cross coupled linear driving force mechanism
 *  \author Austin Ladshaw
 *	\date 01/29/2016
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

#include "CoupledLDF.h"
#include "Material.h"

template<>
InputParameters validParams<CoupledLDF>()
{
	InputParameters params = validParams<LinearDrivingForce>();
	params.addParam<Real>("driving_coef","Coefficient multiplied by the driving_variable");
	params.addCoupledVar("driving_var", "Coupled variable of the conserved quantity");
	return params;
}


CoupledLDF::CoupledLDF(const InputParameters & parameters)
:LinearDrivingForce(parameters),
_drive_coef(getParam<Real>("driving_coef")),
_drive_var(coupledValue("driving_var"))
{
}

Real CoupledLDF::computeQpResidual()
{
	_driving_value = _drive_coef * _drive_var[_qp];
	return LinearDrivingForce::computeQpResidual();
	
}

Real CoupledLDF::computeQpJacobian()
{
	return 0.0;
}

