/*!
 *  \file Aux_LDF.h
 *	\brief Generic auxillary kernel to calculate the value of an aux variable using LDF kinetics
 *  \author Austin Ladshaw
 *	\date 02/04/2016
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

#include "Aux_LDF.h"

template<>
InputParameters validParams<Aux_LDF>()
{
	InputParameters params = validParams<AuxKernel>();
	params.addParam<Real>("ldf_coeff",1.0,"Value of the linear driving force coefficient");
	params.addParam<Real>("driving_value",1.0,"Value that the driving force is pushing the aux variable towards.");
	return params;
}

Aux_LDF::Aux_LDF(const InputParameters & parameters) :
AuxKernel(parameters),
_ldf_coef(getParam<Real>("ldf_coeff")),
_driving_value(getParam<Real>("driving_value"))
{
}

Real Aux_LDF::computeValue()
{
	if (_dt == 0.0)
	{
		if (_ldf_coef <= 1.0 && _ldf_coef >= 0.0)
			return (_u[_qp]*(1.0 - _ldf_coef)) + (_ldf_coef*_driving_value);
		else
		{
			if (_ldf_coef*_driving_value >= _driving_value)
				return _driving_value;
			else
				return _ldf_coef*_driving_value;
		}
	}
	else
	{
		return (_u[_qp] + (_dt*_ldf_coef*_driving_value))/(1.0 + (_dt*_ldf_coef));
	}
}
