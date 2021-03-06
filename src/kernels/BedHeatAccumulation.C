/*!
 *  \file BedHeatAccumulation.h
 *	\brief Time Derivative kernel for the accumulation of heat in a fixed-bed column
 *  \author Austin Ladshaw
 *	\date 11/20/2015
 *	\copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science and was developed for use
 *			   by Idaho National Laboratory and Oak Ridge National Laboratory
 *			   engineers and scientists. Portions Copyright (c) 2015, all
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

#include "BedHeatAccumulation.h"

template<>
InputParameters validParams<BedHeatAccumulation>()
{
  InputParameters params = validParams<TimeDerivative>();
  return params;
}


BedHeatAccumulation::BedHeatAccumulation(const InputParameters & parameters)
:TimeDerivative(parameters),
_heat_retardation(getMaterialProperty<Real>("heat_retardation"))
{
  
}

Real
BedHeatAccumulation::computeQpResidual()
{
  return _heat_retardation[_qp] * TimeDerivative::computeQpResidual();
}

Real
BedHeatAccumulation::computeQpJacobian()
{
  return _heat_retardation[_qp] * TimeDerivative::computeQpJacobian();
}
