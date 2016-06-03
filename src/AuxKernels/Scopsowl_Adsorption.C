/*!
 *  \file Scopsowl_Adsorption.h
 *	\brief Auxillary kernel to calculate adsorption kinetics of a particular gas species in the system
 *  \author Austin Ladshaw
 *	\date 06/03/2016
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

#include "Scopsowl_Adsorption.h"

template<>
InputParameters validParams<Scopsowl_Adsorption>()
{
	InputParameters params = validParams<AuxKernel>();
	params.addParam<unsigned int>("index",0,"Index of the species that we are interested in.");
	return params;
}

Scopsowl_Adsorption::Scopsowl_Adsorption(const InputParameters & parameters) :
AuxKernel(parameters),
_index(getParam<unsigned int>("index")),
_owl_dat(getMaterialProperty< SCOPSOWL_DATA >("owl_data"))
{
	//Forces specific execution behavior of the auxkernel
	_exec_flags.clear();
	_exec_flags.push_back(EXEC_INITIAL);
	_exec_flags.push_back(EXEC_TIMESTEP_END);
}

Real
Scopsowl_Adsorption::computeValue()
{
	Real q = 0.0;
	
	return q;
}