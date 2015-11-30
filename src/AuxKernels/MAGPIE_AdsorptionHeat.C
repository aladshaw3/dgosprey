/*!
 *  \file MAGPIE_AdsorptionHeat.h
 *	\brief Auxillary kernel to calculate heat of adsorption of a particular gas species in the system
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

#include "MAGPIE_AdsorptionHeat.h"

template<>
InputParameters validParams<MAGPIE_AdsorptionHeat>()
{
	InputParameters params = validParams<AuxKernel>();
	params.addParam<unsigned int>("index",0,"Index of the species that we are interested in.");
	params.addCoupledVar("solid_conc","Coupled variable for the solid concentration of interest");
	return params;
}

MAGPIE_AdsorptionHeat::MAGPIE_AdsorptionHeat(const InputParameters & parameters) :
AuxKernel(parameters),
_index(getParam<unsigned int>("index")),
_magpie_dat(getMaterialProperty< MAGPIE_DATA >("magpie_data")),
_solid_conc(coupledValue("solid_conc"))
{
}

Real
MAGPIE_AdsorptionHeat::computeValue()
{
	MAGPIE_DATA magpie_copy;
	magpie_copy = _magpie_dat[_qp];
	
	//Call MAGPIE Simulation for Unperturbed data
	if (_magpie_dat[_qp].gsta_dat[_index].qmax > 0.0)
	{
		//int success = 0;
		//success = MAGPIE( (void *)&magpie_copy );
		//if (success < 0 || success > 3) {mError(simulation_fail);}
		//else success = 0;
		
		double pi = magpie_copy.gpast_dat[_index].y * magpie_copy.sys_dat.PT;
		return _solid_conc[_qp] * Qst(pi,(void *)&magpie_copy,_index);
	}
	else
	{
		return 0.0;
	}
}
