/*!
 *  \file MAGPIE_ConstLDF_Adsorption.h
 *	\brief Auxillary kernel to calculate adsorption based on LDF kinetics with constant coefficients
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

#include "MAGPIE_ConstLDF_Adsorption.h"

template<>
InputParameters validParams<MAGPIE_ConstLDF_Adsorption>()
{
	InputParameters params = validParams<Aux_LDF>();
	params.addParam<unsigned int>("index",0,"Index of the species that we are interested in.");
	return params;
}

MAGPIE_ConstLDF_Adsorption::MAGPIE_ConstLDF_Adsorption(const InputParameters & parameters) :
Aux_LDF(parameters),
_index(getParam<unsigned int>("index")),
_magpie_dat(getMaterialProperty< MAGPIE_DATA >("magpie_data"))
{
}

Real MAGPIE_ConstLDF_Adsorption::computeValue()
{
	MAGPIE_DATA magpie_copy;
	magpie_copy = _magpie_dat[_qp];
	
	//Call MAGPIE Simulation for Unperturbed data
	if (_magpie_dat[_qp].gsta_dat[_index].qmax > 0.0)
	{
		int success = 0;
		success = MAGPIE( (void *)&magpie_copy );
		if (success < 0 || success > 3)
		{
			for (int i=0; i<magpie_copy.sys_dat.N; i++)
			{
				if (i != _index)
				{
					magpie_copy.gpast_dat[i].y = 0.0;
					magpie_copy.sys_dat.Carrier = true;
				}
			}
			success = MAGPIE( (void *)&magpie_copy );
			if (success < 0 || success > 3)
			{
				mError(simulation_fail);
			}
		}
		else success = 0;
		
		_driving_value = magpie_copy.gpast_dat[_index].q;
		
	}
	else
	{
		_driving_value = 0.0;
	}
	
	return Aux_LDF::computeValue();
}