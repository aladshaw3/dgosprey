/*!
 *  \file MAGPIE_Perturbation.h
 *	\brief Auxillary kernel to calculate the perturbed adsorption equilibria of a particular gas species in the system
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

#include "MAGPIE_Perturbation.h"

template<>
InputParameters validParams<MAGPIE_Perturbation>()
{
	InputParameters params = validParams<AuxKernel>();
	params.addParam<unsigned int>("index",0,"Index of the species that we are interested in.");
	return params;
}

MAGPIE_Perturbation::MAGPIE_Perturbation(const InputParameters & parameters) :
AuxKernel(parameters),
_index(getParam<unsigned int>("index")),
_magpie_dat(getMaterialProperty< MAGPIE_DATA >("magpie_data"))
{
	//Forces specific execution behavior of the auxkernel
	_exec_flags.clear();
	_exec_flags.push_back(EXEC_INITIAL);
	_exec_flags.push_back(EXEC_TIMESTEP_END);
	//_exec_flags.push_back(EXEC_TIMESTEP_BEGIN);
}

Real
MAGPIE_Perturbation::computeValue()
{
	MAGPIE_DATA magpie_copy;
	magpie_copy = _magpie_dat[_qp];
	
	
	//Check for adsorption
	if (_magpie_dat[_qp].gsta_dat[_index].qmax > 0.0)
	{
		//perturn the copy's _index y
		double pi = _magpie_dat[_qp].gpast_dat[_index].y * _magpie_dat[_qp].sys_dat.PT;
		double ci = Cstd(pi,_magpie_dat[_qp].sys_dat.T) + sqrt(DBL_EPSILON);
		double yi = Pstd(ci,_magpie_dat[_qp].sys_dat.T) / _magpie_dat[_qp].sys_dat.PT;
		magpie_copy.gpast_dat[_index].y = yi;
	
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
		
		return magpie_copy.gpast_dat[_index].q;
		
	}
	else
	{
		return 0.0;
	}
}

