/*!
 *  \file MAGPIE_Adsorption.h
 *	\brief Auxillary kernel to calculate adsorption equilibria of a particular gas species in the system
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

#include "MAGPIE_Adsorption.h"

template<>
InputParameters validParams<MAGPIE_Adsorption>()
{
	InputParameters params = validParams<AuxKernel>();
	params.addParam<unsigned int>("index",0,"Index of the species that we are interested in.");
	return params;
}

MAGPIE_Adsorption::MAGPIE_Adsorption(const InputParameters & parameters) :
AuxKernel(parameters),
_index(getParam<unsigned int>("index")),
_magpie_dat(getMaterialProperty< MAGPIE_DATA >("magpie_data"))
{
}

Real
MAGPIE_Adsorption::computeValue()
{
	MAGPIE_DATA magpie_copy;
	magpie_copy = _magpie_dat[_qp];
	
	//Call MAGPIE Simulation for Unperturbed data
	if (_magpie_dat[_qp].gsta_dat[_index].qmax > 0.0)
	{
		int success = 0;
		//magpie_copy.sys_dat.Output = true;
		success = MAGPIE( (void *)&magpie_copy );
		if (success < 0 || success > 5)
		{
			mError(simulation_fail);
			std::cout << success << std::endl;
			std::cout << "index = " << _index << "\tq = " << magpie_copy.gpast_dat[_index].q << std::endl;
		}
		else success = 0;
		
		//std::cout << "q = " << magpie_copy.gpast_dat[_index].q << std::endl;
		//std::cout << "qT = " << magpie_copy.sys_dat.qT << std::endl;
		//std::cout << "gama = " << magpie_copy.mspd_dat[_index].gama << std::endl;
		//std::cout << "x = " << magpie_copy.gpast_dat[_index].x << std::endl;
		
		return magpie_copy.gpast_dat[_index].q;
		
		//Temporary override to demonstrate LDF kinetics
		//double k = magpie_copy.gpast_dat[_index].q * 0.008;
		//double qe = magpie_copy.gpast_dat[_index].q;
		//return (_u_old[_qp] + (_dt*k*qe))/(1.0+(_dt*k));
	}
	else
	{
		return 0.0;
	}
}
