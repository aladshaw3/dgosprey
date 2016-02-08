/*!
 *  \file MAGPIE_MaterialLDF_Adsorption.h
 *	\brief Auxillary kernel to calculate adsorption based on LDF kinetics with material property coefficients
 *  \author Austin Ladshaw
 *	\date 02/05/2016
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

#include "MAGPIE_MaterialLDF_Adsorption.h"

template<>
InputParameters validParams<MAGPIE_MaterialLDF_Adsorption>()
{
	InputParameters params = validParams<Aux_LDF>();
	params.addParam<unsigned int>("index",0,"Index of the species that we are interested in.");
	return params;
}

MAGPIE_MaterialLDF_Adsorption::MAGPIE_MaterialLDF_Adsorption(const InputParameters & parameters) :
Aux_LDF(parameters),
_index(getParam<unsigned int>("index")),
_magpie_dat(getMaterialProperty< MAGPIE_DATA >("magpie_data")),
_pellet_diameter(getMaterialProperty<Real>("pellet_diameter")),
_porosity(getMaterialProperty<Real>("porosity")),
_binder_porosity(getMaterialProperty<Real>("binder_porosity")),
_crystal_radius(getMaterialProperty<Real>("crystal_radius")),
_partition_ratio(getMaterialProperty<std::vector<Real> >("partition_ratio")),
_film_transfer(getMaterialProperty<std::vector<Real> >("film_transfer")),
_pore_diffusion(getMaterialProperty<std::vector<Real> >("pore_diffusion")),
_surface_diffusion(getMaterialProperty<std::vector<Real> >("surface_diffusion"))
{
}

Real MAGPIE_MaterialLDF_Adsorption::computeValue()
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
	
	Real _part_rat = 1.69 * 10.16 / 6.66E-5;
	Real _surf_diff = 0.17468;
	
	//Calculate the LDF coefficient
	_ldf_coef = ((_part_rat*_pellet_diameter[_qp])/(6.0*(1.0-_porosity[_qp])*_film_transfer[_qp][_index])) + ((_part_rat*_pellet_diameter[_qp]*_pellet_diameter[_qp])/(60.0*(1.0-_porosity[_qp])*_binder_porosity[_qp]*_pore_diffusion[_qp][_index]));
	
	if (_surface_diffusion[_qp][_index] > 0.0)
		_ldf_coef = _ldf_coef + ((_crystal_radius[_qp]*_crystal_radius[_qp])/(15.0*_surf_diff));
	
	if (_ldf_coef <= 0.0)
		_ldf_coef = 0.0;
	else
		_ldf_coef = 1.0/_ldf_coef;
	
	std::cout << _index << std::endl;
	std::cout <<  "ldf coef = " << _ldf_coef << std::endl;
	
	std::cout << "part coef = " << _partition_ratio[_qp][_index] << std::endl;
	std::cout << "part coef 2 = " << _part_rat << std::endl;
	std::cout << "pellet dia = " << _pellet_diameter[_qp] << std::endl;
	std::cout << "porosity = " << _porosity[_qp] << std::endl;
	std::cout << "film mt = " << _film_transfer[_qp][_index] << std::endl;
	std::cout << "binder pore = " << _binder_porosity[_qp] << std::endl;
	std::cout << "pore diff = " << _pore_diffusion[_qp][_index] << std::endl;
	std::cout << "cry rad = " << _crystal_radius[_qp] << std::endl;
	std::cout << "surf diff = " << _surface_diffusion[_qp][_index] << std::endl;
	std::cout << "surf diff 2 = " << _surf_diff << std::endl;
	
	std::cout << std::endl;
	
	return Aux_LDF::computeValue();
}