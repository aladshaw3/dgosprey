/*!
 *  \file MAGPIE_MaterialLDF_Perturbation.h
 *	\brief Auxillary kernel to calculate adsorption perturbation based on LDF kinetics with material property coefficients
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

#include "MAGPIE_MaterialLDF_Perturbation.h"

template<>
InputParameters validParams<MAGPIE_MaterialLDF_Perturbation>()
{
	InputParameters params = validParams<Aux_LDF>();
	params.addParam<unsigned int>("index",0,"Index of the species that we are interested in.");
	return params;
}

MAGPIE_MaterialLDF_Perturbation::MAGPIE_MaterialLDF_Perturbation(const InputParameters & parameters) :
Aux_LDF(parameters),
_index(getParam<unsigned int>("index")),
_magpie_dat(getMaterialProperty< MAGPIE_DATA >("magpie_data")),
_pellet_diameter(getMaterialProperty<Real>("pellet_diameter")),
_porosity(getMaterialProperty<Real>("porosity")),
_binder_porosity(getMaterialProperty<Real>("binder_porosity")),
_crystal_radius(getMaterialProperty<Real>("crystal_radius")),
_pellet_density(getMaterialProperty<Real>("pellet_density")),
_film_transfer(getMaterialProperty<std::vector<Real> >("film_transfer")),
_pore_diffusion(getMaterialProperty<std::vector<Real> >("pore_diffusion")),
_surface_diffusion(getMaterialProperty<std::vector<Real> >("surface_diffusion"))
{
	//Forces specific execution behavior of the auxkernel
	_exec_flags.clear();
	_exec_flags.push_back(EXEC_INITIAL);
	_exec_flags.push_back(EXEC_TIMESTEP_END);
}

Real MAGPIE_MaterialLDF_Perturbation::computeValue()
{
	MAGPIE_DATA magpie_copy;
	magpie_copy = _magpie_dat[_qp];
	
	//Call MAGPIE Simulation for Unperturbed data
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
		
		_driving_value = magpie_copy.gpast_dat[_index].q;
		
		if (std::isnan(_driving_value))
			_driving_value = 0.0;
		
	}
	else
	{
		_driving_value = 0.0;
	}
	
	//Calculate the partition ratio
	Real _part_rat = (1.0-_porosity[_qp])*_pellet_density[_qp] * q_p((magpie_copy.gpast_dat[_index].y * magpie_copy.sys_dat.PT), (void *)&magpie_copy, _index) * Rstd * magpie_copy.sys_dat.T;
	
	//Calculate each resistance
	Real _kf_res, _Dp_res, _Ds_res;
	if (_magpie_dat[_qp].gsta_dat[_index].qmax > 0.0)
	{
		_kf_res = (6.0*(1.0-_porosity[_qp])*_film_transfer[_qp][_index])/(_part_rat*_pellet_diameter[_qp])/25.0;
		_Dp_res = (60.0*(1.0-_porosity[_qp])*_binder_porosity[_qp]*_pore_diffusion[_qp][_index])/(_part_rat*_pellet_diameter[_qp]*_pellet_diameter[_qp])/25.0;
		_Ds_res = (15.0*_surface_diffusion[_qp][_index])/(_crystal_radius[_qp]*_crystal_radius[_qp])/25.0;
	}
	else
	{
		_ldf_coef = 0.0;
		return Aux_LDF::computeValue();
	}
	
	//Check for errors
	if (std::isnan(_kf_res) || std::isinf(_kf_res))
		_kf_res = sqrt(DBL_MAX);
	if (std::isnan(_Dp_res))
		_Dp_res = sqrt(DBL_MAX);
	if (std::isnan(_Ds_res))
		_Ds_res = sqrt(DBL_MAX);
	
	
	//Calculate the LDF coefficient
	if (_surface_diffusion[_qp][_index] == 0.0)
	{
		_ldf_coef = (_kf_res + (1.0/((1.0/_Dp_res))));
	}
	else
	{
		_ldf_coef = (_kf_res + (1.0/((1.0/_Dp_res)+(1.0/_Ds_res))));
	}
	
	return Aux_LDF::computeValue();
}