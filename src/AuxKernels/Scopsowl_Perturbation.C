/*!
 *  \file Scopsowl_Perturbation.h
 *	\brief Auxillary kernel to calculate perturbed adsorption kinetics of a particular gas species in the system
 *  \author Austin Ladshaw
 *	\date 06/21/2016
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

#include "Scopsowl_Perturbation.h"

template<>
InputParameters validParams<Scopsowl_Perturbation>()
{
	InputParameters params = validParams<AuxKernel>();
	params.addParam<unsigned int>("index",0,"Index of the species that we are interested in.");
	return params;
}

Scopsowl_Perturbation::Scopsowl_Perturbation(const InputParameters & parameters) :
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
Scopsowl_Perturbation::computeValue()
{
	int success = 0;
	Real q = 0.0;
	
	// Initial Conditions
	if (_dt == 0.0)
	{
		_dat[_current_elem->id()] = _owl_dat[_qp];
		
		success = setup_SCOPSOWL_DATA(NULL, default_adsorption, default_retardation, default_pore_diffusion, default_filmMassTransfer, _dat[_current_elem->id()].eval_surfDiff, (void *)&_dat[_current_elem->id()], _dat[_current_elem->id()].gas_dat, &_dat[_current_elem->id()]);
		if (success != 0) {mError(simulation_fail); return -1;}
		
		_dat[_current_elem->id()].total_pressure = _owl_dat[_qp].total_pressure;
		_dat[_current_elem->id()].gas_temperature = _owl_dat[_qp].gas_temperature;
		_dat[_current_elem->id()].gas_velocity = _owl_dat[_qp].gas_velocity;
		
		double pi = _owl_dat[_qp].y[_index] * _owl_dat[_qp].total_pressure;
		double ci = Cstd(pi,_owl_dat[_qp].gas_temperature) + sqrt(DBL_EPSILON);
		double yi = Pstd(ci,_owl_dat[_qp].gas_temperature) / _owl_dat[_qp].total_pressure;
		
		//Set time step
		for (int i=0; i<_owl_dat[_qp].magpie_dat.sys_dat.N; i++)
		{
			_dat[_current_elem->id()].y[i] = _owl_dat[_qp].y[i];
			
			_dat[_current_elem->id()].finch_dat[i].dt = 0.01;
			_dat[_current_elem->id()].finch_dat[i].t = _dat[_current_elem->id()].finch_dat[i].dt + _dat[_current_elem->id()].finch_dat[i].t_old;
			
			if (_dat[_current_elem->id()].SurfDiff == true && _dat[_current_elem->id()].Heterogeneous == true)
			{
				for (int l=0; l<_dat[_current_elem->id()].finch_dat[i].LN; l++)
				{
					_dat[_current_elem->id()].skua_dat[l].finch_dat[i].dt = _dat[_current_elem->id()].finch_dat[i].dt;
					_dat[_current_elem->id()].skua_dat[l].finch_dat[i].t = _dat[_current_elem->id()].finch_dat[i].t;
					_dat[_current_elem->id()].skua_dat[l].t_old = _dat[_current_elem->id()].finch_dat[i].t_old;
					_dat[_current_elem->id()].skua_dat[l].t = _dat[_current_elem->id()].finch_dat[i].t;
				}
			}
		}
		_dat[_current_elem->id()].y[_index] = yi;
		_dat[_current_elem->id()].t_old = _dat[_current_elem->id()].finch_dat[0].t_old;
		_dat[_current_elem->id()].t = _dat[_current_elem->id()].finch_dat[0].t;
		
		//Establish parameters
		if (_owl_dat[_qp].param_dat[_index].Adsorbable == false)
			return 0.0;
		
		//Need to set dat.magpie_dat.sys_dat.qT and dat.param_dat[i].xIC for all
		_dat[_current_elem->id()].magpie_dat.sys_dat.qT = 0.0;
		for (int i=0; i<_dat[_current_elem->id()].magpie_dat.sys_dat.N; i++)
			_dat[_current_elem->id()].magpie_dat.sys_dat.qT += _dat[_current_elem->id()].param_dat[i].qIntegralAvg_old;
		for (int i=0; i<_dat[_current_elem->id()].magpie_dat.sys_dat.N; i++)
		{
			if (_dat[_current_elem->id()].magpie_dat.sys_dat.qT > 0.0)
				_dat[_current_elem->id()].param_dat[i].xIC = _dat[_current_elem->id()].param_dat[i].qIntegralAvg_old/_dat[_current_elem->id()].magpie_dat.sys_dat.qT;
			else
				_dat[_current_elem->id()].param_dat[i].xIC = 0.0;
		}
		
		//Establish ICs, then calculate adsorption
		success = set_SCOPSOWL_ICs(&_dat[_current_elem->id()]);
		if (success != 0) {mError(simulation_fail); return -1;}
		
		//Call Executioner
		success = SCOPSOWL_Executioner(&_dat[_current_elem->id()]);
		if (success != 0) {mError(simulation_fail); return -1;}
		
		q = _dat[_current_elem->id()].param_dat[_index].qIntegralAvg;
		
		//Fix info
		for (int i=0; i<_owl_dat[_qp].magpie_dat.sys_dat.N; i++)
		{
			_dat[_current_elem->id()].y[i] = _owl_dat[_qp].y[i];
			
			_dat[_current_elem->id()].finch_dat[i].dt = _dt;
			_dat[_current_elem->id()].finch_dat[i].t = _dat[_current_elem->id()].finch_dat[i].dt + _dat[_current_elem->id()].finch_dat[i].t_old;
			
			if (_dat[_current_elem->id()].SurfDiff == true && _dat[_current_elem->id()].Heterogeneous == true)
			{
				for (int l=0; l<_dat[_current_elem->id()].finch_dat[i].LN; l++)
				{
					_dat[_current_elem->id()].skua_dat[l].finch_dat[i].dt = _dat[_current_elem->id()].finch_dat[i].dt;
					_dat[_current_elem->id()].skua_dat[l].finch_dat[i].t = _dat[_current_elem->id()].finch_dat[i].t;
					_dat[_current_elem->id()].skua_dat[l].t_old = _dat[_current_elem->id()].finch_dat[i].t_old;
					_dat[_current_elem->id()].skua_dat[l].t = _dat[_current_elem->id()].finch_dat[i].t;
				}
			}
		}
		_dat[_current_elem->id()].t_old = _dat[_current_elem->id()].finch_dat[0].t_old;
		_dat[_current_elem->id()].t = _dat[_current_elem->id()].finch_dat[0].t;

		
	}
	// After Initial Conditions
	else
	{
		if (_owl_dat[_qp].param_dat[_index].Adsorbable == false)
			return 0.0;
		
		//Establish parameters
		_dat[_current_elem->id()].total_pressure = _owl_dat[_qp].total_pressure;
		_dat[_current_elem->id()].gas_temperature = _owl_dat[_qp].gas_temperature;
		_dat[_current_elem->id()].gas_velocity = _owl_dat[_qp].gas_velocity;
		
		double pi = _owl_dat[_qp].y[_index] * _owl_dat[_qp].total_pressure;
		double ci = Cstd(pi,_owl_dat[_qp].gas_temperature) + sqrt(DBL_EPSILON);
		double yi = Pstd(ci,_owl_dat[_qp].gas_temperature) / _owl_dat[_qp].total_pressure;
		
		//Set time step
		for (int i=0; i<_owl_dat[_qp].magpie_dat.sys_dat.N; i++)
		{
			_dat[_current_elem->id()].y[i] = _owl_dat[_qp].y[i];
			
			_dat[_current_elem->id()].finch_dat[i].dt = _dt;
			_dat[_current_elem->id()].finch_dat[i].t = _dat[_current_elem->id()].finch_dat[i].dt + _dat[_current_elem->id()].finch_dat[i].t_old;
			
			if (_dat[_current_elem->id()].SurfDiff == true && _dat[_current_elem->id()].Heterogeneous == true)
			{
				for (int l=0; l<_dat[_current_elem->id()].finch_dat[i].LN; l++)
				{
					_dat[_current_elem->id()].skua_dat[l].finch_dat[i].dt = _dat[_current_elem->id()].finch_dat[i].dt;
					_dat[_current_elem->id()].skua_dat[l].finch_dat[i].t = _dat[_current_elem->id()].finch_dat[i].t;
					_dat[_current_elem->id()].skua_dat[l].t_old = _dat[_current_elem->id()].finch_dat[i].t_old;
					_dat[_current_elem->id()].skua_dat[l].t = _dat[_current_elem->id()].finch_dat[i].t;
				}
			}
		}
		_dat[_current_elem->id()].y[_index] = yi;
		_dat[_current_elem->id()].t_old = _dat[_current_elem->id()].finch_dat[0].t_old;
		_dat[_current_elem->id()].t = _dat[_current_elem->id()].finch_dat[0].t;
		
		//Call Executioner
		success = SCOPSOWL_Executioner(&_dat[_current_elem->id()]);
		if (success != 0) {mError(simulation_fail); return -1;}
		
		q = _dat[_current_elem->id()].param_dat[_index].qIntegralAvg;
		
		//Fix the state
		for (int i=0; i<_owl_dat[_qp].magpie_dat.sys_dat.N; i++)
		{
			_dat[_current_elem->id()].y[i] = _owl_dat[_qp].y[i];
		}
		
		//ReCall Executioner
		//success = SCOPSOWL_Executioner(&_dat[_current_elem->id()]);
		//if (success != 0) {mError(simulation_fail); return -1;}
		
		//Reset for next step
		success = SCOPSOWL_reset(&_dat[_current_elem->id()]);
		if (success != 0) {mError(simulation_fail); return -1;}
	}
	
	return q;
}