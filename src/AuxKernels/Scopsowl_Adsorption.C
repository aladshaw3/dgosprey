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
_owl_dat(getMaterialProperty< SCOPSOWL_DATA >("owl_data")),
_gas_dat(getMaterialProperty< MIXED_GAS >("gas_data"))
{
	//Forces specific execution behavior of the auxkernel
	_exec_flags.clear();
	_exec_flags.push_back(EXEC_INITIAL);
	_exec_flags.push_back(EXEC_TIMESTEP_END);
	
}

Real
Scopsowl_Adsorption::computeValue()
{
	int success = 0;
	Real q = 0.0;
	
	// Initial Conditions
	if (_dt == 0.0)
	{
		_dat[_current_elem->id()] = _owl_dat[_qp];
		_mixed_dat[_current_elem->id()] = _gas_dat[_qp];
		_dat[_current_elem->id()].magpie_dat = _owl_dat[_qp].magpie_dat;
		
		success = setup_SCOPSOWL_DATA(NULL, default_adsorption, default_retardation, default_pore_diffusion, default_filmMassTransfer, _dat[_current_elem->id()].eval_surfDiff, (void *)&_dat[_current_elem->id()], &_mixed_dat[_current_elem->id()], &_dat[_current_elem->id()]);
		if (success != 0) {mError(simulation_fail); return -1;}
		
		_dat[_current_elem->id()].total_pressure = _owl_dat[_qp].total_pressure;
		_dat[_current_elem->id()].gas_temperature = _owl_dat[_qp].gas_temperature;
		_dat[_current_elem->id()].gas_velocity = _owl_dat[_qp].gas_velocity;
		
		for (int i=0; i<_owl_dat[_qp].magpie_dat.sys_dat.N; i++)
		{
			_dat[_current_elem->id()].y[i] = _owl_dat[_qp].y[i];
		}
		
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
		
		q = _dat[_current_elem->id()].param_dat[_index].qIntegralAvg;
		
	}
	// After Initial Conditions
	else
	{
		if (_owl_dat[_qp].param_dat[_index].Adsorbable == false)
			return 0.0;
		_mixed_dat[_current_elem->id()] = _gas_dat[_qp];
		_dat[_current_elem->id()].magpie_dat = _owl_dat[_qp].magpie_dat;
		
		//Establish parameters
		_dat[_current_elem->id()].total_pressure = _owl_dat[_qp].total_pressure;
		_dat[_current_elem->id()].gas_temperature = _owl_dat[_qp].gas_temperature;
		_dat[_current_elem->id()].gas_velocity = _owl_dat[_qp].gas_velocity;
		
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
		_dat[_current_elem->id()].t_old = _dat[_current_elem->id()].finch_dat[0].t_old;
		_dat[_current_elem->id()].t = _dat[_current_elem->id()].finch_dat[0].t;
		
		//_dat[_current_elem->id()].magpie_dat.sys_dat.Output = true;
		
		//Call Executioner
		success = SCOPSOWL_Executioner(&_dat[_current_elem->id()]);
		if (success != 0) {mError(simulation_fail); return -1;}
		
		q = _dat[_current_elem->id()].param_dat[_index].qIntegralAvg;
		
		//Reset for next step
		success = SCOPSOWL_reset(&_dat[_current_elem->id()]);
	}
	
	return q;
}