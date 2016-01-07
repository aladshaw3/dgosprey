/*!
 *  \file AdsorbentProperties.h
 *	\brief Material Properties kernel that will setup and hold all information associated with the adsorbent
 *	\warning THIS KERNEL IS INCOMPLETE! ONLY USED FOR DATA STORAGE FOR PELLET DENSITY AND HEAT CAPACITY!
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

#include "AdsorbentProperties.h"

template<>
// input parameters are the parameters that are constant and not calculated from other parameters
InputParameters validParams<AdsorbentProperties>()
{
  InputParameters params = validParams<Material>();
  
  params.addParam<Real>("binder_fraction","Binder fraction of the adsorbent pellet");
  params.addParam<Real>("binder_porosity","Porosity of the binder material in the adsorbent pellet");
  params.addParam<Real>("crystal_radius","Radius of the adsorbent crystals in the binder matrix (um)");
  params.addParam<Real>("pellet_diameter","Diameter of the adsorbent pellet (cm)");
  params.addParam<Real>("macropore_radius","Nominal pore size of the macropores in the binder material (cm)");
  params.addParam<Real>("pellet_density","Density of the adsorbent pellet (g/cm^3)");
  params.addParam<Real>("pellet_heat_capacity","Pellet heat capacity (J/g/K)");
  params.addCoupledVar("temperature","Coupled variable for temperature");
  params.addCoupledVar("coupled_gases", "Gas concentrations variables being coupled");
  
  return params;
}

AdsorbentProperties::AdsorbentProperties(const InputParameters & parameters)
:Material(parameters),

_binder_fraction(getParam<Real>("binder_fraction")),
_binder_porosity(getParam<Real>("binder_porosity")),
_crystal_radius(getParam<Real>("crystal_radius")),
_pellet_diameter(getParam<Real>("pellet_diameter")),
_macropore_radius(getParam<Real>("macropore_radius")),
_rhos(getParam<Real>("pellet_density")),
_hs(getParam<Real>("pellet_heat_capacity")),
_pellet_density(declareProperty<Real>("pellet_density")),
_pellet_heat_capacity(declareProperty<Real>("pellet_heat_capacity")),
_temperature(coupledValue("temperature"))

{
	unsigned int n = coupledComponents("coupled_gases");
	_index.resize(n);
	_gas_conc.resize(n);
	
	for (unsigned int i = 0; i<_gas_conc.size(); ++i)
	{
		_index[i] = coupled("coupled_gases",i); //may only be useful for compute Jacobian Off Diagonal (see ~/projects/truck/moose/modules/chemical_reactions/ .../ CoupledConvectionReactionSub.C)
		_gas_conc[i] = &coupledValue("coupled_gases",i);
	}
	/*
	 Note: When using _gas_conc[i], it must be appropriately dereferenced as follows...
	 (*_gas_conc[i])[_qp] = concentration of species i at node _qp
	 */
}

void
AdsorbentProperties::computeQpProperties()
{
  //For constant bed properties...
  _pellet_density[_qp] = _rhos;
  _pellet_heat_capacity[_qp] = _hs;
  //Note: some of these may vary with space or temperature or concentration, but for now we assume constant
}
