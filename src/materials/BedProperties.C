/*!
 *  \file AdsorbentProperties.h
 *	\brief Material Properties kernel that will setup and hold all information associated with the fixed-bed
 *	\warning THIS KERNEL IS INCOMPLETE! ONLY USED FOR DATA STORAGE FOR VARIOUS BED PARAMETERS!
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

#include "BedProperties.h"

template<>
// input parameters are the parameters that are constant and not calculated from other parameters
InputParameters validParams<BedProperties>()
{
  InputParameters params = validParams<Material>();
  
  params.addParam<Real>("length", "Bed length (cm)");
  params.addParam<Real>("din", "Wall inner diameter (cm)");
  params.addParam<Real>("dout","Wall outer diameter (cm)");
  params.addParam<Real>("eb","Bed bulk porosity");
  params.addParam<Real>("Kz","Effective Axial Thermal Conductivity (J/hr/cm/K");
  params.addParam<Real>("rhow","Wall material density (g/cm^3)");
  params.addParam<Real>("hw","Wall heat capacity (J/g/K)");
  params.addParam<Real>("Uw","Bed-Wall Heat Transfer Coefficient (J/hr/cm^2/K)");
  params.addParam<Real>("Ua","Wall-External Heat Transfer Coefficient (J/hr/cm^2/K)");
  params.addCoupledVar("temperature","Coupled variable for temperature");
  params.addCoupledVar("coupled_gases", "Gas concentrations variables being coupled");
  
return params;
}

BedProperties::BedProperties(const InputParameters & parameters)
  :Material(parameters),
   
   _length(getParam<Real>("length")),
   _din(getParam<Real>("din")),
   _dout(getParam<Real>("dout")),
   _eb(getParam<Real>("eb")),
   _Kz(getParam<Real>("Kz")),
   _rhow(getParam<Real>("rhow")),
   _hw(getParam<Real>("hw")),
   _Uw(getParam<Real>("Uw")),
   _Ua(getParam<Real>("Ua")),
   _inner_dia(declareProperty<Real>("inner_dia")),
   _outer_dia(declareProperty<Real>("outer_dia")),
   _porosity(declareProperty<Real>("porosity")),
   _conductivity(declareProperty<Real>("conductivity")),
   _wall_density(declareProperty<Real>("wall_density")),
   _wall_heat_capacity(declareProperty<Real>("wall_heat_capacity")),
   _bed_wall_transfer_coeff(declareProperty<Real>("bed_wall_transfer_coeff")),
   _wall_exterior_transfer_coeff(declareProperty<Real>("wall_exterior_transfer_coeff")),
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
BedProperties::computeQpProperties()
{
  //For constant bed properties...
  _inner_dia[_qp] = _din;
  _outer_dia[_qp] = _dout;
  _porosity[_qp] = _eb;
  _conductivity[_qp] = _Kz;
  _wall_density[_qp] = _rhow;
  _wall_heat_capacity[_qp] = _hw;
  _bed_wall_transfer_coeff[_qp] = _Uw;
  _wall_exterior_transfer_coeff[_qp] = _Ua;
  //Note: some of these may vary with space or temperature or concentration, but for now we assume constant 
}
