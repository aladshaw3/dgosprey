/*!
 *  \file FlowProperties.h
 *	\brief Material Properties kernel that will setup and calculate gas flow properties based on physical characteristics
 *	\details This file creates a material property object for the flow properties in the fixed-bed column. The flow properties
 *			are calculated based on some dimensional analysis, empirical relationships, and kinetic theory of gases. Those
 *			properties are then coupled the with mass and energy kernels in the domain to simulate the dynamic and non-linear
 *			behavior in the system.
 *
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

#include "Material.h"
#include "flock.h"

#ifndef FLOWPROPERTIES_H
#define FLOWPROPERTIES_H

#ifndef _gas_const
#define _gas_const 8.3144621		///< Gas Law Constant - J/K/mol
#endif

/// FlowProperties class object forward declaration
class FlowProperties;

template<>
InputParameters validParams<FlowProperties>();

/// FlowProperties class object inherits from Material object
/** This class object inherits from the Material object in the MOOSE framework.
	All public and protected members of this class are required function overrides. The object
	will set up and calculate various flow properties including linear velocities, molecular diffusion,
	mechanical dispersion, gas density, gas viscosity, gas heat capacity, etc. This object also approximates
	the effective retardation coefficient for each species in the mass balance. The evaluation of that 
	parameter is dependent on the solid phase concentration variable, which will either be calculated by
	MAGPIE (see magpie.h) or SCOPSOWL (see scopsowl.h) dependending on whether or not we will consider
	adsorption kinetics in the simulation. */
class FlowProperties : public Material
{
public:
	/// Required constructor for objects in MOOSE
	FlowProperties(const InputParameters & parameters);
	
protected:
	/// Required function override for Material objects in MOOSE
	/** This function computes the material properties when they are needed by other MOOSE objects.*/
	virtual void computeQpProperties();
	
	virtual void initQpStatefulProperties();
	
private:
	std::vector<Real> _molecular_weight;							///< Molecular weights for each gas species (g/mol)
	std::vector<Real> _comp_heat_capacity;							///< Heat capacities for each gas species (J/g/K)
	std::vector<Real> _comp_ref_viscosity;							///< Sutherland's reference viscosity for each gas species (g/cm/s)
	std::vector<Real> _comp_ref_temp;								///< Sutherland's reference temperature for each species (K)
	std::vector<Real> _comp_Sutherland_const;						///< Sutherland's constant for each gas species (K)
	Real _flow_rate;												///< Inlet flow rate for the fixed-bed column (cm^3/hr)
	Real _column_length;											///< Length of the fixed-bed column (cm)
  
	MaterialProperty<Real> & _velocity;								///< MaterialProperty for the linear velocity in the bed (cm/hr)
	MaterialProperty<Real> & _gas_density;							///< MaterialProperty for the gas density (g/cm^3)
	MaterialProperty<Real> & _gas_viscosity;						///< MaterialProperty for the gas viscosity (g/cm/s)
	MaterialProperty<Real> & _gas_heat_capacity;					///< MaterialProperty for the gas heat capacity (J/g/K)
	MaterialProperty<Real> & _gas_molecular_wieght;					///< MaterialProperty for the gas total molecular wieght (g/mol)
	
	const MaterialProperty<Real> & _inner_dia;						///< Coupled material property for bed inner diameter
	const MaterialProperty<Real> & _porosity;						///< Coupled material property for bed bulk porosity
	const MaterialProperty<Real> & _pellet_density;					///< Coupled material property for adsorbent pellet density
	const MaterialProperty<Real> & _pellet_heat_capacity;			///< Coupled material property for adsorbent heat capacity
	const MaterialProperty<Real> & _pellet_diameter;				///< Coupled material property for the adsorbent pellet diameter
	const MaterialProperty<Real> & _binder_porosity;				///< MaterialProperty for the binder porosity
	const MaterialProperty<Real> & _pore_size;						///< MaterialProperty for the macropore radius (cm)
	
	MaterialProperty<Real> & _heat_retardation;						///< MaterialProperty for energy balance retardation coefficient
	MaterialProperty<std::vector<Real> > & _molecular_diffusion;	///< MaterialProperty for each species' molecular diffusion (cm^2/s)
	MaterialProperty<std::vector<Real> > & _dispersion;			///< MaterialProperty for each species' dispersion coefficient (cm^2/hr)
	MaterialProperty<std::vector<Real> > & _retardation;		///< MaterialProperty for each species' retardation coefficient
	MaterialProperty< MIXED_GAS > & _mixed_gas;					///< MaterialProperty for the MIXED_GAS struct in egret.h
	MaterialProperty< MIXED_GAS > & _mixed_gas_old;					///< Old MaterialProperty for the MIXED_GAS struct in egret.h
	
	MaterialProperty<std::vector<Real> > & _film_transfer;			///< MaterialProperty for the film mass transfer coeff (cm/hr)
	MaterialProperty<std::vector<Real> > & _pore_diffusion;			///< MaterialProperty for the pore diffusion (cm^2/hr)
  	
  	VariableValue & _temperature;					///< Reference to the coupled column temperature
	VariableValue & _total_pressure;				///< Reference to the coupled column pressure
	std::vector<unsigned int> _index;				///< Indices for the gas species in the system
	std::vector<VariableValue *> _gas_conc;			///< Pointer list to the coupled gases
	std::vector<VariableValue *> _solid_conc;		///< Pointer list to the coupled adsorption concentrations
	std::vector<VariableValue *> _solid_perturb;	///< Pointer list to the coupled adsorption perturbations
	
};


#endif //FLOWPROPERTIES_H
