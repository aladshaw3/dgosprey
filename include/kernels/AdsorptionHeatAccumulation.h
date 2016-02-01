/*!
 *  \file AdsorptionHeatAccumulation.h
 *	\brief Standard kernel for the heat of adsorption and its effect on the system temperature
 *	\details This file creates a standard MOOSE kernel for the transfer of energy as heat between the
 *			bulk gas temperature of the fixed-bed and the adsorbent material in the column. The
 *			heat transfer is based on the heat of adsorption and the amount currently adsorbed.
 *			It is coupled to the heat of the gas in the column.
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

#ifndef AdsorptionHeatAccumulation_H
#define AdsorptionHeatAccumulation_H

#include "Kernel.h"

/// AdsorptionHeatAccumulation class object forward declarationss
class AdsorptionHeatAccumulation;

template<>
InputParameters validParams<AdsorptionHeatAccumulation>();

/// AdsorptionHeatAccumulation class object inherits from Kernel object
/** This class object inherits from the Kernel object in the MOOSE framework.
	All public and protected members of this class are required function overrides.
	The kernel interfaces the material properties for the bulk bed porosity and the
	pellet density, as well as coupling with the heat of adsorption as it changes in time,
	in order to form a residuals and Jacobians for the gas temperature variable. */
class AdsorptionHeatAccumulation : public Kernel
{
public:
	/// Required constructor for objects in MOOSE
	AdsorptionHeatAccumulation(const InputParameters & parameters);
	
protected:
	/// Required residual function for standard kernels in MOOSE
	/** This function returns a residual contribution for this object.*/
	virtual Real computeQpResidual();
	/// Required Jacobian function for standard kernels in MOOSE
	/** This function returns a Jacobian contribution for this object. The Jacobian being
		computed is the associated diagonal element in the overall Jacobian matrix for the
		system and is used in preconditioning of the linear sub-problem. */
	virtual Real computeQpJacobian();
	
private:
	const MaterialProperty<Real> & _porosity;			///< Reference to the bed bulk porosity material property
	const MaterialProperty<Real> & _pellet_density;		///< Reference to the pellet density material property
	std::vector<VariableValue *> _solid_heat;			///< Pointer list to the coupled heats of adsorption at the current time
	std::vector<VariableValue *> _solid_heat_old;		///< Pointer list to the coupled heats of adsorption at the previous time
};

#endif