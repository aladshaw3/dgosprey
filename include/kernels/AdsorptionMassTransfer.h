/*!
 *  \file AdsorptionMassTransfer.h
 *	\brief Standard kernel for the transfer of mass via adsorption
 *	\details This file creates a standard MOOSE kernel for the transfer of mass between the
 *			bulk gas of the fixed-bed and the adsorbent material in the column. The
 *			mass transfer is based on the amount of material in the bed and the solid adsorption variables.
 *
 *  \author Austin Ladshaw
 *	\date 01/29/2016
 *	\copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science and was developed for use
 *			   by Idaho National Laboratory and Oak Ridge National Laboratory
 *			   engineers and scientists. Portions Copyright (c) 2016, all
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

#ifndef AdsorptionMassTransfer_H
#define AdsorptionMassTransfer_H

#include "Kernel.h"

/// AdsorptionHeatAccumulation class object forward declarationss
class AdsorptionMassTransfer;

template<>
InputParameters validParams<AdsorptionMassTransfer>();

/// AdsorptionMassTransfer class object inherits from Kernel object
/** This class object inherits from the Kernel object in the MOOSE framework.
	All public and protected members of this class are required function overrides.
	The kernel interfaces the material properties for the bulk bed porosity and the
	pellet density, as well as coupling with adsorption as it changes in time,
	in order to form a residuals and Jacobians for the gas concentration variable. */
class AdsorptionMassTransfer : public Kernel
{
public:
	/// Required constructor for objects in MOOSE
	AdsorptionMassTransfer(const InputParameters & parameters);
	
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
	VariableValue & _solid;								///< Pointer to coupled adsorption at the current time
	VariableValue & _solid_old;							///< Pointer to coupled adsorption at the previous time
};

#endif