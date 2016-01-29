/*!
 *  \file CoupledLDF.h
 *	\brief Advanced kernel for a cross coupled linear driving force mechanism
 *	\details This file creates a standard MOOSE kernel for a coupled linear driving force type of mechanism that
 *			can be added to the non-linear residuals. It contains a boolean argument to determine whether
 *			the driving force is gaining or losing, a coefficient for the rate of the driving force, and
 *			a driving value to where the non-linear coupled variable is heading towards.
 *
 *			This file inherits from LinearDrivingForce.h
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

#include "LinearDrivingForce.h"

#ifndef COUPLEDLDF_H
#define COUPLEDLDF_H

/// CoupledLDF class object forward declarations
class CoupledLDF;

template<>
InputParameters validParams<CoupledLDF>();

/// CoupledLDF class object inherits from LinearDrivingForce object
/** This class object inherits from the LinearDrivingForce object in DGOSPREY.
	All public and protected members of this class are required function overrides.
	The kernel has several protected members including: a boolean for gaining or
	losing mechanisms, a coefficient for the rate or strength of the driving force,
	a driving value to where the coupled non-linear variable is driving toward, and
	the coupled non-linear variable.
 
	Additionally, this object couples the driving value to other non-linear variables
 
	\note To create a specific linear driving force kernel, inherit from this class
	and use other non-linear variables or material properties to change the protected
	member values to reflect the physics for your problem. */
class CoupledLDF : public LinearDrivingForce
{
public:
	/// Required constructor for objects in MOOSE
	CoupledLDF(const InputParameters & parameters);
	
protected:
	/// Required residual function for standard kernels in MOOSE
	/** This function returns a residual contribution for this object.*/
	virtual Real computeQpResidual();
	/// Required Jacobian function for standard kernels in MOOSE
	/** This function returns a Jacobian contribution for this object. The Jacobian being
		computed is the associated diagonal element in the overall Jacobian matrix for the
		system and is used in preconditioning of the linear sub-problem. */
	virtual Real computeQpJacobian();
	
	Real _drive_coef;				///< Coefficient for relationship between coupled variables
	VariableValue & _drive_var;		///< Reference to the coupled non-linear variable that is driving
	
private:
	
};
#endif //COUPLEDLDF_H
