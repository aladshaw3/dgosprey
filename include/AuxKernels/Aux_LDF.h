/*!
 *  \file Aux_LDF.h
 *	\brief Generic auxillary kernel to calculate the value of an aux variable using LDF kinetics
 *	\details This file is responsible for calculating the value of the aux variable based on an implicit
 *			integration of the linear driving force expression. It's intended use will be to create a 
 *			generic class that can be inherited by a more specific class to have certain values overriden
 *			that may be coulped to other non-linear variables in the simuation. Coupling between this aux
 *			kernel and other non-linear variables should be done "loosely" as the intent will ultimately
 *			be to couple multi-scale physical phenomena. DO NOT try to fully couple this with non-linear
 *			variables. The convergence of the overall system may suffer.
 *
 *  \author Austin Ladshaw
 *	\date 02/04/2016
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

#include "AuxKernel.h"

#ifndef AUX_LDF_H
#define AUX_LDF_H

/// Aux_LDF class object forward declaration
class Aux_LDF;

template<>
InputParameters validParams<Aux_LDF>();

/// Aux_LDF class inherits from AuxKernel
/** This class object creates an AuxKernel for use in the MOOSE framework. The AuxKernel will
	calculate the result of the linear driving force function, integrated implicitly, for the
	aux variable it is associated with. It contains two parameters: (i) the ldf coefficient and
	(ii) the driving value. Inherit from this base class to alter the parameters and change the
	behavior of this kernel to fit your particular problem. */
class Aux_LDF : public AuxKernel
{
public:
	/// Standard MOOSE public constructor
	Aux_LDF(const InputParameters & parameters);
	
protected:
	/// Required MOOSE function override
	/** This is the function that is called by the MOOSE framework when a calculation of the AuxVariable
		is needed. You are required to override this function for any inherited AuxKernel. */
	virtual Real computeValue();
	
	Real _ldf_coef;				/// Value of the driving force coefficient
	Real _driving_value;		/// Value that the driving force is pushing the aux variable towards
	
private:
	
};

#endif
