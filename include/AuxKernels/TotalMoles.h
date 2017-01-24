/*!
 *  \file TotalMoles.h
 *	\brief Aux Kernel to calculate total moles of a gas in the bed.
 *	\details This file is responsible for calculating the total moles of a gas in the bed.
 *
 *  \author Austin Ladshaw
 *	\date 01/24/2017
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

#ifndef TOTALMOLES_H
#define TOTALMOLES_H

/// TotalMoles class object forward declaration
class TotalMoles;

template<>
InputParameters validParams<TotalMoles>();

/// TotalMoles class inherits from AuxKernel
/** This class object creates an AuxKernel for use in the MOOSE framework. The AuxKernel will
	calculate the Total moles of a gas in the bed. */
class TotalMoles : public AuxKernel
{
public:
	/// Standard MOOSE public constructor
	TotalMoles(const InputParameters & parameters);
	
protected:
	/// Required MOOSE function override
	/** This is the function that is called by the MOOSE framework when a calculation of the AuxVariable
		is needed. You are required to override this function for any inherited AuxKernel. */
	virtual Real computeValue();
	
	const MaterialProperty<Real> & _porosity;			///< Reference to the bed bulk porosity material property
	const MaterialProperty<Real> & _pellet_density;		///< Reference to the pellet density material property
	const MaterialProperty<Real> & _bed_length;			///< Reference to the bed length
	const MaterialProperty<Real> & _inner_dia;			///< Reference to the bed inner diameter
	const VariableValue & _solid;						///< Pointer to coupled adsorption at the current time
	const VariableValue & _gas;							///< Pointer to coupled adsorption at the previous time
	
private:
	
};

#endif
