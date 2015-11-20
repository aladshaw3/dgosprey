/*!
 *  \file MAGPIE_Perturbation.h
 *	\brief Auxillary kernel to calculate the perturbed adsorption equilibria of a particular gas species in the system
 *	\details This file is responsible for calculating the perturbed adsorption equilibria of a particular species
 *			in the system. The MAGPIE object is stored as a material property whose constants are set
 *			in the corresponding material property file (see MagpieAdsorbateProperties.h). That information
 *			is then used to call the MAGPIE routine to calculate the mixed gas perturbed adsorption for a specific
 *			species of interest.
 *
 *			The perturbation is used to approximate the strength of adsorption via first order finite difference.
 *			That adsorption strength is then loosely coupled to the gaseous species non-linear variable through
 *			a retardation coefficient in the mass transport equations. We use loose coupling to improve the efficiency
 *			of the solutions for this multi-scale mass transfer problem. Full coupling would result in significant 
 *			losses in efficiency, or even complete failure to converge. DO NOT TRY FULL COUPLING!
 *
 *			Unfortunately, the material property system has recently changed in MOOSE, making this operation
 *			much less efficient. Under the new system, all material properties are declared as constants
 *			when outside of their respective material property files. This means that in order for me to
 *			call the MAGPIE subroutine, which edits values in the MAGPIE object, I have to create a copy
 *			of the entire object and have the subroutine act on that copy.
 *
 *	\note	We will only use this kernel to approximate the retardation effect of adsorption IF we are neglecting
 *			the micro-scale kinetics of adsorption/mass transfer into the adsorbent pellets. Kinetics coupling
 *			will be accomplished in another kernel.
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

#include "AuxKernel.h"
#include "flock.h"

#ifndef MAGPIE_Perturbation_H
#define MAGPIE_Perturbation_H

/// Magpie Perturbation class object forward declaration
class MAGPIE_Perturbation;

template<>
InputParameters validParams<MAGPIE_Perturbation>();

/// Magpie Perturbation class inherits from AuxKernel
/** This class object creates an AuxKernel for use in the MOOSE framework. The AuxKernel will
	calculate the perturbed equilibria for a given species in the gas phase based on parameters,
	variables, and constants set in the MAGPIE object. Those values include temperature, pressure,
	concentration, and associated equilibrium energy constants. The return value is the adsorption
	perturbation value in mol/kg. */
class MAGPIE_Perturbation : public AuxKernel
{
public:
	/// Standard MOOSE public constructor
	MAGPIE_Perturbation(const InputParameters & parameters);
	
protected:
	/// Required MOOSE function override
	/** This is the function that is called by the MOOSE framework when a calculation of the AuxVariable
		is needed. You are required to override this function for any inherited AuxKernel. */
	virtual Real computeValue();
	
private:
	unsigned int _index;										///< Index of the gaseous species to calculate equilibria for
	const MaterialProperty< MAGPIE_DATA > & _magpie_dat;		///< Material Property holding the MAGPIE data structure
	
};

#endif
