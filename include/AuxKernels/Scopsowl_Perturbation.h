/*!
 *  \file Scopsowl_Perturbation.h
 *	\brief Auxillary kernel to calculate perturbed adsorption kinetics of a particular gas species in the system
 *	\details This file is responsible for calculating the perturbed adsorption kinetics of a particular species
 *			in the system. The SCOPSOWL object is stored as a material property whose constants are set
 *			in the corresponding material property file (see ScopsowlProperties.h). That information
 *			is then used to call the SCOPSOWL routine to calculate the mixed gas adsorption kinetics for
 *			a specific species of interest.
 *
 *			Unfortunately, the material property system has recently changed in MOOSE, making this operation
 *			much less efficient. Under the new system, all material properties are declared as constants
 *			when outside of their respective material property files. This means that in order for me to
 *			call the SCOPSOWL subroutine, which edits values in the MAGPIE object, I have to create a copy
 *			of the entire object and have the subroutine act on that copy.
 *
 *			In addition, although SCOPSOWL calculates the kinetics of all species in the mixture simultaneously,
 *			MOOSE only allows kernels to act on and calculate properties and values of a single variable at a
 *			time. This means that to get the kinetics of all species in MOOSE, you have to redundantly perform
 *			these calculations again for each other species of interest. This is highly inefficient, but under
 *			time constraints it is the fastest way to finish this project.
 *
 *  \author Austin Ladshaw
 *	\date 06/21/2016
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

#include "AuxKernel.h"
#include "flock.h"

#ifndef Scopsowl_Perturbation_H
#define Scopsowl_Perturbation_H

/// Scopsowl Perturbation class object forward declaration
class Scopsowl_Perturbation;

template<>
InputParameters validParams<Scopsowl_Perturbation>();

/// dataStore function for SCOPSOWL_DATA
/** This function is now REQUIRED by the MOOSE system to utilize the restart and multiapp functionality of
	MOOSE. Although this program does not use these new features, our source code must include these functions
	or the code will NO LONGER COMPILE.
 
	NOTE: Currently, these functions are blank, because we do not use them for anything. However, we can finalize
	these functions by adding in separate dataStore functions for each individual data type in the data structure.
	(See DataIO.h in moose/framework/include/restart/ for the individual functions based on type.)*/
template<>
inline void
dataStore(std::ostream & stream, SCOPSOWL_DATA & data, void * context)
{
	
}

/// dataLoad function for SCOPSOWL_DATA
/** This function is now REQUIRED by the MOOSE system to utilize the restart and multiapp functionality of
	MOOSE. Although this program does not use these new features, our source code must include these functions
	or the code will NO LONGER COMPILE.
 
	NOTE: Currently, these functions are blank, because we do not use them for anything. However, we can finalize
	these functions by adding in separate dataLoad functions for each individual data type in the data structure.
	(See DataIO.h in moose/framework/include/restart/ for the individual functions based on type.)*/
template<>
inline void
dataLoad(std::istream & stream, SCOPSOWL_DATA & data, void * context)
{
	
}

/// dataStore function for MIXED_GAS
/** This function is now REQUIRED by the MOOSE system to utilize the restart and multiapp functionality of
	MOOSE. Although this program does not use these new features, our source code must include these functions
	or the code will NO LONGER COMPILE.
 
	NOTE: Currently, these functions are blank, because we do not use them for anything. However, we can finalize
	these functions by adding in separate dataStore functions for each individual data type in the data structure.
	(See DataIO.h in moose/framework/include/restart/ for the individual functions based on type.)*/
template<>
inline void
dataStore(std::ostream & stream, MIXED_GAS & data, void * context)
{
	
}

/// dataLoad function for MIXED_GAS
/** This function is now REQUIRED by the MOOSE system to utilize the restart and multiapp functionality of
	MOOSE. Although this program does not use these new features, our source code must include these functions
	or the code will NO LONGER COMPILE.
 
	NOTE: Currently, these functions are blank, because we do not use them for anything. However, we can finalize
	these functions by adding in separate dataLoad functions for each individual data type in the data structure.
	(See DataIO.h in moose/framework/include/restart/ for the individual functions based on type.)*/
template<>
inline void
dataLoad(std::istream & stream, MIXED_GAS & data, void * context)
{
	
}

/// Scopsowl Perturbation class inherits from AuxKernel
/** This class object creates an AuxKernel for use in the MOOSE framework. The AuxKernel will
	calculate the perturbed adsorption kinetics for a given species in the gas phase based on parameters,
	variables, and constants set in the SCOPSOWL object. Those values include temperature, pressure,
	concentration, and associated equilibrium energy constants. The return value is the adsorption
	value in mol/kg. */
class Scopsowl_Perturbation : public AuxKernel
{
public:
	/// Standard MOOSE public constructor
	Scopsowl_Perturbation(const InputParameters & parameters);
	
protected:
	/// Required MOOSE function override
	/** This is the function that is called by the MOOSE framework when a calculation of the AuxVariable
		is needed. You are required to override this function for any inherited AuxKernel. */
	virtual Real computeValue();
	
private:
	unsigned int _index;									///< Index of the gaseous species to calculate equilibria for
	Real _dt0;												///< Initial dt for perturbation synchronization (hrs)
	const MaterialProperty< SCOPSOWL_DATA > & _owl_dat;		///< Material Property holding the SCOPSOWL data structure
	const MaterialProperty< MIXED_GAS > & _gas_dat;			///< Material Property holding the MIXED_GAS data structure
	std::map< unsigned int, SCOPSOWL_DATA > _dat;			///< Map for holding material property info for SCOPSOWL
	std::map< unsigned int, MIXED_GAS > _mixed_dat;			///< Map for holding material property info for MIXED_GAS
	
};


#endif /* Scopsowl_Perturbation_h */
