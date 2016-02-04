/*!
 *  \file MAGPIE_ConstLDF_Adsorption.h
 *	\brief Auxillary kernel to calculate adsorption based on LDF kinetics with constant coefficients
 *	\details This file is responsible for calculating the adsorption based on linear driving force kinetics
 *			implicitly for the aux variable object. That calculation is based on assuming a constant ldf 
 *			coefficient, but updates the driving value at every iteration based on a MAGPIE simulation
 *			that estimates the new equilibrium point for that aux variable. Remember, it is intended that
 *			this kernel will be loosely coupled to the non-linear variables. Otherwise, this gives poor
 *			performance and may not converge. We use loose coupling because of the multiscale nature of
 *			the physics; we are coupling macro-scale transport to micro-scale equilibria and kinetics.
 *
 *			Unfortunately, the material property system has recently changed in MOOSE, making this operation
 *			much less efficient. Under the new system, all material properties are declared as constants
 *			when outside of their respective material property files. This means that in order for me to
 *			call the MAGPIE subroutine, which edits values in the MAGPIE object, I have to create a copy
 *			of the entire object and have the subroutine act on that copy.
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

#include "Aux_LDF.h"
#include "flock.h"

#ifndef MAGPIE_ConstLDF_Adsorption_H
#define MAGPIE_ConstLDF_Adsorption_H

/// MAGPIE_ConstLDF class object forward declaration
class MAGPIE_ConstLDF_Adsorption;

template<>
InputParameters validParams<MAGPIE_ConstLDF_Adsorption>();

/// MAGPIE_ConstLDF class inherits from AuxKernel
/** This class object inherits from Aux_LDF to calculated the adsorption of an aux variable
	based on a constant linear driving force parameter and a MAGPIE simulation. The MAGPIE
	simulation is used to override the driving value of the base class at every iteration, thus
	coupling the kinetics to the transport problem. NOTE: This coupling should be done loosely
	to avoid poor convergence behavior between the multiple scales of the problem. */
class MAGPIE_ConstLDF_Adsorption : public Aux_LDF
{
public:
	/// Standard MOOSE public constructor
	MAGPIE_ConstLDF_Adsorption(const InputParameters & parameters);
	
protected:
	/// Required MOOSE function override
	/** This is the function that is called by the MOOSE framework when a calculation of the AuxVariable
		is needed. You are required to override this function for any inherited AuxKernel. */
	virtual Real computeValue();
	
private:
	unsigned int _index;									///< Index of the gaseous species to calculate equilibria for
	const MaterialProperty< MAGPIE_DATA > & _magpie_dat;	///< Material Property holding the MAGPIE data structure
	
};

#endif