/*!
 *  \file ParameterizedAdsorptionEquil.h
 *	\brief Standard kernel for  simulating multispecies adsorption equilibria
 *	\details This file creates a standard MOOSE kernel for the call the MAGPIE adsorption equilibria function
 *
 *
 *  \author Austin Ladshaw, Alexander Wiechert
 *	\date 06/12/2017
 *	\copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Alexander Wiechert for PhD research in the area
 *             of adsorption and surface science and was developed for use
 *			   by Idaho National Laboratory and Oak Ridge National Laboratory
 *			   engineers and scientists. Portions Copyright (c) 2017, all
 *             rights reserved.
 *
 *			   Alexander Wiechert does not claim any ownership or copyright to the
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

#include "CoupledLangmuirForcingFunction.h"
#include "flock.h"
#include "DataStruct_StoreLoad.h"

#ifndef ParameterizedAdsorptionEquil_h
#define ParameterizedAdsorptionEquil_h

/// ParameterizedAdsorptionEquil class object forward declarationss
class ParameterizedAdsorptionEquil;

template<>
InputParameters validParams<ParameterizedAdsorptionEquil>();

/// ParameterizedAdsorptionEquil class object inherits from Kernel object
/** This class object inherits from the Kernel object in the MOOSE framework.
	All public and protected members of this class are required function overrides.
	The kernel interfaces the two non-linear variables to couple a Langmuir forcing
	function between given objects. */
class ParameterizedAdsorptionEquil : public CoupledLangmuirForcingFunction
{
public:
    /// Required constructor for objects in MOOSE
    ParameterizedAdsorptionEquil(const InputParameters & parameters);
    
protected:
    /// Required residual function for standard kernels in MOOSE
    /** This function returns a residual contribution for this object.*/
    virtual Real computeQpResidual();
    
    /// Required Jacobian function for standard kernels in MOOSE
    /** This function returns a Jacobian contribution for this object. The Jacobian being
     computed is the associated diagonal element in the overall Jacobian matrix for the
     system and is used in preconditioning of the linear sub-problem. */
    virtual Real computeQpJacobian();
    
    /// Not Required, but aids in the preconditioning step
    /** This function returns the off diagonal Jacobian contribution for this object. By
     returning a non-zero value we will hopefully improve the convergence rate for the
     cross coupling of the variables. */
    virtual Real computeQpOffDiagJacobian(unsigned int jvar);
    
private:
    unsigned int _index;									///< Index of the gaseous species to calculate equilibria for
    const MaterialProperty< MAGPIE_DATA > & _magpie_dat;	///< Material Property holding the MAGPIE data structure
	const VariableValue & _q_est;							///< Coupled adsorption variable 
  
};



#endif /* CoupledLangmuirForcingFunction_h */
