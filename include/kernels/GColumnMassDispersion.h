/*!
 *  \file GColumnMassDispersion.h
 *	\brief Kernel for use with the corresponding DGColumnMassDispersion object
 *	\details This file creates a standard MOOSE kernel that is to be used in conjunction with the DGColumnMassDispersion kernel
 *			for the discontinous Galerkin formulation of mass dispersion in a fixed-bed adsorber. This kernel is coupled with
 *			material properties, then uses that information to override the Diffusion parameter of the base class and call its
 *			methods.
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

#ifndef GCOLUMNMASSDISPERSION_H
#define GCOLUMNMASSDISPERSION_H

#include "GAnisotropicDiffusion.h"

/// GColumnHeatDispersion class object forward declarations
class GColumnMassDispersion;

template<>
InputParameters validParams<GColumnMassDispersion>();

/// GColumnMassDispersion class object inherits from GAnisotropicDiffusion object
/** This class object inherits from the GAnisotropicDiffusion object in DGOSPREY.
	It must be used in conjunction with the DGColumnMassDispersion object to complete
	the physical description of diffusion for DG methods in MOOSE. The dispersion and
	molecular diffusion material properties are coupled with this object and is used to 
	form/override the diffusion tensor of the base class. Then the base class methods 
	are called to form the residuals and Jacobian elements. */
class GColumnMassDispersion : public GAnisotropicDiffusion
{
public:
	/// Required constructor for objects in MOOSE
	GColumnMassDispersion(const InputParameters & parameters);
	
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
	unsigned int _index;												///< Index of the species of interest for this kernel
	const MaterialProperty<std::vector<Real> > & _dispersion;			///< Reference to the dispersion material property
	const MaterialProperty<std::vector<Real> > & _molecular_diffusion;	///< Reference to the molecular diffusion material property
};

#endif //GCOLUMNMASSDISPERSION_H
