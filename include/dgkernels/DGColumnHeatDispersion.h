/*!
 *  \file DGColumnHeatDispersion.h
 *	\brief Discontinous Galerkin kernel for energy dispersion in a fixed-bed column
 *	\details This file creates a discontinous Galerkin kernel for the dispersive heat transfer in a fixed-bed column.
 *			The dispersion portion of the energy transport equations involves the thermal conductivity in the system.
 *			That parameter is calculated in a material properties file and passed into this object for use the construction
 *			of residuals and Jacobians.
 *
 *	\note Any DG kernel under DGOSPREY will have a cooresponding G kernel (usually of same name) that must be included
 *		with the DG kernel in the input file. This is because the DG finite element method breaks into several different
 *		residual pieces, only a handful of which are handled by the DG kernel system and the other parts must be handled
 *		by the standard Galerkin system. This my be due to some legacy code in MOOSE. I am not sure if it is possible to
 *		lump all of these actions into a single DG kernel.
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

#include "DGAnisotropicDiffusion.h"

#ifndef DGCOLUMNHEATDISPERSION_H
#define DGCOLUMNHEATDISPERSION_H

/// DGColumnHeatDispersion class object forward declarations
class DGColumnHeatDispersion;

template<>
InputParameters validParams<DGColumnHeatDispersion>();

/// DGColumnHeatDispersion class object inherits from DGAnisotropicDiffusion object
/** This class object inherits from the DGAnisotropicDiffusion object in OSPREY.
	All public and protected members of this class are required function overrides. The object
	will provide residuals and Jacobians for the discontinous Galerkin formulation of the heat
	dispersion physics in a fixed-bed column. Parameters for this kernel are given as material
	properties and will be used to override the inherited classes diffusion tensor.
 
	\note As a reminder, any DGKernel in MOOSE was be accompanied by the equivalent GKernel in
	order to provide the full residuals and Jacobians for the system. */
class DGColumnHeatDispersion : public DGAnisotropicDiffusion
{
public:
	/// Required constructor for objects in MOOSE
	DGColumnHeatDispersion(const InputParameters & parameters);
	
protected:
	/// Required residual function for DG kernels in MOOSE
	/** This function returns a residual contribution for this object.*/
	virtual Real computeQpResidual(Moose::DGResidualType type);
	/// Required Jacobian function for DG kernels in MOOSE
	/** This function returns a Jacobian contribution for this object. The Jacobian being
		computed is the associated diagonal element in the overall Jacobian matrix for the
		system and is used in preconditioning of the linear sub-problem. */
	virtual Real computeQpJacobian(Moose::DGJacobianType type);
	
private:
	const MaterialProperty<Real> & _conductivity;		///< Reference to the thermal conductivity material property
};

#endif
