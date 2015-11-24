/*!
 *  \file BedHeatAccumulation.h
 *	\brief Time Derivative kernel for the accumulation of heat in a fixed-bed column
 *	\details This file creates a time derivative kernel to be used in the energy transport equations for adsorption
 *			in a fixed-bed column. It combines the retardation coefficient from a material property with the standard
 *			time derivative kernel object in MOOSE.
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

#ifndef BEDHEATACCUMULATION_H
#define BEDHEATACCUMULATION_H

#include "TimeDerivative.h"

/// BedHeatAccumulation class object forward declarations
class BedHeatAccumulation;

template<>
InputParameters validParams<BedHeatAccumulation>();

/// BedHeatAccumulation class object inherits from TimeDerivative object
/** This class object inherits from the TimeDerivative object.
	All public and protected members of this class are required function overrides.
	The flux BC uses the velocity and diffusivity in the system to apply a boundary
	condition based on whether or not material is leaving or entering the boundary. */
class BedHeatAccumulation : public TimeDerivative
{
public:
  
  BedHeatAccumulation(const InputParameters & parameters);
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  
private:
	const MaterialProperty<Real> & _heat_retardation;
};

#endif // BEDHEATACCUMULATION_H
