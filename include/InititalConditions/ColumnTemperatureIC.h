/*!
 *  \file ColumnTemperatureIC.h
 *	\brief Initial Condition kernel for initial temperature in a fixed-bed column
 *	\details This file creates an initial condition for the temperature in the bed. The initial condition for temperature
 *			is assumed a constant value at all points in the bed. However, this can be modified later to include spatially
 *			varying initial conditions for temperature.
 *
 *	\note If you want to have spatially varying initial conditions, you will need to modify the virtual value function of
 *		this kernel. Otherwise, it is assumed that the non-linear variable is initially constant at all points in the domain.
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

#ifndef COLUMNTEMPERATUREIC_H
#define	COLUMNTEMPERATUREIC_H

#include "InitialCondition.h"

/// ColumnTemperatureIC class object forward declarations
class ColumnTemperatureIC;

template<> InputParameters validParams<ColumnTemperatureIC>();

class ColumnTemperatureIC : public InitialCondition
{
public:
	ColumnTemperatureIC(const InputParameters & parameters);
	virtual Real value(const Point & p);
  
private:
	Real _TC_IC;
};

#endif //COLUMNTEMPERATUREIC_H
