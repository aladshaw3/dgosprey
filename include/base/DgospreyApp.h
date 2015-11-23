/*!
 *  \file DgospreyApp.h
 *	\brief Registration object for creating a registering DGOSPREY kernels
 *	\details This file is responsible for registering all DGOSPREY kernels in the MOOSE
 *			framework. Any additional kernel developed under DGOSPREY must be included
 *			and registered in this structure. This structure is required by MOOSE in
 *			order to have the DGOSPREY objects interact with the underlying MOOSE solvers
 *			and finite element framework.
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

#ifndef DGOSPREYAPP_H
#define DGOSPREYAPP_H

#include "MooseApp.h"

/// DgospreyApp class object forward declaration
class DgospreyApp;

template<>
InputParameters validParams<DgospreyApp>();

/// DgospreyApp inherits from the MooseApp object
/** This object defines the required constructors, destructors, and functions that must
	be a part of every MooseApp based object. All MooseApp objects must be created in this
	way and override these functions. */
class DgospreyApp : public MooseApp
{
public:
	/// DgospreyApp constructor (required)
	DgospreyApp(InputParameters parameters);
	/// DgospreyApp destructor (requried)
	virtual ~DgospreyApp();

	/// Function to the DgospreyApp into MOOSE (required)
	static void registerApps();
	/// Function to register kernels/objects created in DGOSPREY into the application (required)
	/** This is the function where the user must register all the kernels and other modules that are to
		be used in DGOSPREY. Each time a new kernel or other object is created in DGOSPREY, it must be
		registered here prior to building and running the application. Otherwise, the new functionallity
		added will not show up or be utilized. */
	static void registerObjects(Factory & factory);
	/// Function to associate syntax with the DgospreyApp (required?)
	/** I don't know what this is or does or what is actually being registered. */
	static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* DGOSPREYAPP_H */
