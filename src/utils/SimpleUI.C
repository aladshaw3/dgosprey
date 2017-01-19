/*!
 *  \file SimpleUI.h
 *	\brief Utilities kernel to provide a simpler user interface through yaml input files
 *  \author Austin Ladshaw
 *	\date 01/19/2017
 *	\copyright This kernel was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science and was developed for use
 *			   by Idaho National Laboratory and Oak Ridge National Laboratory
 *			   engineers and scientists. Portions Copyright (c) 2017, all
 *             rights reserved.
 *
 *			   Austin Ladshaw does not claim any ownership or copyright to the
 *			   MOOSE framework in which these kernels are constructed, only
 *			   the kernels themselves. The MOOSE framework copyright is held
 *			   by the Battelle Energy Alliance, LLC (c) 2010, all rights reserved.
 */

#include "SimpleUI.h"

/// Function to return true if the file extension of the given data is .yml
bool isYamlFile(char argv[])
{
	std::string arg = argv;
	if((arg.substr(arg.find_last_of(".") + 1) == "yml") || (arg.substr(arg.find_last_of(".") + 1) == "yaml"))
		return true;
	else
		return false;
}