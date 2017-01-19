/*!
 *  \file SimpleUI.h
 *	\brief Utilities kernel to provide a simpler user interface through yaml input files
 *	\details This file is responsible for providing users with a very simple, optional
 *				user interface that utilizes the yaml file structure. It is coupled with
 *				other data base files, stored in the project directory, that holds parameter
 *				information for common systems that the user may wish to model. The intent
 *				is to give an easier way for the uneducated user to be capable of utilizing
 *				this software without any advanced understanding or manual. 
 *
 *				The simple user interface will read input from a yaml (.yml) file and then
 *				read in the necessary data base information to create a standard MOOSE
 *				input file (.i) that the DGOSPREY program will run. As such, this interface
 *				will have significantly less control over the options that the lower level
 *				interface will have.
 *
 *
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

#include "yaml_wrapper.h"

#ifndef SimpleUI_h
#define SimpleUI_h

/// Function to return true if the file extension of the given data is .yml
bool isYamlFile(char argv[]);

#endif /* SimpleUI_h */