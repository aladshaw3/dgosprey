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

/// Create a blank digital MOOSE input file for DGOSPREY
void createMooseBlank(YamlWrapper *moose_input)
{
	moose_input->addDocKey("GlobalParams");
	
	moose_input->addDocKey("Problem");
	moose_input->getDocument("Problem").addPair("coord_type","RZ");
	
	moose_input->addDocKey("Mesh");
	moose_input->getDocument("Mesh").addPair("type","GeneratedMesh");
	moose_input->getDocument("Mesh").addPair("dim","2");
	moose_input->getDocument("Mesh").addPair("xmin","0.0");
	moose_input->getDocument("Mesh").addPair("ymin","0.0");
	
	moose_input->addDocKey("Variables");
	
	moose_input->addDocKey("AuxVariables");
	
	moose_input->addDocKey("ICs");
	
	moose_input->addDocKey("Kernels");
	
	moose_input->addDocKey("DGKernels");
	
	moose_input->addDocKey("AuxKernels");
	
	moose_input->addDocKey("BCs");
	
	moose_input->addDocKey("Materials");
	moose_input->getDocument("Materials").addHeadKey("BedMaterials");
	moose_input->getDocument("Materials").addHeadKey("FlowMaterials");
	moose_input->getDocument("Materials").addHeadKey("AdsorbentMaterials");
	moose_input->getDocument("Materials").addHeadKey("AdsorbateMaterials");
	
	moose_input->addDocKey("Postprocessors");
	
	moose_input->addDocKey("Executioner");
	moose_input->getDocument("Executioner").addHeadKey("TimeStepper");
	
	moose_input->addDocKey("Preconditioning");
	
	moose_input->addDocKey("Outputs");
}

/// Execute the simple user interface to read yaml files and create DGOSPREY input files
void exec_SimpleUI(const char *file)
{
	//Declarations
	int success = 0;
	SimpleUI sui;
	
	// Create a blank digital MOOSE input file
	createMooseBlank(&sui.getMooseInput());
	
	//Read the input file
	success = sui.readInputFile(file);
	if (success != 0)
	{
		mError(read_error);
		return; //May return blank moose doc
	}
	
	sui.DisplayInput();
	sui.DisplayOutput();
	
}

//Default constructor
SimpleUI::SimpleUI()
{
}

//Default destructor
SimpleUI::~SimpleUI()
{
}

//Function to read a given input file
int SimpleUI::readInputFile(const char *file)
{
	return this->yaml_input.executeYamlRead(file);
}

//Return reference to the YamlWrapper object for the input file
YamlWrapper& SimpleUI::getYamlInput()
{
	return this->yaml_input.getYamlWrapper();
}

//Return reference to the YamlWrapper object to create the MOOSE input file
YamlWrapper& SimpleUI::getMooseInput()
{
	return this->moose_input;
}

//Function to display input file contents
void SimpleUI::DisplayInput()
{
	this->yaml_input.DisplayContents();
}

//Function to display MOOSE input file contents that we are making
void SimpleUI::DisplayOutput()
{
	this->moose_input.DisplayContents();
}
