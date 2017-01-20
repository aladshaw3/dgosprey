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

/// Execute the simple user interface to read yaml files and create DGOSPREY input files
int exec_SimpleUI(const char *file)
{
	//Declarations
	int success = 0;
	SimpleUI sui;
	
	// Create a blank digital MOOSE input file with default and required information
	sui.createMooseBlank();
	
	//Read the input file
	success = sui.readInputFile(file);
	if (success != 0)
	{
		mError(read_error);
		return -1;
	}
	
	//create example file
	sui.createExample();
	
	//sui.DisplayInput();
	sui.DisplayOutput();
	
	return -1;
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

//Function to create a blank MOOSE input file
void SimpleUI::createMooseBlank()
{
	this->moose_input.addDocKey("GlobalParams");
	
	this->moose_input.addDocKey("Problem");
	this->moose_input.getDocument("Problem").addPair("coord_type","RZ");
	
	this->moose_input.addDocKey("Mesh");
	this->moose_input.getDocument("Mesh").addPair("type","GeneratedMesh");
	this->moose_input.getDocument("Mesh").addPair("dim","2");
	this->moose_input.getDocument("Mesh").addPair("nx","10");
	this->moose_input.getDocument("Mesh").addPair("ny","40");
	this->moose_input.getDocument("Mesh").addPair("xmin","0.0");
	this->moose_input.getDocument("Mesh").addPair("ymin","0.0");
	
	this->moose_input.addDocKey("Variables");
	this->moose_input.getDocument("Variables").addHeadKey("wall_temp");
	this->moose_input.getDocument("Variables").getHeader("wall_temp").addPair("order","CONSTANT");
	this->moose_input.getDocument("Variables").getHeader("wall_temp").addPair("family","MONOMIAL");
	this->moose_input.getDocument("Variables").addHeadKey("column_temp");
	this->moose_input.getDocument("Variables").getHeader("column_temp").addPair("order","CONSTANT");
	this->moose_input.getDocument("Variables").getHeader("column_temp").addPair("family","MONOMIAL");
	
	this->moose_input.addDocKey("AuxVariables");
	this->moose_input.getDocument("AuxVariables").addHeadKey("total_pressure");
	this->moose_input.getDocument("AuxVariables").getHeader("total_pressure").addPair("order","CONSTANT");
	this->moose_input.getDocument("AuxVariables").getHeader("total_pressure").addPair("family","MONOMIAL");
	this->moose_input.getDocument("AuxVariables").getHeader("total_pressure").addPair("initial_condition","101.35");
	this->moose_input.getDocument("AuxVariables").addHeadKey("ambient_temp");
	this->moose_input.getDocument("AuxVariables").getHeader("ambient_temp").addPair("order","CONSTANT");
	this->moose_input.getDocument("AuxVariables").getHeader("ambient_temp").addPair("family","MONOMIAL");
	
	this->moose_input.addDocKey("ICs");
	
	this->moose_input.addDocKey("Kernels");
	this->moose_input.getDocument("Kernels").addHeadKey("wallAccum");
	this->moose_input.getDocument("Kernels").getHeader("wallAccum").addPair("type","WallHeatAccumulation");
	this->moose_input.getDocument("Kernels").getHeader("wallAccum").addPair("variable","wall_temp");
	this->moose_input.getDocument("Kernels").addHeadKey("wall_bed_trans");
	this->moose_input.getDocument("Kernels").getHeader("wall_bed_trans").addPair("type","BedWallHeatTransfer");
	this->moose_input.getDocument("Kernels").getHeader("wall_bed_trans").addPair("variable","wall_temp");
	this->moose_input.getDocument("Kernels").getHeader("wall_bed_trans").addPair("coupled","column_temp");
	this->moose_input.getDocument("Kernels").addHeadKey("wall_amb_trans");
	this->moose_input.getDocument("Kernels").getHeader("wall_amb_trans").addPair("type","WallAmbientHeatTransfer");
	this->moose_input.getDocument("Kernels").getHeader("wall_amb_trans").addPair("variable","wall_temp");
	this->moose_input.getDocument("Kernels").getHeader("wall_amb_trans").addPair("coupled","ambient_temp");
	this->moose_input.getDocument("Kernels").addHeadKey("columnAccum");
	this->moose_input.getDocument("Kernels").getHeader("columnAccum").addPair("type","BedHeatAccumulation");
	this->moose_input.getDocument("Kernels").getHeader("columnAccum").addPair("variable","column_temp");
	this->moose_input.getDocument("Kernels").addHeadKey("columnConduction");
	this->moose_input.getDocument("Kernels").getHeader("columnConduction").addPair("type","GColumnHeatDispersion");
	this->moose_input.getDocument("Kernels").getHeader("columnConduction").addPair("variable","column_temp");
	this->moose_input.getDocument("Kernels").addHeadKey("columnAdvection");
	this->moose_input.getDocument("Kernels").getHeader("columnAdvection").addPair("type","GColumnHeatAdvection");
	this->moose_input.getDocument("Kernels").getHeader("columnAdvection").addPair("variable","column_temp");
	
	this->moose_input.addDocKey("DGKernels");
	this->moose_input.getDocument("DGKernels").addHeadKey("DGcolumnConduction");
	this->moose_input.getDocument("DGKernels").getHeader("DGcolumnConduction").addPair("type","DGColumnHeatDispersion");
	this->moose_input.getDocument("DGKernels").getHeader("DGcolumnConduction").addPair("variable","column_temp");
	this->moose_input.getDocument("DGKernels").addHeadKey("DGcolumnAdvection");
	this->moose_input.getDocument("DGKernels").getHeader("DGcolumnAdvection").addPair("type","DGColumnHeatAdvection");
	this->moose_input.getDocument("DGKernels").getHeader("DGcolumnAdvection").addPair("variable","column_temp");
	
	this->moose_input.addDocKey("AuxKernels");
	this->moose_input.getDocument("AuxKernels").addHeadKey("column_pressure");
	this->moose_input.getDocument("AuxKernels").getHeader("column_pressure").addPair("type","TotalColumnPressure");
	this->moose_input.getDocument("AuxKernels").getHeader("column_pressure").addPair("variable","total_pressure");
	this->moose_input.getDocument("AuxKernels").getHeader("column_pressure").addPair("temperature","column_temp");
	
	this->moose_input.addDocKey("BCs");
	this->moose_input.getDocument("BCs").addHeadKey("Heat_Gas_Flux");
	this->moose_input.getDocument("BCs").getHeader("Heat_Gas_Flux").addPair("variable","column_temp");
	this->moose_input.getDocument("BCs").getHeader("Heat_Gas_Flux").addPair("boundary","'top bottom'");
	this->moose_input.getDocument("BCs").addHeadKey("Heat_Wall_Flux");
	this->moose_input.getDocument("BCs").getHeader("Heat_Wall_Flux").addPair("variable","column_temp");
	this->moose_input.getDocument("BCs").getHeader("Heat_Wall_Flux").addPair("boundary","'right left'");
	this->moose_input.getDocument("BCs").getHeader("Heat_Wall_Flux").addPair("wall_temp","wall_temp");
	
	
	
	this->moose_input.addDocKey("Materials");
	this->moose_input.getDocument("Materials").addHeadKey("BedMaterials");
	this->moose_input.getDocument("Materials").getHeader("BedMaterials").addPair("type","BedProperties");
	this->moose_input.getDocument("Materials").getHeader("BedMaterials").addPair("block","0");
	this->moose_input.getDocument("Materials").getHeader("BedMaterials").addPair("temperature","column_temp");
	
	this->moose_input.getDocument("Materials").addHeadKey("FlowMaterials");
	this->moose_input.getDocument("Materials").getHeader("FlowMaterials").addPair("type","FlowProperties");
	this->moose_input.getDocument("Materials").getHeader("FlowMaterials").addPair("block","0");
	this->moose_input.getDocument("Materials").getHeader("FlowMaterials").addPair("temperature","column_temp");
	this->moose_input.getDocument("Materials").getHeader("FlowMaterials").addPair("total_pressure","total_pressure");
	
	this->moose_input.getDocument("Materials").addHeadKey("AdsorbentMaterials");
	this->moose_input.getDocument("Materials").getHeader("AdsorbentMaterials").addPair("type","AdsorbentProperties");
	this->moose_input.getDocument("Materials").getHeader("AdsorbentMaterials").addPair("block","0");
	this->moose_input.getDocument("Materials").getHeader("AdsorbentMaterials").addPair("temperature","column_temp");
	
	this->moose_input.getDocument("Materials").addHeadKey("AdsorbateMaterials");
	this->moose_input.getDocument("Materials").getHeader("AdsorbateMaterials").addPair("type","MagpieAdsorbateProperties");
	this->moose_input.getDocument("Materials").getHeader("AdsorbateMaterials").addPair("block","0");
	this->moose_input.getDocument("Materials").getHeader("AdsorbateMaterials").addPair("temperature","column_temp");
	this->moose_input.getDocument("Materials").getHeader("AdsorbateMaterials").addPair("total_pressure","total_pressure");
	
	this->moose_input.addDocKey("Postprocessors");
	
	this->moose_input.addDocKey("Executioner");
	this->moose_input.getDocument("Executioner").addPair("type","Transient");
	this->moose_input.getDocument("Executioner").addPair("scheme","implicit-euler");
	this->moose_input.getDocument("Executioner").addPair("nl_rel_tol","1e-6");
	this->moose_input.getDocument("Executioner").addPair("nl_abs_tol","1e-6");
	this->moose_input.getDocument("Executioner").addPair("nl_rel_step_tol","1e-10");
	this->moose_input.getDocument("Executioner").addPair("nl_abs_step_tol","1e-10");
	this->moose_input.getDocument("Executioner").addPair("l_tol","1e-6");
	this->moose_input.getDocument("Executioner").addPair("l_max_its","100");
	this->moose_input.getDocument("Executioner").addPair("nl_max_its","10");
	this->moose_input.getDocument("Executioner").addPair("solve_type","newton");
	this->moose_input.getDocument("Executioner").addPair("line_search","none");
	this->moose_input.getDocument("Executioner").addPair("start_time","0.0");
	this->moose_input.getDocument("Executioner").addPair("petsc_options_iname","'-pc_type -pc_hypre_type -ksp_gmres_restart'");
	this->moose_input.getDocument("Executioner").addPair("petsc_options_value","'hypre boomeramg 100'");
	
	this->moose_input.getDocument("Executioner").addHeadKey("TimeStepper");
	
	this->moose_input.addDocKey("Preconditioning");
	
	this->moose_input.addDocKey("Outputs");
	this->moose_input.getDocument("Outputs").addPair("exodus","true");
	this->moose_input.getDocument("Outputs").addPair("csv","true");
	this->moose_input.getDocument("Outputs").addPair("print_linear_residuals","true");
}

//Function to create a MOOSE example input file
void SimpleUI::createExample()
{
	this->moose_input.getDocument("GlobalParams").addPair("initial_dt","0.01");
	this->moose_input.getDocument("GlobalParams").addPair("length","22.86");
	
	this->moose_input.getDocument("Mesh").addPair("xmax","0.8636");
	this->moose_input.getDocument("Mesh").addPair("ymax","22.86");
	
	this->moose_input.getDocument("Variables").addHeadKey("Kr");
	this->moose_input.getDocument("Variables").getHeader("Kr").addPair("order","CONSTANT");
	this->moose_input.getDocument("Variables").getHeader("Kr").addPair("family","MONOMIAL");
	
	this->moose_input.getDocument("Variables").addHeadKey("Xe");
	this->moose_input.getDocument("Variables").getHeader("Xe").addPair("order","CONSTANT");
	this->moose_input.getDocument("Variables").getHeader("Xe").addPair("family","MONOMIAL");
	
	this->moose_input.getDocument("Variables").addHeadKey("He");
	this->moose_input.getDocument("Variables").getHeader("He").addPair("order","CONSTANT");
	this->moose_input.getDocument("Variables").getHeader("He").addPair("family","MONOMIAL");
	
	this->moose_input.getDocument("Variables").getHeader("wall_temp").addPair("initial_condition","253.15");
	this->moose_input.getDocument("Variables").getHeader("column_temp").addPair("initial_condition","253.15");
	
	this->moose_input.getDocument("AuxVariables").getHeader("ambient_temp").addPair("initial_condition","253.15");
	
	this->moose_input.getDocument("AuxVariables").addHeadKey("He_Adsorbed");
	this->moose_input.getDocument("AuxVariables").getHeader("He_Adsorbed").addPair("order","CONSTANT");
	this->moose_input.getDocument("AuxVariables").getHeader("He_Adsorbed").addPair("family","MONOMIAL");
	this->moose_input.getDocument("AuxVariables").getHeader("He_Adsorbed").addPair("initial_condition","0.0");
	
	this->moose_input.getDocument("AuxVariables").addHeadKey("Kr_Adsorbed");
	this->moose_input.getDocument("AuxVariables").getHeader("Kr_Adsorbed").addPair("order","CONSTANT");
	this->moose_input.getDocument("AuxVariables").getHeader("Kr_Adsorbed").addPair("family","MONOMIAL");
	this->moose_input.getDocument("AuxVariables").getHeader("Kr_Adsorbed").addPair("initial_condition","0.0");
	
	this->moose_input.getDocument("AuxVariables").addHeadKey("Xe_Adsorbed");
	this->moose_input.getDocument("AuxVariables").getHeader("Xe_Adsorbed").addPair("order","CONSTANT");
	this->moose_input.getDocument("AuxVariables").getHeader("Xe_Adsorbed").addPair("family","MONOMIAL");
	this->moose_input.getDocument("AuxVariables").getHeader("Xe_Adsorbed").addPair("initial_condition","0.0");
	
	this->moose_input.getDocument("AuxVariables").addHeadKey("He_Perturb");
	this->moose_input.getDocument("AuxVariables").getHeader("He_Perturb").addPair("order","CONSTANT");
	this->moose_input.getDocument("AuxVariables").getHeader("He_Perturb").addPair("family","MONOMIAL");
	this->moose_input.getDocument("AuxVariables").getHeader("He_Perturb").addPair("initial_condition","0.0");
	
	this->moose_input.getDocument("AuxVariables").addHeadKey("Kr_Perturb");
	this->moose_input.getDocument("AuxVariables").getHeader("Kr_Perturb").addPair("order","CONSTANT");
	this->moose_input.getDocument("AuxVariables").getHeader("Kr_Perturb").addPair("family","MONOMIAL");
	this->moose_input.getDocument("AuxVariables").getHeader("Kr_Perturb").addPair("initial_condition","0.0");
	
	this->moose_input.getDocument("AuxVariables").addHeadKey("Xe_Perturb");
	this->moose_input.getDocument("AuxVariables").getHeader("Xe_Perturb").addPair("order","CONSTANT");
	this->moose_input.getDocument("AuxVariables").getHeader("Xe_Perturb").addPair("family","MONOMIAL");
	this->moose_input.getDocument("AuxVariables").getHeader("Xe_Perturb").addPair("initial_condition","0.0");
	
	this->moose_input.getDocument("AuxVariables").addHeadKey("He_AdsorbedHeat");
	this->moose_input.getDocument("AuxVariables").getHeader("He_AdsorbedHeat").addPair("order","CONSTANT");
	this->moose_input.getDocument("AuxVariables").getHeader("He_AdsorbedHeat").addPair("family","MONOMIAL");
	this->moose_input.getDocument("AuxVariables").getHeader("He_AdsorbedHeat").addPair("initial_condition","0.0");
	
	this->moose_input.getDocument("AuxVariables").addHeadKey("Kr_AdsorbedHeat");
	this->moose_input.getDocument("AuxVariables").getHeader("Kr_AdsorbedHeat").addPair("order","CONSTANT");
	this->moose_input.getDocument("AuxVariables").getHeader("Kr_AdsorbedHeat").addPair("family","MONOMIAL");
	this->moose_input.getDocument("AuxVariables").getHeader("Kr_AdsorbedHeat").addPair("initial_condition","0.0");
	
	this->moose_input.getDocument("AuxVariables").addHeadKey("Xe_AdsorbedHeat");
	this->moose_input.getDocument("AuxVariables").getHeader("Xe_AdsorbedHeat").addPair("order","CONSTANT");
	this->moose_input.getDocument("AuxVariables").getHeader("Xe_AdsorbedHeat").addPair("family","MONOMIAL");
	this->moose_input.getDocument("AuxVariables").getHeader("Xe_AdsorbedHeat").addPair("initial_condition","0.0");
	
	this->moose_input.getDocument("ICs").addHeadKey("Kr_IC");
	this->moose_input.getDocument("ICs").getHeader("Kr_IC").addPair("type","ConcentrationIC");
	this->moose_input.getDocument("ICs").getHeader("Kr_IC").addPair("variable","Kr");
	this->moose_input.getDocument("ICs").getHeader("Kr_IC").addPair("initial_mole_frac","0.0");
	this->moose_input.getDocument("ICs").getHeader("Kr_IC").addPair("initial_press","101.35");
	this->moose_input.getDocument("ICs").getHeader("Kr_IC").addPair("initial_temp","253.15");
	
	this->moose_input.getDocument("ICs").addHeadKey("Xe_IC");
	this->moose_input.getDocument("ICs").getHeader("Xe_IC").addPair("type","ConcentrationIC");
	this->moose_input.getDocument("ICs").getHeader("Xe_IC").addPair("variable","Xe");
	this->moose_input.getDocument("ICs").getHeader("Xe_IC").addPair("initial_mole_frac","0.0");
	this->moose_input.getDocument("ICs").getHeader("Xe_IC").addPair("initial_press","101.35");
	this->moose_input.getDocument("ICs").getHeader("Xe_IC").addPair("initial_temp","253.15");
	
	this->moose_input.getDocument("ICs").addHeadKey("He_IC");
	this->moose_input.getDocument("ICs").getHeader("He_IC").addPair("type","ConcentrationIC");
	this->moose_input.getDocument("ICs").getHeader("He_IC").addPair("variable","He");
	this->moose_input.getDocument("ICs").getHeader("He_IC").addPair("initial_mole_frac","1.0");
	this->moose_input.getDocument("ICs").getHeader("He_IC").addPair("initial_press","101.35");
	this->moose_input.getDocument("ICs").getHeader("He_IC").addPair("initial_temp","253.15");
	
	this->moose_input.getDocument("Kernels").addHeadKey("Kr_accum");
	this->moose_input.getDocument("Kernels").getHeader("Kr_accum").addPair("type","BedMassAccumulation");
	this->moose_input.getDocument("Kernels").getHeader("Kr_accum").addPair("variable","Kr");
	this->moose_input.getDocument("Kernels").getHeader("Kr_accum").addPair("index","0");
	this->moose_input.getDocument("Kernels").addHeadKey("Kr_masstrans");
	this->moose_input.getDocument("Kernels").getHeader("Kr_masstrans").addPair("type","AdsorptionMassTransfer");
	this->moose_input.getDocument("Kernels").getHeader("Kr_masstrans").addPair("variable","Kr");
	this->moose_input.getDocument("Kernels").getHeader("Kr_masstrans").addPair("solid_conc","Kr_Adsorbed");
	this->moose_input.getDocument("Kernels").addHeadKey("Kr_diff");
	this->moose_input.getDocument("Kernels").getHeader("Kr_diff").addPair("type","GColumnMassDispersion");
	this->moose_input.getDocument("Kernels").getHeader("Kr_diff").addPair("variable","Kr");
	this->moose_input.getDocument("Kernels").getHeader("Kr_diff").addPair("index","0");
	this->moose_input.getDocument("Kernels").addHeadKey("Kr_adv");
	this->moose_input.getDocument("Kernels").getHeader("Kr_adv").addPair("type","GColumnMassAdvection");
	this->moose_input.getDocument("Kernels").getHeader("Kr_adv").addPair("variable","Kr");
	
	this->moose_input.getDocument("Kernels").addHeadKey("Xe_accum");
	this->moose_input.getDocument("Kernels").getHeader("Xe_accum").addPair("type","BedMassAccumulation");
	this->moose_input.getDocument("Kernels").getHeader("Xe_accum").addPair("variable","Xe");
	this->moose_input.getDocument("Kernels").getHeader("Xe_accum").addPair("index","1");
	this->moose_input.getDocument("Kernels").addHeadKey("Xe_masstrans");
	this->moose_input.getDocument("Kernels").getHeader("Xe_masstrans").addPair("type","AdsorptionMassTransfer");
	this->moose_input.getDocument("Kernels").getHeader("Xe_masstrans").addPair("variable","Xe");
	this->moose_input.getDocument("Kernels").getHeader("Xe_masstrans").addPair("solid_conc","Xe_Adsorbed");
	this->moose_input.getDocument("Kernels").addHeadKey("Xe_diff");
	this->moose_input.getDocument("Kernels").getHeader("Xe_diff").addPair("type","GColumnMassDispersion");
	this->moose_input.getDocument("Kernels").getHeader("Xe_diff").addPair("variable","Xe");
	this->moose_input.getDocument("Kernels").getHeader("Xe_diff").addPair("index","1");
	this->moose_input.getDocument("Kernels").addHeadKey("Xe_adv");
	this->moose_input.getDocument("Kernels").getHeader("Xe_adv").addPair("type","GColumnMassAdvection");
	this->moose_input.getDocument("Kernels").getHeader("Xe_adv").addPair("variable","Xe");
	
	this->moose_input.getDocument("Kernels").addHeadKey("He_accum");
	this->moose_input.getDocument("Kernels").getHeader("He_accum").addPair("type","BedMassAccumulation");
	this->moose_input.getDocument("Kernels").getHeader("He_accum").addPair("variable","He");
	this->moose_input.getDocument("Kernels").getHeader("He_accum").addPair("index","2");
	this->moose_input.getDocument("Kernels").addHeadKey("He_masstrans");
	this->moose_input.getDocument("Kernels").getHeader("He_masstrans").addPair("type","AdsorptionMassTransfer");
	this->moose_input.getDocument("Kernels").getHeader("He_masstrans").addPair("variable","He");
	this->moose_input.getDocument("Kernels").getHeader("He_masstrans").addPair("solid_conc","He_Adsorbed");
	this->moose_input.getDocument("Kernels").addHeadKey("He_diff");
	this->moose_input.getDocument("Kernels").getHeader("He_diff").addPair("type","GColumnMassDispersion");
	this->moose_input.getDocument("Kernels").getHeader("He_diff").addPair("variable","He");
	this->moose_input.getDocument("Kernels").getHeader("He_diff").addPair("index","2");
	this->moose_input.getDocument("Kernels").addHeadKey("He_adv");
	this->moose_input.getDocument("Kernels").getHeader("He_adv").addPair("type","GColumnMassAdvection");
	this->moose_input.getDocument("Kernels").getHeader("He_adv").addPair("variable","He");
	
	this->moose_input.getDocument("Kernels").addHeadKey("column_AdsHeat");
	this->moose_input.getDocument("Kernels").getHeader("column_AdsHeat").addPair("type","AdsorptionHeatAccumulation");
	this->moose_input.getDocument("Kernels").getHeader("column_AdsHeat").addPair("variable","column_temp");
	this->moose_input.getDocument("Kernels").getHeader("column_AdsHeat").addPair("solid_heats","'Kr_AdsorbedHeat Xe_AdsorbedHeat He_AdsorbedHeat'");
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
