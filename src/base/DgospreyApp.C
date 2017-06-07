/*!
 *  \file DgospreyApp.h
 *	\brief Registration object for creating a registering DGOSPREY kernels
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

#include "DgospreyApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

#include "LinearDrivingForce.h"
#include "BedProperties.h"
#include "AdsorbentProperties.h"
#include "FlowProperties.h"
#include "BedMassAccumulation.h"
#include "BedWallHeatTransfer.h"
#include "WallAmbientHeatTransfer.h"
#include "WallHeatAccumulation.h"
#include "BedHeatAccumulation.h"
#include "TotalColumnPressure.h"
#include "TotalPressureIC.h"
#include "ColumnTemperatureIC.h"
#include "ConcentrationIC.h"
#include "DGAdvection.h"
#include "DGFluxBC.h"
#include "GAdvection.h"
#include "DGAnisotropicDiffusion.h"
#include "GAnisotropicDiffusion.h"
#include "DGFluxLimitedBC.h"

#include "DGColumnMassAdvection.h"
#include "DGColumnMassDispersion.h"
#include "DGMassFluxBC.h"
#include "DGMassFluxLimitedBC.h"
#include "GColumnMassAdvection.h"
#include "GColumnMassDispersion.h"

#include "DGColumnHeatAdvection.h"
#include "GColumnHeatAdvection.h"
#include "DGColumnHeatDispersion.h"
#include "GColumnHeatDispersion.h"
#include "DGHeatFluxBC.h"
#include "DGHeatFluxLimitedBC.h"
#include "DGColumnWallHeatFluxBC.h"
#include "DGColumnWallHeatFluxLimitedBC.h"

#include "MagpieAdsorbateProperties.h"
#include "MAGPIE_Adsorption.h"
#include "MAGPIE_Perturbation.h"
#include "MAGPIE_AdsorptionHeat.h"
#include "AdsorptionHeatAccumulation.h"
#include "AdsorptionMassTransfer.h"

#include "Aux_LDF.h"
#include "MAGPIE_ConstLDF_Adsorption.h"
#include "MAGPIE_ConstLDF_Perturbation.h"
#include "MAGPIE_MaterialLDF_Adsorption.h"
#include "MAGPIE_MaterialLDF_Perturbation.h"

#include "ScopsowlProperties.h"
#include "Scopsowl_Adsorption.h"
#include "Scopsowl_Perturbation.h"

#include "TotalMoles.h"
#include "CoupledCoeffTimeDerivative.h"
#include "CoupledLinearForcingFunction.h"
#include "CoupledLinearLDF.h"

#include "ThermodynamicProperties.h"
#include "GasFlowProperties.h"
#include "KineticProperties.h"
#include "HeatofAdsorption.h"
#include "SolidMassTransfer.h"
#include "SolidHeatTransfer.h"

#include "WallTemperature.h"
#include "CoupledLangmuirForcingFunction.h"

template<>
InputParameters validParams<DgospreyApp>()
{
  InputParameters params = validParams<MooseApp>();
  params.set<bool>("use_legacy_output_syntax") = false;
  return params;
}

DgospreyApp::DgospreyApp(InputParameters parameters) :
    MooseApp(parameters)
{

  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  DgospreyApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  DgospreyApp::associateSyntax(_syntax, _action_factory);
}

DgospreyApp::~DgospreyApp()
{
}

extern "C" void DgospreyApp__registerApps() { DgospreyApp::registerApps(); }
void
DgospreyApp::registerApps()
{
	registerApp(DgospreyApp);
}

void
DgospreyApp::registerObjects(Factory & factory)
{
	registerMaterial(BedProperties);
	registerMaterial(AdsorbentProperties);
	registerMaterial(FlowProperties);
	registerMaterial(MagpieAdsorbateProperties);
	registerMaterial(ScopsowlProperties);
	registerMaterial(ThermodynamicProperties);
	registerMaterial(GasFlowProperties);
	registerMaterial(KineticProperties);
	
	registerKernel(LinearDrivingForce);
	registerKernel(BedMassAccumulation);
	registerKernel(BedWallHeatTransfer);
	registerKernel(WallAmbientHeatTransfer);
	registerKernel(WallHeatAccumulation);
	registerKernel(BedHeatAccumulation);
	registerKernel(AdsorptionMassTransfer);
	registerKernel(CoupledCoeffTimeDerivative);
	registerKernel(CoupledLinearForcingFunction);
	registerKernel(CoupledLinearLDF);
	registerKernel(HeatofAdsorption);
	registerKernel(SolidMassTransfer);
	registerKernel(SolidHeatTransfer);
    registerKernel(CoupledLangmuirForcingFunction);
	
	registerAux(TotalColumnPressure);
	registerAux(MAGPIE_Adsorption);
	registerAux(MAGPIE_Perturbation);
	registerAux(MAGPIE_AdsorptionHeat);
	registerAux(Aux_LDF);
	registerAux(MAGPIE_ConstLDF_Adsorption);
	registerAux(MAGPIE_ConstLDF_Perturbation);
	registerAux(MAGPIE_MaterialLDF_Adsorption);
	registerAux(MAGPIE_MaterialLDF_Perturbation);
	registerAux(Scopsowl_Adsorption);
	registerAux(Scopsowl_Perturbation);
	registerAux(TotalMoles);
	registerAux(WallTemperature);

	registerInitialCondition(TotalPressureIC);
	registerInitialCondition(ColumnTemperatureIC);
	registerInitialCondition(ConcentrationIC);

	registerDGKernel(DGAdvection);
	registerBoundaryCondition(DGFluxBC);
	registerKernel(GAdvection);
	registerDGKernel(DGAnisotropicDiffusion);
	registerKernel(GAnisotropicDiffusion);
	registerBoundaryCondition(DGFluxLimitedBC);

	registerDGKernel(DGColumnMassAdvection);
	registerDGKernel(DGColumnMassDispersion);
	registerBoundaryCondition(DGMassFluxBC);
	registerBoundaryCondition(DGMassFluxLimitedBC);
	registerKernel(GColumnMassAdvection);
	registerKernel(GColumnMassDispersion);

	registerDGKernel(DGColumnHeatAdvection);
	registerKernel(GColumnHeatAdvection);
	registerDGKernel(DGColumnHeatDispersion);
	registerKernel(GColumnHeatDispersion);
	registerKernel(AdsorptionHeatAccumulation);
	registerBoundaryCondition(DGHeatFluxBC);
	registerBoundaryCondition(DGHeatFluxLimitedBC);
	registerBoundaryCondition(DGColumnWallHeatFluxBC);
	registerBoundaryCondition(DGColumnWallHeatFluxLimitedBC);
}

void
DgospreyApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
	//Register Actions Here
}
