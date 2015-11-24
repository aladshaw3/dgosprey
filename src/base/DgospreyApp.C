#include "DgospreyApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"

#include "LinearDrivingForce.h"
#include "BedProperties.h"
#include "AdsorbentProperties.h"
#include "FlowProperties.h"
#include "RetardedTimeDerivative.h"
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
#include "MAGPIE_HeatAccumulation.h"


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
	registerKernel(LinearDrivingForce);
	registerMaterial(BedProperties);
	registerMaterial(AdsorbentProperties);
	registerMaterial(FlowProperties);
	registerMaterial(MagpieAdsorbateProperties);
	registerKernel(RetardedTimeDerivative);
	registerKernel(BedWallHeatTransfer);
	registerKernel(WallAmbientHeatTransfer);
	registerKernel(WallHeatAccumulation);
	registerKernel(BedHeatAccumulation);

	registerAux(TotalColumnPressure);
	registerAux(MAGPIE_Adsorption);
	registerAux(MAGPIE_Perturbation);
	registerAux(MAGPIE_AdsorptionHeat);

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
	registerKernel(MAGPIE_HeatAccumulation);
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
