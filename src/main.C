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
#include "MooseInit.h"
#include "Moose.h"
#include "MooseApp.h"
#include "AppFactory.h"
#include "SimpleUI.h"

// Create a performance log
PerfLog Moose::perf_log("Dgosprey");

// Begin the main program.
int main(int argc, char *argv[])
{
	/**	Simplified command line interface for higher level control */
	
	//Initialize MPI environment
	int pid;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	
	//Only valid for three input arguments (simple controls)
	if (argc == 3)
	{
		//Parse 3rd argument to determine file extension
		if (pid == 0)
		{
			std::cout << "\nDo stuff before MOOSE...\n";
			std::cout << argc << std::endl;
			std::cout << argv[0] << std::endl;
			std::cout << argv[1] << std::endl;
			std::cout << argv[2] << std::endl;
		}
	}
	
	// Initialize MPI, solvers and MOOSE
	MooseInit init(argc, argv);

	// Register this application's MooseApp and any it depends on
	DgospreyApp::registerApps();

	// This creates dynamic memory that we're responsible for deleting
	MooseApp * app = AppFactory::createApp("DgospreyApp", argc, argv);

	app->legacyUoInitializationDefault() = false;
	app->legacyUoAuxComputationDefault() = false;

	// Execute the application
	app->run();

	// Free up the memory we created earlier
	delete app;

	return 0;
}
