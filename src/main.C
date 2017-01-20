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
	int success = 0;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	
	//Only valid for three input arguments (simple controls)
	if (argc == 3)
	{
		//Parse 3rd argument to determine file extension
		if (isYamlFile(argv[2]) == true)
		{
			//Use a single node to read the yaml input file
			if (pid == 0)
				success = exec_SimpleUI(argv[2]);
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(&success,1,MPI_INT,0,MPI_COMM_WORLD);
			
			if (success == 0)
			{
				//Run MOOSE immediately, using the input file created
				if (pid == 0)
				{
					std::cout << "\nCreated a MOOSE input file from given Yaml file!\n\n";
					std::cout << "Executing DGOSPREY application with the created input file...\n\n";
				}
			}
			else if (success == 1)
			{
				//Do not run MOOSE, send message of completed input file creation
				if (pid == 0)
				{
					std::cout << "\nCreated a MOOSE input file from given Yaml file!\n\n";
					std::cout << "You may open that file and edit it or run it with DGSOSPREY application...\n\n";
				}
				MPI_Finalize();
				return 0;
			}
			else
			{
				//Some error has occured
				if (pid == 0)
					std::cout << "\nERROR!!! Failed to create MOOSE input file from given Yaml file!\n\n";
				MPI_Finalize();
				return 0;
			}
		}//If not a yaml file, then continue with low level MOOSE interface
	}
	//Synchronization step to ensure that other processors don't start before pid 0  is finished
	MPI_Barrier(MPI_COMM_WORLD);
	
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
