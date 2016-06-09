/*!
 *  \file scopsowl.h scopsowl.cpp
 *	\brief Simultaneously Coupled Objects for Pore and Surface diffusion Operations With Linear systems
 *	\details This file contains structures and functions associated with modeling
 *			adsorption in commercial, bi-porous adsorbents such as zeolites and
 *			mordenites. The pore diffusion and mass transfer equations are coupled
 *			with adsorption and surface diffusion through smaller crystals embedded
 *			in a binder matrix. However, you can also direct this simulation to treat
 *			the adsorbent as homogeneous (instead of heterogeneous) in order to model
 *			an even greater variety of gaseous adsorption kinetic problems. This object
 *			is coupled with either MAGPIE, SKUA, or BOTH depending on the type of
 *			simulation requested.
 *
 *  \author Austin Ladshaw
 *	\date 01/29/2015
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science. Copyright (c) 2015, all
 *             rights reserved.
 */

#include "egret.h"							//EGRET handles the estimation of gas-phase properties
#include "skua.h"							//SKUA couples MAGPIE and Surface kinetics together

#ifndef SCOPSOWL_HPP_
#define SCOPSOWL_HPP_

#ifndef Dp
#define Dp(Dm,ep) (ep*ep*Dm)					///< Estimate of Pore Diffusivity (cm^2/s)
#endif

#ifndef Dk
#define Dk(rp,T,MW) (9700.0*rp*pow((T/MW),0.5))	///< Estimate of Knudsen Diffusivity (cm^2/s)
#endif

#ifndef avgDp
#define avgDp(Dp,Dk) (pow(((1/Dp)+(1/Dk)),-1.0))	///< Estimate of Average Pore Diffusion (cm^2/s)
#endif

/// Data structure for the species' parameters in SCOPSOWL
/** C-style object that holds information on all species for a particular SCOPSOWL
	simulation. Initial conditions, kinetic parameters, and interim matrix objects
	are stored here for use in various SCOSPSOWL functions. */
typedef struct
{
	Matrix<double> qAvg;					///< Average adsorbed amount for a species at each node (mol/kg)
	Matrix<double> qAvg_old;				///< Old Average adsorbed amount for a species at each node (mol/kg)
	
	Matrix<double> Qst;						///< Heat of adsorption for all nodes (J/mol)
	Matrix<double> Qst_old;					///< Old Heat of adsorption for all nodes (J/mol)
	
	Matrix<double> dq_dc;					///< Storage vector for current adsorption slope/strength (dq/dc) (L/kg)
	
	double xIC;								///< Initial conditions for adsorbed molefractions
	
	double qIntegralAvg;					///< Integral average of adsorption over the entire pellet (mol/kg)
	double qIntegralAvg_old;				///< Old Integral average of adsorption over the entire pellet (mol/kg)
	
	double QstAvg;							///< Integral average heat of adsorption (J/mol)
	double QstAvg_old;						///< Old integral average heat of adsorption (J/mol)
	
	double qo;								///< Boundary value of adsorption if using Dirichlet BCs (mol/kg)
	double Qsto;							///< Boundary value of adsorption heat if using Dirichlet BCs (J/mol)
	double dq_dco;							///< Boundary value of adsorption slope for Dirichelt BCs (L/kg)
	
	double pore_diffusion;					///< Value for constant pore diffusion (cm^2/hr)
	double film_transfer;					///< Value for constant film mass transfer (cm/hr)
	
	double activation_energy;				///< Activation energy for surface diffusion (J/mol)
	double ref_diffusion;					///< Reference state surface diffusivity (um^2/hr)
	double ref_temperature;					///< Reference temperature for empirical adjustments (K)
	double affinity;						///< Affinity parameter used in empirical adjustments (-)
	double ref_pressure;					///< Reference pressure used in empirical adjustments (kPa)
	
	bool Adsorbable;						///< True = species can adsorb; False = species cannot adsorb
	
	std::string speciesName;				///< String to hold the name of each species
	
}SCOPSOWL_PARAM_DATA;

/// Primary data structure for SCOPSOWL simulations
/** C-style object holding necessary information to run a SCOPSOWL simulation.
	SCOPSOWL is a multi-scale problem involving PDE solution for the macro-scale
	adsorbent pellet and the micro-scale adsorbent crystals. As such, each SCOPSOWL
	simulation involves multiple SKUA simulations at the nodes in the macro-scale
	domain. Alternatively, if the user wishes to specify that the adsorbent is
	homogeneous, then you can run SCOPSOWL as a single-scale problem. Additionally,
	you can simplfy the model by assuming that the micro-scale diffusion is very
	fast, and therefore replace each SKUA simulation with a simpler MAGPIE evaluation.
	Details on running SCOPSOWL with the various options will be discussed in the
	SCOPSOWL_SCENARIOS function. */
typedef struct
{
	unsigned long int total_steps;			///< Running total of all calculation steps
	int coord_macro;						///< Coordinate system for large pellet
	int coord_micro;						///< Coordinate system for small crystal (if any)
	int level = 2;							///< Level of coupling between the different scales (default = 2)
	double sim_time;						///< Stopping time for the simulation (hrs)
	double t_old;							///< Old time of the simulations (hrs)
	double t;								///< Current time of the simulations (hrs)
	double t_counter = 0.0;					///< Counter for the time output
	double t_print;							///< Print output at every t_print time (hrs)
	
	bool Print2File = true;					///< True = results to .txt; False = no printing
	bool Print2Console = true;				///< True = results to console; False = no printing
	bool SurfDiff = true;					///< True = includes SKUA simulation if Heterogeneous; False = only uses MAGPIE
	bool Heterogeneous = true;				///< True = pellet is made of binder and crystals, False = all one phase
	
	double gas_velocity;					///< Superficial Gas Velocity arount pellet (cm/s)
	double total_pressure;					///< Gas phase total pressure (kPa)
	double gas_temperature;					///< Gas phase temperature (K)
	double pellet_radius;					///< Nominal radius of the pellet - macroscale domain (cm)
	double crystal_radius;					///< Nominal radius of the crystal - microscale domain (um)
	double char_macro;						///< Characteristic size for macro scale (cm or cm^2) - only if pellet is not spherical
	double char_micro;						///< Characteristic size for micro scale (um or um^2) - only if crystal is not spherical
	double binder_fraction;					///< Volume of binder per total volume of pellet (-)
	double binder_porosity;					///< Volume of pores per volume of binder (-)
	double binder_poresize;					///< Nominal radius of the binder pores (cm)
	double pellet_density;					///< Mass of the pellet per volume of pellet (kg/L)
	
	bool DirichletBC = false;				///< True = Dirichlet BC; False = Neumann BC
	bool NonLinear = true;					///< True = Non-linear solver; False = Linear solver
	std::vector<double> y;					///< Outside mole fractions of each component (-)
	std::vector<double> tempy;				///< Temporary place holder for gas mole fractions in other locations (-)
	
	FILE *OutputFile;											///< Output file pointer to the output file for postprocesses
	double (*eval_ads) (int i, int l, const void *user_data);	///< Function pointer for evaluating adsorption (mol/kg)
	double (*eval_retard) (int i, int l, const void *user_data);///< Function pointer for evaluating retardation (-)
	double (*eval_diff) (int i, int l, const void *user_data);	///< Function pointer for evaluating pore diffusion (cm^2/hr)
	double (*eval_surfDiff) (int i, int l, const void *user_data);///< Function pointer for evaluating surface diffusion (um^2/hr)
	double (*eval_kf) (int i, const void *user_data);			///< Function pointer for evaluating film mass transfer (cm/hr)
	
	const void *user_data;						///< Data structure for users info to calculate parameters
	MIXED_GAS *gas_dat;							///< Pointer to the MIXED_GAS data structure (may or may not be used)
	MAGPIE_DATA magpie_dat;						///< Data structure for a magpie problem (to be used if not using skua)
	std::vector<FINCH_DATA> finch_dat;			///< Data structure for pore adsorption kinetics for all species (u in mol/L)
	std::vector<SCOPSOWL_PARAM_DATA> param_dat;	///< Data structure for parameter info for all species
	
	std::vector<SKUA_DATA> skua_dat;	///< Data structure holding a skua object for all nodes (each skua has an object for each species)
	
}SCOPSOWL_DATA;

/// Function to print out the main header for the output file
void print2file_species_header(FILE *Output, SCOPSOWL_DATA *owl_dat, int i);

/// Function to print out the time and space header for the output file
void print2file_SCOPSOWL_time_header(FILE *Output, SCOPSOWL_DATA *owl_dat, int i);

/// Function to call the species and time header functions
void print2file_SCOPSOWL_header(SCOPSOWL_DATA *owl_dat);

/// Function to print out the old time results to the output file
void print2file_SCOPSOWL_result_old(SCOPSOWL_DATA *owl_dat);

/// Function to print out the new time results to the output file
void print2file_SCOPSOWL_result_new(SCOPSOWL_DATA *owl_dat);

/// Default function for evaluating adsorption and adsorption strength
/** This function is called in the preprocesses and postprocesses to estimate the strength
	of adsorption in the macro-scale problem from perturbations. It will use perturbations
	in either the MAGPIE simulation or SKUA simulation, depending on the type of problem
	the user is solving.
 
	\param i index for the ith species in the system
	\param l index for the lth node in the macro-scale domain
	\param user_data pointer for the SCOSPOWL_DATA structure*/
double default_adsorption(int i, int l, const void *user_data);

/// Default function for evaluating retardation coefficient
/** This function is called in the preprocesses and postprocesses to estimate the retardation
	coefficient for the simulation. It is recalculated at every time step to keep track of
	all changing conditions in the simulation.
 
	\param i index for the ith species in the system
	\param l index for the lth node in the macro-scale domain
	\param user_data pointer for the SCOSPOWL_DATA structure*/
double default_retardation(int i, int l, const void *user_data);

/// Default function for evaluating pore diffusivity
/** This function is called during the evaluation of non-linear residuals to more accurately
	represent non-linearities in the pore diffusion behavior. The pore diffusion is calculated
	based on kinetic theory of gases (see egret.h) and is adjusted according to the Knudsen
	Diffusion model and the porosity of the binder material.
 
	\param i index for the ith species in the system
	\param l index for the lth node in the macro-scale domain
	\param user_data pointer for the SCOSPOWL_DATA structure*/
double default_pore_diffusion(int i, int l, const void *user_data);

/// Default function for evaluating surface diffusion for HOMOGENEOUS pellets
/** This function is ONLY used if the pellet is determined to be homogeneous. Otherwise,
	this is replaced by the surface diffusion function for the SKUA simulation. The diffusivity
	is calculated based on the Arrhenius rate expression and then adjusted by the outside partial
	pressure of the adsorbing species.
 
	\param i index for the ith species in the system
	\param l index for the lth node in the macro-scale domain
	\param user_data pointer for the SCOSPOWL_DATA structure*/
double default_surf_diffusion(int i, int l, const void *user_data);

/// Default function for evaluating effective diffusivity for HOMOGENEOUS pellets
/** This function is ONLY used if the pellet is determined to be homogeneous. Otherwise,
	this is replaced by the pore diffusion function. The effective diffusivity is determined
	by the combination of pore diffusivity and surface diffusivity with adsorption strength
	in an homogeneous pellet.
 
	\param i index for the ith species in the system
	\param l index for the lth node in the macro-scale domain
	\param user_data pointer for the SCOSPOWL_DATA structure*/
double default_effective_diffusion(int i, int l, const void *user_data);

/// Constant pore diffusion function for homogeneous or heterogeneous pellets
/** This function should be used if the user wants to specify a constant pore diffusivity.
	The value of pore diffusion is then set equal to the value of pore_diffusion in the
	SCOPSOWL_PARAM_DATA structure.
 
	\param i index for the ith species in the system
	\param l index for the lth node in the macro-scale domain
	\param user_data pointer for the SCOSPOWL_DATA structure*/
double const_pore_diffusion(int i, int l, const void *user_data);

/// Default function for evaluating the film mass transfer coefficient
/** This function is called during the setup of the boundary conditions and is used to
	estimate the film mass transfer coefficient for the macro-scale problem. The coefficient
	is calculated according to the kinetic theory of gases (see egret.h).
 
	\param i index for the ith species in the system
	\param user_data pointer for the SCOSPOWL_DATA structure*/
double default_filmMassTransfer(int i, const void *user_data);

/// Constant film mass transfer coefficient function
/** This function is used when the user wants to specify a constant value for film
	mass transfer. The value of that coefficient is then set equal to the value of
	film_transfer in the SCOPSOWL_PARAM_DATA structure.
 
	\param i index for the ith species in the system
	\param user_data pointer for the SCOSPOWL_DATA structure*/
double const_filmMassTransfer(int i, const void *user_data);

/// Setup function to allocate memory and setup function pointers for the SCOPSOWL simulation
/** This function sets up the memory and function pointers used in SCOPSOWL simulations. User can
	provide NULL in place of functions for the function pointers and the setup will automatically
	use just the default settings. However, the user is required to pass the necessary data structure
	pointers for MIXED_GAS and SCOPSOWL_DATA.
 
	\param file pointer to the output file to print out results
	\param eval_sorption pointer to the adsorption evaluation function
	\param eval_retardation pointer to the retardation evaluation function
	\param eval_pore_diff pointer to the pore diffusion function
	\param eval_filmMT pointer to the film mass transfer function
	\param eval_surface_diff pointer to the surface diffusion function (required)
	\param user_data pointer to the user's data structure used for the parameter functions
	\param gas_data pointer to the MIXED_GAS structure used to evaluate kinetic gas theory
	\param owl_data pointer to the SCOPSOWL_DATA structure*/
int setup_SCOPSOWL_DATA(FILE *file,
						double (*eval_sorption) (int i, int l, const void *user_data),
						double (*eval_retardation) (int i, int l, const void *user_data),
						double (*eval_pore_diff) (int i, int l, const void *user_data),
						double (*eval_filmMT) (int i, const void *user_data),
						double (*eval_surface_diff) (int i, int l, const void *user_data),
						const void *user_data,MIXED_GAS *gas_data,SCOPSOWL_DATA *owl_data);

/// SCOPSOWL executioner function to solve a time step
/** This function will call the preprocess, solver, and postprocess functions to evaluate
	a single time step in a simulation. All simulation conditions must be set prior to
	calling this function. This function will typically be the one called from other
	simulations that will involve a SCOPSOWL evaluation to resolve kinetic coupling.
 
	\param owl_dat pointer to the SCOPSOWL_DATA structure (must be initialized)*/
int SCOPSOWL_Executioner(SCOPSOWL_DATA *owl_dat);

/// Function to set the initial conditions for a SCOPSOWL simulation
/** This function will setup the initial conditions of the simulation based on the initial
	temperature, pressure, and adsorbed molefractions. It assumes that the initial conditions
	are constant throughout the domain of the problem. This function should only be called
	once during a simulation.
 
	\param owl_dat pointer to the SCOPSOWL_DATA structure (must be initialized)*/
int set_SCOPSOWL_ICs(SCOPSOWL_DATA *owl_dat);


/// Function to set the timestep of the SCOPSOWL simulation
/** This function is used to set the next time step to be used in the SCOPSOWL simulation.
	A constant time step based on the size of the pellet discretization will be used. Users
	may want to use a custom time step to ensure that coupled-multi-scale systems are all
	in sync.
 
	\param owl_dat pointer to the SCOPSOWL_DATA structure (must be initialized)*/
int set_SCOPSOWL_timestep(SCOPSOWL_DATA *owl_dat);

/// Function to perform all preprocess SCOPSOWL operations
/** This function will update the boundary conditions and simulation conditions based on
	the current temperature, pressure, and gas phase composition, which may all vary in
	time. Since this function is called by the SCOPSOWL_Executioner, it does not need
	to be called explicitly by the user.
 
	\param owl_dat pointer to the SCOPSOWL_DATA structure (must be initialized)*/
int SCOPSOWL_preprocesses(SCOPSOWL_DATA *owl_dat);

/// Function to set the values of all non-linear system parameters during simulation
/** This is the function override for the FINCH setparams function (see finch.h). It will
	update the values of non-linear parameters in the residuals so that all variables in
	a species' system are fully coupled.
 
	\param user_data pointer to the SCOPSOWL_DATA structure (must be initialized)*/
int set_SCOPSOWL_params(const void *user_data);

/// Function to perform all postprocess SCOPSOWL operations
/** This function will update the retardation coefficients based on newly obtained
	simulation results for the current time step and calculate the average and total
	amount of adsorption of each species in the domain. Additionally, this function
	will call the print functions to store simulation results in the output file.
 
	\param owl_dat pointer to the SCOPSOWL_DATA structure (must be initialized)*/
int SCOPSOWL_postprocesses(SCOPSOWL_DATA *owl_dat);

/// Function to reset all stateful information to prepare for next simulation
/** This function will update the stateful information used in SCOPSOWL to prepare the
	system for the next time step in the simulation. However, because updating the
	states erases the old state, the user must be absolutely sure that the simulation
	is ready to be updated. For just running standard simulations, this is not an
	issue, but in coupling with other simulations it is very important.
 
	\param owl_dat pointer to the SCOPSOWL_DATA structure (must be initialized)*/
int SCOPSOWL_reset(SCOPSOWL_DATA *owl_dat);

/// Function to progress the SCOPSOWL simulation through time till complete
/** This function will call the initial conditions, then progressively call the
	executioner, time step, and reset functions to propagate the simulation in
	time. As such, this function is primarily used when running a SCOPSOWL
	simulation by itself and not when coupling it to an other problem.
 
	\param owl_dat pointer to the SCOPSOWL_DATA structure (must be initialized)*/
int SCOPSOWL(SCOPSOWL_DATA *owl_dat);

#endif
