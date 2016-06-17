/*!
 *  \file skua.h skua.cpp
 *	\brief Surface Kinetics for Uptake by Adsorption
 *  \author Austin Ladshaw
 *	\date 01/26/2015
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science. Copyright (c) 2015, all
 *             rights reserved.
 */

#include "skua.h"

//Function to print out the species header
void print2file_species_header(FILE *Output, SKUA_DATA *skua_dat, int i)
{
	const char *name = skua_dat->param_dat[i].speciesName.c_str();
	fprintf(Output, "%s\t", name);
	if (skua_dat->finch_dat[i].Dirichlet == true)
		fprintf(Output,"-\t");
	for (int l=0; l<skua_dat->finch_dat[i].LN+4; l++)
	{
		fprintf(Output,"-\t");
	}
}

//Print time header for skua output
void print2file_SKUA_time_header(FILE *Output, SKUA_DATA *skua_dat, int i)
{
	fprintf(Output,"Time\t");
  	fprintf(Output,"q[0] (mol/kg)\t");
  	for (int l=0; l<skua_dat->finch_dat[i].LN; l++)
  	{
      	if (skua_dat->finch_dat[i].Dirichlet == false && l == skua_dat->finch_dat[i].LN-1)
          	break;
      	fprintf(Output,"q[%i] (mol/kg)\t",l+1);
  	}
  	fprintf(Output,"q Total (mol/kg)\tq Average (mol/kg)\t");
	fprintf(Output,"x[%i] (-)\tQst[%i] (J/mol)\t",i,i);
}

//Print out the skua header to file
void print2file_SKUA_header(SKUA_DATA *skua_dat)
{
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		if (skua_dat->param_dat[i].Adsorbable == true)
		{
			print2file_species_header(skua_dat->OutputFile, skua_dat, i);
			print2file_tab(skua_dat->OutputFile, &skua_dat->finch_dat[i]);
		}
	}
	print2file_newline(skua_dat->OutputFile, NULL);
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		if (skua_dat->param_dat[i].Adsorbable == true)
		{
			print2file_dim_header(skua_dat->OutputFile, &skua_dat->finch_dat[i]);
			print2file_tab(skua_dat->OutputFile, &skua_dat->finch_dat[i]);
			print2file_tab(skua_dat->OutputFile, &skua_dat->finch_dat[i]);
			print2file_tab(skua_dat->OutputFile, &skua_dat->finch_dat[i]);
			print2file_tab(skua_dat->OutputFile, &skua_dat->finch_dat[i]);
			print2file_tab(skua_dat->OutputFile, &skua_dat->finch_dat[i]);
		}
	}
	print2file_newline(skua_dat->OutputFile, NULL);
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		if (skua_dat->param_dat[i].Adsorbable == true)
		{
			print2file_SKUA_time_header(skua_dat->OutputFile, skua_dat, i);
			print2file_tab(skua_dat->OutputFile, &skua_dat->finch_dat[i]);
		}
	}
	print2file_newline(skua_dat->OutputFile, NULL);
}

//Print old results to file
void print2file_SKUA_results_old(SKUA_DATA *skua_dat)
{
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		if (skua_dat->param_dat[i].Adsorbable == true)
		{
			print2file_result_old(skua_dat->OutputFile, &skua_dat->finch_dat[i]);
			fprintf(skua_dat->OutputFile, "%.6g\t%.6g\t",skua_dat->param_dat[i].xn,skua_dat->param_dat[i].Qstn);
			print2file_tab(skua_dat->OutputFile, &skua_dat->finch_dat[i]);
		}
	}
	print2file_newline(skua_dat->OutputFile, NULL);
}

//Print new results to file
void print2file_SKUA_results_new(SKUA_DATA *skua_dat)
{
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		if (skua_dat->param_dat[i].Adsorbable == true)
		{
			print2file_result_new(skua_dat->OutputFile, &skua_dat->finch_dat[i]);
			fprintf(skua_dat->OutputFile, "%.6g\t%.6g\t",skua_dat->param_dat[i].xnp1,skua_dat->param_dat[i].Qstnp1);
			print2file_tab(skua_dat->OutputFile, &skua_dat->finch_dat[i]);
		}
	}
	print2file_newline(skua_dat->OutputFile, NULL);
}

//Default function for computing Dc
double default_Dc(int i, int l, const void *data)
{
	double Dc = 0.0;
	SKUA_DATA *dat = (SKUA_DATA *) data;
	
	if (i < 0)
	{
		mError(indexing_error);
		return -1;
	}
	Dc = D_o(dat->param_dat[i].ref_diffusion, dat->param_dat[i].activation_energy, dat->magpie_dat.sys_dat.T);
	Dc = D_inf(Dc, dat->param_dat[i].ref_temperature, dat->param_dat[i].affinity, dat->param_dat[i].ref_pressure, dat->magpie_dat.sys_dat.T);
	
	return Dc;
}

//Default function for computing kf (Film Mass Transfer: um/hr)
double default_kf(int i, const void *data)
{
	return 0.0;
}

//Function to set Dc to a constant based on ref_diffusion
double const_Dc(int i, int l, const void *data)
{
	SKUA_DATA *dat = (SKUA_DATA *) data;
	return dat->param_dat[i].ref_diffusion;
}

//Simple Darken Dc
double simple_darken_Dc(int i, int l, const void *data)
{
	double Dc = 0.0;
	SKUA_DATA *dat = (SKUA_DATA *) data;
	
	if (i < 0)
	{
		mError(indexing_error);
		return -1;
	}
	
	//Boundary Dc
	if (l < 0)
	{
		Dc = D_c(D_o(dat->param_dat[i].ref_diffusion, dat->param_dat[i].activation_energy, dat->magpie_dat.sys_dat.T),dat->finch_dat[i].uo/(dat->magpie_dat.gsta_dat[i].qmax+1.0E-6));
	}
	//Interior Dc
	else
	{
		Dc = D_c(D_o(dat->param_dat[i].ref_diffusion, dat->param_dat[i].activation_energy, dat->magpie_dat.sys_dat.T),dat->finch_dat[i].unp1(l,0)/(dat->magpie_dat.gsta_dat[i].qmax+1.0E-6));
	}
	Dc = D_inf(Dc, dat->param_dat[i].ref_temperature, dat->param_dat[i].affinity, dat->param_dat[i].ref_pressure, dat->magpie_dat.sys_dat.T);
	
	return Dc;
}

//Average value Darken Dc {DO NOT USE}
/*
 double average_darken_Dc(int i, int l, const void *data)
 {
 double Dc = 0.0;
 SKUA_DATA *dat = (SKUA_DATA *) data;
 
 if (i < 0)
 {
 mError(indexing_error);
 return -1;
 }
 
 //Boundary Dc
 if (l < 0)
 {
 Dc = D_c(D_o(dat->param_dat[i].ref_diffusion, dat->param_dat[i].activation_energy, dat->magpie_dat.sys_dat.T),dat->finch_dat[i].uAvg/(dat->magpie_dat.gsta_dat[i].qmax+1.0E-6));
 }
 //Interior Dc
 else
 {
 Dc = D_c(D_o(dat->param_dat[i].ref_diffusion, dat->param_dat[i].activation_energy, dat->magpie_dat.sys_dat.T),dat->finch_dat[i].uAvg/(dat->magpie_dat.gsta_dat[i].qmax+1.0E-6));
 }
 
 return Dc;
 }
 */

//Theoretical Darken Dc (Much more computationally intense)
double theoretical_darken_Dc(int i, int l, const void *data)
{
	double Dc = 0.0;
	SKUA_DATA *dat = (SKUA_DATA *) data;
	
	if (i < 0)
	{
		mError(indexing_error);
		return -1;
	}
	
	if (dat->param_dat[i].Adsorbable == true)
	{
		//Boundary Dc
		if (l < 0)
		{
			Dc = D_o(dat->param_dat[i].ref_diffusion, dat->param_dat[i].activation_energy, dat->magpie_dat.sys_dat.T);
			Dc = Dc * q_p(dat->y[i]*dat->magpie_dat.sys_dat.PT, (void *)&dat->magpie_dat, i) / dq_dp(dat->y[i]*dat->magpie_dat.sys_dat.PT, (void *)&dat->magpie_dat, i);
		}
		//Interior Dc
		else
		{
			
			dat->magpie_dat.sys_dat.Recover = true;
			dat->magpie_dat.sys_dat.Carrier = true;
			dat->magpie_dat.gpast_dat[i].x = 1.0;
			dat->magpie_dat.sys_dat.qT = dat->finch_dat[i].unp1(l,0);
			for (int j=0; j<dat->magpie_dat.sys_dat.N; j++)
			{
				if (j!=i)
					dat->magpie_dat.gpast_dat[j].x = 0.0;
			}
			
			//Call MAGPIE to determine equilibrium adsorption
			double success = 0;
			if (dat->magpie_dat.sys_dat.qT > 0.0)
			{
				success = MAGPIE((void *)&dat->magpie_dat);
				if (success < 0 || success > 5) {mError(simulation_fail); return -1;}
				else success = 0;
			}
			else
			{
				dat->magpie_dat.gpast_dat[i].y = 0.0;
			}
			dat->magpie_dat.sys_dat.Recover = false;
			dat->magpie_dat.sys_dat.Carrier = false;
			dat->total_steps = dat->total_steps + dat->magpie_dat.sys_dat.total_eval;
			
			
			Dc = D_o(dat->param_dat[i].ref_diffusion, dat->param_dat[i].activation_energy, dat->magpie_dat.sys_dat.T);
			Dc = Dc * q_p(dat->magpie_dat.gpast_dat[i].y*dat->magpie_dat.sys_dat.PT, (void *)&dat->magpie_dat, i) / dq_dp(dat->magpie_dat.gpast_dat[i].y*dat->magpie_dat.sys_dat.PT, (void *)&dat->magpie_dat, i);
		}
	}
	else
	{
		Dc = 0.0;
	}
	Dc = D_inf(Dc, dat->param_dat[i].ref_temperature, dat->param_dat[i].affinity, dat->param_dat[i].ref_pressure, dat->magpie_dat.sys_dat.T);
	
	return Dc;
}

//Based on average qT and y_eff to use the theoretical Darken Dc expression {DO NOT USE}
/*
 double average_theoretical_darken_Dc(int i, int l, const void *data)
 {
 double Dc = 0.0;
 SKUA_DATA *dat = (SKUA_DATA *) data;
 
 if (i < 0)
 {
 mError(indexing_error);
 return -1;
 }
 
 if (dat->param_dat[i].Adsorbable == true)
 {
 //Boundary Dc
 if (l < 0)
 {
 Dc = D_o(dat->param_dat[i].ref_diffusion, dat->param_dat[i].activation_energy, dat->magpie_dat.sys_dat.T);
 Dc = Dc * q_p(dat->param_dat[i].y_eff*dat->magpie_dat.sys_dat.PT, (void *)&dat->magpie_dat, i) / dq_dp(dat->param_dat[i].y_eff*dat->magpie_dat.sys_dat.PT, (void *)&dat->magpie_dat, i);
 }
 //Interior Dc
 else
 {
 Dc = D_o(dat->param_dat[i].ref_diffusion, dat->param_dat[i].activation_energy, dat->magpie_dat.sys_dat.T);
 Dc = Dc * q_p(dat->param_dat[i].y_eff*dat->magpie_dat.sys_dat.PT, (void *)&dat->magpie_dat, i) / dq_dp(dat->param_dat[i].y_eff*dat->magpie_dat.sys_dat.PT, (void *)&dat->magpie_dat, i);
 }
 }
 else
 {
 Dc = 0.0;
 }
 
 return Dc;
 }
 */

//Estimate crystal diffusivity with a constant Darken model {DO NOT USE}
/*
 double const_darken_Dc(int i, int l, const void *data)
 {
 return theoretical_darken_Dc(i, -1, data);
 }
 */

//Estimation of Film Mass Transfer based on gas phase properties
double empirical_kf(int i, const void *data)
{
	int success = 0;
	SKUA_DATA *dat = (SKUA_DATA *) data;
	if (dat->gas_dat != NULL)
	{
		success = set_variables(dat->magpie_dat.sys_dat.PT, dat->magpie_dat.sys_dat.T, dat->gas_velocity, dat->pellet_radius*2.0E-4, dat->y, dat->gas_dat);
		if (success != 0) {mError(simulation_fail); return -1;}
		//Calculate Properties
		success = calculate_properties(dat->gas_dat);
		if (success != 0) {mError(simulation_fail); return -1;}
		return FilmMTCoeff(dat->gas_dat->species_dat[i].molecular_diffusion, dat->gas_dat->char_length, dat->gas_dat->Reynolds, dat->gas_dat->species_dat[i].Schmidt) * 3.6E+7;
	}
	else
	{
		std::cout << "\nWARNING!!! Calling a Function without the correct pointer initialized!\n" << std::endl;
		return default_kf(i, data);
	}
}

//Function to give the constant film mass transfer rate {NOT USED}
double const_kf(int i, const void *data)
{
	SKUA_DATA *dat = (SKUA_DATA *) data;
	return dat->param_dat[i].film_transfer;
}

//Check the input mole fractions from files in SKUA for errors
int molefractionCheck(SKUA_DATA *skua_dat)
{
	int success = 0;
	double sumx = 0.0, sumy = 0.0;
	
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		sumx = sumx + skua_dat->param_dat[i].xIC;
		sumy = sumy + skua_dat->y[i];
	}
	if (sumy > (1.0 + 1.0E-6) || sumx > (1.0 + 1.0E-6))
	{
		mError(invalid_molefraction);
		success = -1;
	}
	
	return success;
}

//Function to setup remaining memory and structures after initializing parameter information
int setup_SKUA_DATA(FILE *file, double (*eval_Dc) (int i, int l, const void *user_data),
					double (*eval_Kf) (int i, const void *user_data), const void *user_data,
					MIXED_GAS *gas_data, SKUA_DATA *skua_dat)
{
	int success = 0;
	
	if (file == NULL)
		skua_dat->Print2File = false;
	else
	{
		skua_dat->Print2File = true;
		skua_dat->OutputFile = file;
	}
	if ((*eval_Dc) == NULL)
		skua_dat->eval_diff = (*default_Dc);
	else
		skua_dat->eval_diff = (*eval_Dc);
	if ((*eval_Kf) == NULL)
		skua_dat->eval_kf = (*default_kf);
	else
		skua_dat->eval_kf = (*eval_Kf);
	skua_dat->user_data = user_data;
	if (gas_data == NULL)
		skua_dat->DirichletBC = true;
	else
		skua_dat->gas_dat = gas_data;
	
	//Loop to setup the FINCH data
	skua_dat->t = 0.0;
	skua_dat->t_old	= 0.0;
	skua_dat->t_counter = 0.0;
	skua_dat->magpie_dat.sys_dat.Output = false;
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		skua_dat->finch_dat[i].t = skua_dat->t;
		skua_dat->finch_dat[i].t_old = skua_dat->t_old;
		skua_dat->finch_dat[i].L = skua_dat->pellet_radius;			//um
		skua_dat->finch_dat[i].T = skua_dat->sim_time;				//hrs
		skua_dat->finch_dat[i].Dirichlet = skua_dat->DirichletBC;	//-
		skua_dat->finch_dat[i].LN = 10;						//-
		skua_dat->finch_dat[i].SteadyState = false;			//-
		skua_dat->finch_dat[i].CheckMass = true;				//-
		skua_dat->finch_dat[i].Iterative = skua_dat->NonLinear;		//-
		skua_dat->finch_dat[i].NormTrack = skua_dat->Print2Console;	//-
		skua_dat->finch_dat[i].nl_method = FINCH_Picard;					//-
		skua_dat->finch_dat[i].d = skua_dat->coord;							//-
		skua_dat->finch_dat[i].s = skua_dat->char_measure;						//-
		skua_dat->finch_dat[i].t = 0.0;						//hrs
		skua_dat->finch_dat[i].t_old = 0.0;					//hrs
		skua_dat->finch_dat[i].Ro = 1.0;						//-
		skua_dat->finch_dat[i].RIC = 1.0;						//-
		skua_dat->finch_dat[i].vo = 0.0;						//um/hr
		skua_dat->finch_dat[i].vIC = 0.0;						//um/hr
		skua_dat->finch_dat[i].kIC = 0.0;						//1/hr
		skua_dat->finch_dat[i].ko = 0.0;						//1/hr
		skua_dat->param_dat[i].film_transfer = 0.0;
		
		success = setup_FINCH_DATA(NULL,NULL,NULL,NULL,NULL,set_SKUA_params,NULL,NULL,NULL,NULL,NULL,NULL,&skua_dat->finch_dat[i],user_data);
		if (success != 0) {mError(simulation_fail); return -1;}
	}
	
	return success;
}
//Execution order for the SKUA routine (Note: the time step must be set prior to calling this routine)
int SKUA_Executioner(SKUA_DATA *skua_dat)
{
	int success = 0;
	
	//Perform Preprocess Actions
	success = SKUA_preprocesses(skua_dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Loop for all components
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		
		if (skua_dat->param_dat[i].Adsorbable == true)
		{
			//Solve the system
			success = (*skua_dat->finch_dat[i].solve) ((void *)&skua_dat->finch_dat[i]);
			if (success != 0) {mError(simulation_fail); return -1;}
			
			//Check for negative concentrations
			success = check_Mass(&skua_dat->finch_dat[i]);
			if (success != 0) {mError(simulation_fail); return -1;}
			
			//Form Totals
			success = uTotal(&skua_dat->finch_dat[i]);
			if (success != 0) {mError(simulation_fail); return -1;}
			
			//Form Averages
			success = uAverage(&skua_dat->finch_dat[i]);
			if (success != 0) {mError(simulation_fail); return -1;}
		}
		
	}
	
	//Perform Postprocess Actions
	success = SKUA_postprocesses(skua_dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	return success;
}

//Set the initial conditions
int set_SKUA_ICs(SKUA_DATA *skua_dat)
{
	int success = 0;
	
	//Set MAGPIE parameters to perform reverse evaluations
	skua_dat->magpie_dat.sys_dat.Recover = true;
	skua_dat->magpie_dat.sys_dat.Carrier = false;
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		if (skua_dat->param_dat[i].Adsorbable == false)
		{
			skua_dat->magpie_dat.sys_dat.Carrier = true;
			break;
		}
	}
	
	//Set the initial conditions
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		skua_dat->param_dat[i].ref_pressure = skua_dat->y[i]*skua_dat->magpie_dat.sys_dat.PT;
		if (skua_dat->param_dat[i].Adsorbable == true)
			skua_dat->magpie_dat.gpast_dat[i].x = skua_dat->param_dat[i].xIC;
		else
			skua_dat->magpie_dat.gpast_dat[i].x = 0.0;
		skua_dat->param_dat[i].xn = skua_dat->magpie_dat.gpast_dat[i].x;
		skua_dat->param_dat[i].xnp1 = skua_dat->param_dat[i].xn;
		skua_dat->finch_dat[i].uAvg = skua_dat->magpie_dat.sys_dat.qT * skua_dat->magpie_dat.gpast_dat[i].x;
		
	}
	
	//Call MAGPIE to determine equilibrium adsorption
	success = MAGPIE((void *)&skua_dat->magpie_dat);
	if (success < 0 || success > 5) {mError(simulation_fail); return -1;}
	else success = 0;
	skua_dat->total_steps = skua_dat->total_steps + skua_dat->magpie_dat.sys_dat.total_eval;
	
	//Use MAGPIE result to form initial condition in adsorbent
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		skua_dat->finch_dat[i].Sn.ConstantICFill(0.0);
		skua_dat->finch_dat[i].Snp1.ConstantICFill(0.0);
		skua_dat->finch_dat[i].vn.ConstantICFill(0.0);
		skua_dat->finch_dat[i].vnp1.ConstantICFill(0.0);
		skua_dat->finch_dat[i].kn.ConstantICFill(0.0);
		skua_dat->finch_dat[i].knp1.ConstantICFill(0.0);
		skua_dat->finch_dat[i].Rn.ConstantICFill(1.0);
		skua_dat->finch_dat[i].Rnp1.ConstantICFill(1.0);
		skua_dat->finch_dat[i].un.ConstantICFill(skua_dat->magpie_dat.sys_dat.qT*skua_dat->param_dat[i].xIC);
		skua_dat->finch_dat[i].unm1 = skua_dat->finch_dat[i].un;
		skua_dat->finch_dat[i].unp1 = skua_dat->finch_dat[i].un;
		skua_dat->param_dat[i].y_eff = skua_dat->magpie_dat.gpast_dat[i].y;
		if (skua_dat->param_dat[i].Adsorbable == true)
		{
			skua_dat->param_dat[i].Qstn = Qst((skua_dat->magpie_dat.gpast_dat[i].y*skua_dat->magpie_dat.sys_dat.PT),(void *)&skua_dat->magpie_dat,i);
		}
		else
		{
			skua_dat->param_dat[i].Qstn = 0.0;
		}
		skua_dat->param_dat[i].Qstnp1 = skua_dat->param_dat[i].Qstn;
		if (skua_dat->param_dat[i].Adsorbable == true)
		{
			for (int l=0; l<skua_dat->finch_dat[i].LN; l++)
			{
				skua_dat->finch_dat[i].Dn(l,0) = (*skua_dat->eval_diff) (i,l,skua_dat->user_data);
				skua_dat->finch_dat[i].Dnp1(l,0) = skua_dat->finch_dat[i].Dn(l,0);
			}
		}
		else
		{
			skua_dat->finch_dat[i].Dn.ConstantICFill(0.0);
			skua_dat->finch_dat[i].Dnp1.ConstantICFill(0.0);
		}
		skua_dat->finch_dat[i].Do = skua_dat->finch_dat[i].Dn(skua_dat->finch_dat[i].LN-1,0);
		skua_dat->finch_dat[i].uo = skua_dat->finch_dat[i].un(skua_dat->finch_dat[i].LN-1,0);
		skua_dat->finch_dat[i].kfn = (*skua_dat->eval_kf) (i,skua_dat->user_data);
		if (skua_dat->finch_dat[i].kfn < 0.0)
		{
			mError(simulation_fail);
			return -1;
		}
		
		//Form Totals
		success = uTotal(&skua_dat->finch_dat[i]);
		if (success != 0) {mError(simulation_fail); return -1;}
		
		//Form Averages
		success = uAverage(&skua_dat->finch_dat[i]);
		if (success != 0) {mError(simulation_fail); return -1;}
		
		//Update averages and totals
		skua_dat->finch_dat[i].uT_old = skua_dat->finch_dat[i].uT;
		skua_dat->finch_dat[i].uAvg_old = skua_dat->finch_dat[i].uAvg;
	}
	skua_dat->qTnp1 = skua_dat->qTn;
	
	//Reset the MAGPIE parameters before continuing
	skua_dat->magpie_dat.sys_dat.Recover = false;
	skua_dat->magpie_dat.sys_dat.Carrier = false;
	
	return success;
}

//Setup the time step
int set_SKUA_timestep(SKUA_DATA *skua_dat)
{
	int success = 0;
	
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		skua_dat->finch_dat[i].dt = skua_dat->finch_dat[i].dz / 4.0;
		skua_dat->finch_dat[i].t = skua_dat->finch_dat[i].dt + skua_dat->finch_dat[i].t_old;
	}
	skua_dat->t_old = skua_dat->finch_dat[0].t_old;
	skua_dat->t = skua_dat->finch_dat[0].t;
	
	return success;
}

//Perform Preprocess Actions
int SKUA_preprocesses(SKUA_DATA *skua_dat)
{
	int success = 0;
	
	//First, establish the boundary conditions by solving a MAGPIE problem
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		double p_out = skua_dat->y[i]*skua_dat->magpie_dat.sys_dat.PT;
		if (p_out >= skua_dat->param_dat[i].ref_pressure) skua_dat->param_dat[i].ref_pressure = p_out;
		if (skua_dat->param_dat[i].Adsorbable == true)
		{
			skua_dat->magpie_dat.gpast_dat[i].y = skua_dat->y[i];
		}
		else
		{
			skua_dat->magpie_dat.gpast_dat[i].y = 0.0;
		}
	}
	
	//Call MAGPIE to determine equilibrium adsorption
	success = MAGPIE((void *)&skua_dat->magpie_dat);
	if (success < 0 || success > 3) {mError(simulation_fail); return -1;}
	else success = 0;
	skua_dat->total_steps = skua_dat->total_steps + skua_dat->magpie_dat.sys_dat.total_eval;
	
	//Loop to set bounary values and parameters
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		//Bounary value is adsorbed amount
		if (skua_dat->param_dat[i].Adsorbable == true)
		{
			skua_dat->finch_dat[i].uo = skua_dat->magpie_dat.gpast_dat[i].q;
			skua_dat->finch_dat[i].kfnp1 = (*skua_dat->eval_kf) (i,skua_dat->user_data);
			if (skua_dat->finch_dat[i].kfnp1 < 0.0)
			{
				mError(simulation_fail);
				return -1;
			}
		}
		else
		{
			skua_dat->finch_dat[i].uo = 0.0;
			skua_dat->finch_dat[i].kfnp1 = 0.0;
		}
	}
	
	return success;
}

//Set the parameters for the FINCH problem
int set_SKUA_params(const void *user_data)
{
	int success = 0;
	
	//In this file, we force an Implicit Only FINCH solution
	SKUA_DATA *dat = (SKUA_DATA *) user_data;
	for (int i=0; i<dat->magpie_dat.sys_dat.N; i++)
	{
		dat->finch_dat[i].CN = false;
		dat->finch_dat[i].beta = 1.0;
		if (dat->param_dat[i].Adsorbable == true)
		{
			dat->finch_dat[i].Do = (*dat->eval_diff) (i,-1,dat->user_data);
			//Loop over the nodes
			for (int l=0; l<dat->finch_dat[i].LN; l++)
			{
				dat->finch_dat[i].Dnp1(l,0) = (*dat->eval_diff) (i,l,dat->user_data);
			}
		}
		else
		{
			dat->finch_dat[i].Do = 0.0;
			dat->finch_dat[i].Dnp1.ConstantICFill(0.0);
		}
	}
	
	return success;
}

//Perform Postprocess Actions
int SKUA_postprocesses(SKUA_DATA *skua_dat)
{
	int success = 0;
	
	//Update some SKUA information
	skua_dat->qTnp1 = 0.0;
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		skua_dat->qTnp1 = skua_dat->qTnp1 + skua_dat->finch_dat[i].uAvg;
		skua_dat->total_steps = skua_dat->total_steps + skua_dat->finch_dat[i].total_iter;
		skua_dat->finch_dat[i].total_iter = 0;
	}
	if (skua_dat->qTnp1 >=0.0)
		skua_dat->magpie_dat.sys_dat.qT = skua_dat->qTnp1;
	else
		skua_dat->magpie_dat.sys_dat.qT = 0.0;
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		if (skua_dat->qTnp1 > 0.0)
			skua_dat->param_dat[i].xnp1 = skua_dat->finch_dat[i].uAvg / skua_dat->qTnp1;
		else
			skua_dat->param_dat[i].xnp1 = 0.0;
		skua_dat->magpie_dat.gpast_dat[i].x = skua_dat->param_dat[i].xnp1;
	}
	
	//Setup for another reverse MAGPIE evaluation
	skua_dat->magpie_dat.sys_dat.Recover = true;
	skua_dat->magpie_dat.sys_dat.Carrier = false;
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		if (skua_dat->param_dat[i].Adsorbable == false)
		{
			skua_dat->magpie_dat.sys_dat.Carrier = true;
			break;
		}
	}
	
	//Call MAGPIE to determine equilibrium adsorption
	success = MAGPIE((void *)&skua_dat->magpie_dat);
	if (success < 0 || success > 5)
	{
		mError(simulation_fail);
		std::cout << success << std::endl;
		for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
		{
			std::cout << "x[" << i << "] = " << skua_dat->magpie_dat.gpast_dat[i].x << "\t";
			std::cout << "q[" << i << "] = " << skua_dat->magpie_dat.gpast_dat[i].q << "\t";
			std::cout << "y[" << i << "] = " << skua_dat->magpie_dat.gpast_dat[i].y << "\n";
		}
		return -1;
	}
	else success = 0;
	skua_dat->total_steps = skua_dat->total_steps + skua_dat->magpie_dat.sys_dat.total_eval;
	
	//Calculate the new heat of adsorption
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		if (skua_dat->param_dat[i].Adsorbable == true)
		{
			skua_dat->param_dat[i].Qstnp1 = Qst((skua_dat->magpie_dat.gpast_dat[i].y*skua_dat->magpie_dat.sys_dat.PT),(void *)&skua_dat->magpie_dat,i);
		}
		else
		{
			skua_dat->param_dat[i].Qstnp1 = 0.0;
		}
		skua_dat->param_dat[i].y_eff = skua_dat->magpie_dat.gpast_dat[i].y;
	}
	
	//Print results
	if (skua_dat->Print2File == true)
	{
		skua_dat->t_counter = skua_dat->t_counter + skua_dat->finch_dat[0].dt;
		if (skua_dat->t_counter > skua_dat->t_print)
		{
			print2file_SKUA_results_new(skua_dat);
			skua_dat->t_counter = 0.0;
		}
	}
	
	//Reset the MAGPIE parameters before continuing
	skua_dat->magpie_dat.sys_dat.Recover = false;
	skua_dat->magpie_dat.sys_dat.Carrier = false;
	
	return success;
}

//Reset the values for next time step
int SKUA_reset(SKUA_DATA *skua_dat)
{
	int success = 0;
	
	for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
	{
		success = (*skua_dat->finch_dat[i].resettime) ((void *)&skua_dat->finch_dat[i]);
		if (success != 0) {mError(simulation_fail); return -1;}
		skua_dat->param_dat[i].Qstn = skua_dat->param_dat[i].Qstnp1;
		skua_dat->param_dat[i].xn = skua_dat->param_dat[i].xnp1;
	}
	skua_dat->qTn = skua_dat->qTnp1;
	if (skua_dat->t_old == 0.0)
		skua_dat->t_old = 0.0;
	else
		skua_dat->t_old = skua_dat->t;
	
	return success;
}

//Calling the SKUA routine
int SKUA(SKUA_DATA *skua_dat)
{
	int success = 0;
	
	//Print to file
	if (skua_dat->Print2File == true)
	{
		print2file_SKUA_header(skua_dat);
	}
	
	//Set Initial Conditions
	success = set_SKUA_ICs(skua_dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Print ICs
	if (skua_dat->Print2File == true)
	{
		print2file_SKUA_results_old(skua_dat);
	}
	
	//Solve a series of time steps implicitly
	do
	{
		if (skua_dat->finch_dat[0].Update == true)
		{
			success = SKUA_reset(skua_dat);
			if (success != 0) {mError(simulation_fail); return -1;}
		}
		success = set_SKUA_timestep(skua_dat);
		if (success != 0) {mError(simulation_fail); return -1;}
		std::cout << "Evaluating time: " << skua_dat->t << " hrs..." << std::endl;
		success = SKUA_Executioner(skua_dat);
		if (success == 0)
		{
			std::cout << "Simulation Successful!\n" << std::endl;
			for (int i=0; i<skua_dat->magpie_dat.sys_dat.N; i++)
				skua_dat->finch_dat[i].Update = true;
		}
		else {mError(simulation_fail); skua_dat->finch_dat[0].Update = false; return -1;}
		skua_dat->total_steps++;
	} while (skua_dat->t < skua_dat->sim_time);
	
	return success;
}
