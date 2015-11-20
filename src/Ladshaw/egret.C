//----------------------------------------
//  Created by Austin Ladshaw on 1/29/15
//  Copyright (c) 2015
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

/*
 *			EGRET = Estimation of Gas-phase pRopErTies
 *
 *		This file is responsible for estimating various temperature, pressure, and concentration
 *		dependent parameters to be used in other models for gas phase adsorption, mass transfer,
 *		and or mass transport. The goal of this file is to eliminate redundancies in code such
 *		that the higher level programs operate more efficiently and cleanly. Calculations made
 *		here are based on kinetic theory of gases, ideal gas law, and some emperical models that
 *		were developed to account for changes in density and viscosity with changes in temperature
 *		between standard temperatures and up to 1000 K.
 */

#include "egret.h"

//Function to initialize the memory for gas data
int initialize_data(int N, MIXED_GAS *gas_dat)
{
	int success = 0;
	if (N <= 0)
	{
		mError(invalid_components);
		return -1;
	}
	gas_dat->N = N;
	gas_dat->species_dat.resize(N);
	gas_dat->binary_diffusion.set_size(N, N);
	gas_dat->molefraction.resize(N);
	return success;
}

//Function to set values of variables used in calculations
int set_variables(double PT, double T, double us, double L, std::vector<double> &y, MIXED_GAS *gas_dat)
{
	int success = 0;
	gas_dat->total_pressure = PT;
	gas_dat->gas_temperature = T;
	gas_dat->velocity = fabs(us);
	gas_dat->char_length = L;
	if (gas_dat->molefraction.size() != y.size())
		gas_dat->molefraction.resize(y.size());
	double ysum = 0.0;
	for (int i=0; i<y.size(); i++)
	{
		ysum = ysum + y[i];
		gas_dat->molefraction[i] = y[i];
		if (y[i] < -1.0E-6 && gas_dat->CheckMolefractions == true)
		{
			mError(invalid_molefraction);
			return -1;
		}
	}
	if ( (ysum > (1.0 + 1.0E-6) || ysum < (1.0 - 1.0E-6) ) && gas_dat->CheckMolefractions == true)
	{
		mError(invalid_gas_sum);
		return -1;
	}
	return success;
}

//Function to calculate the properties of the gas phase
int calculate_properties(MIXED_GAS *gas_dat)
{
	int success = 0;
	
	gas_dat->total_dyn_vis = 0.0;
	gas_dat->total_molecular_weight = 0.0;
	gas_dat->total_specific_heat = 0.0;
	
	//Check to see if there is only one gas species and quit early if so
	if (gas_dat->N <= 0)
	{
		mError(invalid_components);
		return -1;
	}
	else if (gas_dat->N == 1)
	{
		gas_dat->species_dat[0].density =  CE3(gas_dat->total_pressure*gas_dat->species_dat[0].molecular_weight,gas_dat->gas_temperature);
		gas_dat->species_dat[0].dynamic_viscosity = Mu(gas_dat->species_dat[0].Sutherland_Viscosity, gas_dat->species_dat[0].Sutherland_Temp, gas_dat->species_dat[0].Sutherland_Const, gas_dat->gas_temperature);
		gas_dat->binary_diffusion.edit(0, 0, D_ii(gas_dat->species_dat[0].density, gas_dat->species_dat[0].dynamic_viscosity));
		gas_dat->total_molecular_weight = gas_dat->species_dat[0].molecular_weight;
		gas_dat->total_specific_heat = gas_dat->species_dat[0].specific_heat;
		gas_dat->total_density = gas_dat->species_dat[0].density;
		gas_dat->species_dat[0].molecular_diffusion = gas_dat->binary_diffusion(0,0);
		gas_dat->total_dyn_vis = gas_dat->species_dat[0].dynamic_viscosity;
		gas_dat->kinematic_viscosity = Nu(gas_dat->total_dyn_vis,gas_dat->total_density);
		gas_dat->Reynolds = ReNum(gas_dat->velocity, gas_dat->char_length, gas_dat->kinematic_viscosity);
		gas_dat->species_dat[0].Schmidt = ScNum(gas_dat->kinematic_viscosity, gas_dat->species_dat[0].molecular_diffusion);
	}
	else
	{
		//Loop for all gas species to establish densities and viscosities
		for (int i=0; i<gas_dat->N; i++)
		{
			gas_dat->species_dat[i].density = CE3(gas_dat->total_pressure*gas_dat->species_dat[i].molecular_weight,gas_dat->gas_temperature);
			gas_dat->species_dat[i].dynamic_viscosity = Mu(gas_dat->species_dat[i].Sutherland_Viscosity, gas_dat->species_dat[i].Sutherland_Temp, gas_dat->species_dat[i].Sutherland_Const, gas_dat->gas_temperature);
			
			//Inner Loop for Binary Diffusion Tensor
			for (int j=0; j<=i; j++)
			{
				gas_dat->binary_diffusion.edit(i,j,1.0);
				if (j==i)
					gas_dat->binary_diffusion.edit(i,j,D_ii(gas_dat->species_dat[i].density, gas_dat->species_dat[i].dynamic_viscosity));
				else
				{
					//Calculate upper triangular portion
					gas_dat->binary_diffusion.edit(i, j, D_ij(gas_dat->species_dat[i].molecular_weight, gas_dat->species_dat[j].molecular_weight, gas_dat->species_dat[i].density, gas_dat->species_dat[j].density, gas_dat->species_dat[i].dynamic_viscosity, gas_dat->species_dat[j].dynamic_viscosity));
					
					//Enforce symmetry of the matrix
					gas_dat->binary_diffusion.edit(j, i, gas_dat->binary_diffusion(i,j));
					
				}
			}
			
			//Additive Properties
			gas_dat->total_molecular_weight = gas_dat->total_molecular_weight + (gas_dat->molefraction[i]*gas_dat->species_dat[i].molecular_weight);
			gas_dat->total_specific_heat = gas_dat->total_specific_heat + (gas_dat->molefraction[i]*gas_dat->species_dat[i].specific_heat);
			
		}
		
		//Calculate total density
		gas_dat->total_density = CE3(gas_dat->total_pressure*gas_dat->total_molecular_weight, gas_dat->gas_temperature);
		
		//Secondary loop to evaluate the total dynamic viscosity and molecular diffusion
		for (int i=0; i<gas_dat->N; i++)
		{
			double muT_sum = 0.0;
			double Dm_sum = 0.0;
			for (int j=0; j<gas_dat->N; j++)
			{
				if (j!=i)
				{
					//Evaluate summation terms
					Dm_sum = Dm_sum + (gas_dat->molefraction[j]/gas_dat->binary_diffusion(i,j));
					muT_sum = muT_sum + (gas_dat->molefraction[j]/Dp_ij(gas_dat->binary_diffusion(i,j),gas_dat->total_pressure));
				}
			}
			if (Dm_sum == 0.0)
				gas_dat->species_dat[i].molecular_diffusion = gas_dat->binary_diffusion(i,i);
			else
				gas_dat->species_dat[i].molecular_diffusion = (1.0 - gas_dat->molefraction[i]) / Dm_sum;
			if (gas_dat->molefraction[i] != 0.0)
				gas_dat->total_dyn_vis = gas_dat->total_dyn_vis + (gas_dat->species_dat[i].dynamic_viscosity/(1.0+((113.65*PSI(gas_dat->gas_temperature)*gas_dat->species_dat[i].dynamic_viscosity*gas_dat->gas_temperature)/(gas_dat->molefraction[i]*gas_dat->species_dat[i].molecular_weight))*muT_sum));
		}
		
		//Calculate remaining properties
		gas_dat->kinematic_viscosity = Nu(gas_dat->total_dyn_vis,gas_dat->total_density);
		gas_dat->Reynolds = ReNum(gas_dat->velocity, gas_dat->char_length, gas_dat->kinematic_viscosity);
		for (int i=0; i<gas_dat->N; i++)
			gas_dat->species_dat[i].Schmidt = ScNum(gas_dat->kinematic_viscosity, gas_dat->species_dat[i].molecular_diffusion);
	}
	
	return success;
}