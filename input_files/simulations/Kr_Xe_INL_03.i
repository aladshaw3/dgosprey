[GlobalParams]

 [] #END GlobalParams

 [Problem]

 	coord_type = RZ

 [] #END Problem

[Mesh]

	type = GeneratedMesh
	dim = 2
	nx = 10
 	ny = 40
 	xmin = 0.0
	xmax = 0.8636 #cm
 	ymin = 0.0
	ymax = 22.86 #cm

 [] # END Mesh

[Variables]

	[./Kr]
		order = CONSTANT
		family = MONOMIAL
	[../]

	[./Xe]
		order = CONSTANT
		family = MONOMIAL
	[../]

	[./He]
		order = CONSTANT
		family = MONOMIAL
	[../]

 	[./wall_temp]
 		order = CONSTANT
 		family = MONOMIAL
 		initial_condition = 295.15
 	[../]

	[./column_temp]
 		order = CONSTANT
 		family = MONOMIAL
 		initial_condition = 295.15
	[../]

 [] #END Variables

[AuxVariables]

	[./total_pressure]
		order = CONSTANT
		family = MONOMIAL
		initial_condition = 101.35
	[../]

 	[./ambient_temp]
 		order = CONSTANT
 		family = MONOMIAL
 		initial_condition = 295.15
 	[../]

	[./He_Adsorbed]
		order = CONSTANT
		family = MONOMIAL
		initial_condition = 0.0
	[../]

	[./Kr_Adsorbed]
		order = CONSTANT
		family = MONOMIAL
		initial_condition = 0.0
	[../]

	[./Xe_Adsorbed]
		order = CONSTANT
		family = MONOMIAL
		initial_condition = 0.0
	[../]

	[./He_Perturb]
		order = CONSTANT
		family = MONOMIAL
		initial_condition = 0.0
	[../]

	[./Kr_Perturb]
		order = CONSTANT
		family = MONOMIAL
		initial_condition = 0.0
	[../]

	[./Xe_Perturb]
		order = CONSTANT
		family = MONOMIAL
		initial_condition = 0.0
	[../]

	[./Kr_AdsorbedHeat]
		order = CONSTANT
		family = MONOMIAL
		initial_condition = 0.0
	[../]

	[./Xe_AdsorbedHeat]
		order = CONSTANT
		family = MONOMIAL
		initial_condition = 0.0
	[../]

	[./He_AdsorbedHeat]
		order = CONSTANT
		family = MONOMIAL
		initial_condition = 0.0
	[../]

 [] #END AuxVariables

[ICs]

	[./Kr_IC]
		type = ConcentrationIC
		variable = Kr
		initial_mole_frac = 0.0
		initial_press = 101.35
		initial_temp = 295.15
	[../]

	[./Xe_IC]
		type = ConcentrationIC
		variable = Xe
		initial_mole_frac = 0.0
 		initial_press = 101.35
 		initial_temp = 295.15
	[../]

	[./He_IC]
		type = ConcentrationIC
		variable = He
		initial_mole_frac = 1.0
 		initial_press = 101.35
 		initial_temp = 295.15
	[../]

 [] #END ICs

[Kernels]

 	[./accumKr]
		type = BedMassAccumulation
 		variable = Kr
		index = 0
 	[../]

	[./diffKr]
		type = GColumnMassDispersion
		variable = Kr
		index = 0
	[../]

	[./advKr]
		type = GColumnMassAdvection
		variable = Kr
	[../]

 	[./accumXe]
		type = BedMassAccumulation
 		variable = Xe
		index = 1
 	[../]

	[./diffXe]
		type = GColumnMassDispersion
		variable = Xe
		index = 1
	[../]

	[./advXe]
		type = GColumnMassAdvection
		variable = Xe
	[../]

 	[./accumHe]
		type = BedMassAccumulation
 		variable = He
		index = 2
 	[../]

	[./diffHe]
		type = GColumnMassDispersion
		variable = He
		index = 2
	[../]

	[./advHe]
		type = GColumnMassAdvection
		variable = He
	[../]

 	[./wallAccum]
 		type = WallHeatAccumulation
 		variable = wall_temp
 	[../]
 	[./wall_bed_trans]
 		type = BedWallHeatTransfer
 		variable = wall_temp
 		coupled = column_temp
 	[../]
 	[./wall_amb_trans]
 		type = WallAmbientHeatTransfer
 		variable = wall_temp
 		coupled = ambient_temp
 	[../]

	[./columnAccum]
		type = BedHeatAccumulation
		variable = column_temp
	[../]
	[./columnConduction]
		type = GColumnHeatDispersion
		variable =column_temp
	[../]
	[./columnAdvection]
		type = GColumnHeatAdvection
		variable =column_temp
	[../]
	[./columnAdsHeat]
		type = AdsorptionHeatAccumulation
		variable = column_temp
		solid_heats = 'Kr_AdsorbedHeat Xe_AdsorbedHeat He_AdsorbedHeat'
	[../]


 [] #END Kernels

[DGKernels]

	[./dg_disp_Kr]
		type = DGColumnMassDispersion
		variable = Kr
		index = 0
	[../]

 	[./dg_adv_Kr]
		type = DGColumnMassAdvection
		variable = Kr
	[../]

	[./dg_disp_Xe]
		type = DGColumnMassDispersion
		variable = Xe
		index = 1
	[../]

 	[./dg_adv_Xe]
		type = DGColumnMassAdvection
		variable = Xe
	[../]

	[./dg_disp_He]
		type = DGColumnMassDispersion
		variable = He
		index = 2
	[../]

	[./dg_adv_He]
		type = DGColumnMassAdvection
		variable = He
	[../]

	[./dg_disp_heat]
		type = DGColumnHeatDispersion
		variable = column_temp
	[../]

	[./dg_adv_heat]
		type = DGColumnHeatAdvection
		variable = column_temp
	[../]

 [] #END DGKernels

[AuxKernels]

	[./column_pressure]
		type = TotalColumnPressure
		variable = total_pressure
		temperature = column_temp
		coupled_gases = 'Kr Xe He'
	[../]

	[./krypton_adsorption]
		type = MAGPIE_Adsorption
		variable = Kr_Adsorbed
		index = 0
		execute_on = 'initial timestep_end'
	[../]

	[./xenon_adsorption]
		type = MAGPIE_Adsorption
		variable = Xe_Adsorbed
		index = 1
		execute_on = 'initial timestep_end'
	[../]

	[./helium_adsorption]
		type = MAGPIE_Adsorption
		variable = He_Adsorbed
		index = 2
		execute_on = 'initial timestep_end'
	[../]

	[./krypton_perturbation]
		type = MAGPIE_Perturbation
		variable = Kr_Perturb
		index = 0
		execute_on = 'initial timestep_end'
	[../]

	[./xenon_perturbation]
		type = MAGPIE_Perturbation
		variable = Xe_Perturb
		index = 1
		execute_on = 'initial timestep_end'
	[../]

	[./helium_perturbation]
		type = MAGPIE_Perturbation
		variable = He_Perturb
		index = 2
		execute_on = 'initial timestep_end'
	[../]

	[./krypton_adsorption_heat]
		type = MAGPIE_AdsorptionHeat
		variable = Kr_AdsorbedHeat
		solid_conc = Kr_Adsorbed
		index = 0
	[../]

	[./xenon_adsorption_heat]
		type = MAGPIE_AdsorptionHeat
		variable = Xe_AdsorbedHeat
		solid_conc = Xe_Adsorbed
		index = 1
	[../]

	[./helium_adsorption_heat]
		type = MAGPIE_AdsorptionHeat
		variable = He_AdsorbedHeat
		solid_conc = He_Adsorbed
		index = 2
	[../]

 [] #END AuxKernels

[BCs]

 	[./Kr_Flux]
		type = DGMassFluxLimitedBC
 		variable = Kr
 		boundary = 'top bottom'
 		input_temperature = 295.15
 		input_pressure = 101.35
 		input_molefraction = 0.000153814
 		index = 0
 	[../]

 	[./Xe_Flux]
		type = DGMassFluxLimitedBC
 		variable = Xe
 		boundary = 'top bottom'
 		input_temperature = 295.15
 		input_pressure = 101.35
 		input_molefraction = 0.001007334
 		index = 1
 	[../]

 	[./He_Flux]
		type = DGMassFluxLimitedBC
 		variable = He
 		boundary = 'top bottom'
 		input_temperature = 295.15
 		input_pressure = 101.35
 		input_molefraction = 0.998838852
 		index = 2
 	[../]

	[./Heat_Gas_Flux]
 		type = DGHeatFluxLimitedBC
 		variable = column_temp
 		boundary = 'top bottom'
 		input_temperature = 295.15
 	[../]
 
	[./Heat_Wall_Flux]
		type = DGColumnWallHeatFluxLimitedBC
		variable = column_temp
		boundary = 'right left'
		wall_temp = wall_temp
	[../]

 [] #END BCs

[Materials]

	[./BedMaterials]
		type = BedProperties
		block = 0
		length = 22.86
		inner_diameter = 1.7272
		outer_diameter = 1.905
		bulk_porosity = 0.798				#not known
		axial_conductivity = 6.292E-05      #not known
		wall_density = 7.7
		wall_heat_capacity = 0.5
		wall_heat_trans_coef = 9.0
		extern_heat_trans_coef = 9.0
		temperature = column_temp
		coupled_gases = 'Kr Xe He'
	[../]

	[./FlowMaterials]
		type = FlowProperties
		block = 0
		molecular_wieght = '83.8 131.29 4.0026'
		comp_heat_capacity = '0.25 0.16 5.1916'
		comp_ref_viscosity = '0.00023219 0.00021216 0.0001885'
		comp_ref_temp = '273.15 273.15 273.15'
		comp_Sutherland_const = '266.505 232.746 80.0'
		flow_rate = 2994.06
		column_length = 22.86
		temperature = column_temp
 		total_pressure = total_pressure
		coupled_gases = 'Kr Xe He'
		coupled_adsorption = 'Kr_Adsorbed Xe_Adsorbed He_Adsorbed'
		coupled_perturbation = 'Kr_Perturb Xe_Perturb He_Perturb'
	[../]

	[./AdsorbentMaterials]
		type = AdsorbentProperties
		block = 0
		binder_fraction = 0.175				#not known
		binder_porosity = 0.27				#not known
		crystal_radius = 1.5				#not known
		pellet_diameter = 0.236				#not known
		macropore_radius = 3.5e-6			#not Known
		pellet_density = 1.69				#not Known
		pellet_heat_capacity = 1.045		#not known
		temperature = column_temp
		coupled_gases = 'Kr Xe He'
	[../]

	[./AdsorbateMaterials]
		type = MagpieAdsorbateProperties
		block = 0
		temperature = column_temp
		total_pressure = total_pressure
		coupled_gases = 'Kr Xe He'
		number_sites = '2 3 0'
		maximum_capacity = '1.716 1.479 0' #mol/kg
		molar_volume = '20.785 25.412 0' #cm^3/mol
		enthalpy_site_1 = '-44696.86 -18455.18 0'
		enthalpy_site_2 = '-65465.52 -35511.74 0'
		enthalpy_site_3 = '0 -53315.13 0'
		enthalpy_site_4 = '0 0 0'
		enthalpy_site_5 = '0 0 0'
		enthalpy_site_6 = '0 0 0'

		entropy_site_1 = '-170.45 -23.25 0'
		entropy_site_2 = '-248.55 -62.45 0'
		entropy_site_3 = '0 -100.10 0'
		entropy_site_4 = '0 0 0'
		entropy_site_5 = '0 0 0'
		entropy_site_6 = '0 0 0'
	[../]

 [] #END Materials

[Postprocessors]

	[./Kr_exit]
		type = SideAverageValue
		boundary = 'top'
		variable = Kr
		execute_on = timestep_end
	[../]

	[./Xe_exit]
		type = SideAverageValue
		boundary = 'top'
		variable = Xe
		execute_on = timestep_end
	[../]

	[./He_exit]
		type = SideAverageValue
		boundary = 'top'
		variable = He
		execute_on = timestep_end
	[../]

	[./temp_exit]
		type = SideAverageValue
		boundary = 'top'
		variable = column_temp
		execute_on = timestep_end
	[../]

	[./press_exit]
		type = SideAverageValue
		boundary = 'top'
		variable = total_pressure
		execute_on = timestep_end
	[../]

 	[./wall_temp]
 		type = SideAverageValue
 		boundary = 'right'
 		variable = wall_temp
		execute_on = timestep_end
 	[../]

	[./Kr_solid]
		type = ElementAverageValue
		variable = Kr_Adsorbed
		execute_on = timestep_end
	[../]

	[./Kr_heat]
		type = ElementAverageValue
		variable = Kr_AdsorbedHeat
		execute_on = timestep_end
	[../]
 
	[./Xe_solid]
		type = ElementAverageValue
		variable = Xe_Adsorbed
		execute_on = timestep_end
	[../]
 
	[./Xe_heat]
		type = ElementAverageValue
		variable = Xe_AdsorbedHeat
		execute_on = timestep_end
	[../]

 [] #END Postprocessors

[Executioner]

 	type = Transient
	scheme = implicit-euler

	# NOTE: The default tolerances are far to strict and cause the program to crawl
 	nl_rel_tol = 1e-6
 	nl_abs_tol = 1e-4
 	nl_rel_step_tol = 1e-10
 	nl_abs_step_tol = 1e-10
 	l_tol = 1e-6
 	l_max_its = 100

	solve_type = pjfnk
    line_search = bt    # Options: default shell none basic l2 bt cp
	start_time = 0.0
	end_time = 10.0
    petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
    petsc_options_value = 'hypre boomeramg 100'

	[./TimeStepper]
		#Need to write a custom TimeStepper to enforce a maximum allowable dt
		#type = ConstantDT
		type = SolutionTimeAdaptiveDT
		dt = 0.01
	[../]

 [] #END Executioner

[Preconditioning]

[] #END Preconditioning

[Outputs]

	exodus = true
	csv = true
	print_linear_residuals = true

 [] #END Outputs
