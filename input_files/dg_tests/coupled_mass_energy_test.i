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
	xmax = 2.54
 	ymin = 0.0
 	ymax = 12.7

 [] # END Mesh

[Variables]

	[./N2]
		order = CONSTANT
		family = MONOMIAL
	[../]

	[./O2]
		order = CONSTANT
		family = MONOMIAL
	[../]

	[./H2O]
		order = CONSTANT
		family = MONOMIAL
	[../]

 	[./wall_temp]
 		order = CONSTANT
 		family = MONOMIAL
 		initial_condition = 298.15
 	[../]

	[./column_temp]
 		order = CONSTANT
 		family = MONOMIAL
 		initial_condition = 298.15
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
 		initial_condition = 298.15
 	[../]
 
	[./H2O_Adsorbed]
		order = CONSTANT
		family = MONOMIAL
		initial_condition = 0.0
	[../]

 [] #END AuxVariables

[ICs]

	[./N2_IC]
		type = DGConcentrationIC
		variable = N2
		initial_mole_frac = 0.79
		initial_press = 101.35
		initial_temp = 298.15
	[../]

	[./O2_IC]
		type = DGConcentrationIC
		variable = O2
		initial_mole_frac = 0.21
 		initial_press = 101.35
 		initial_temp = 298.15
	[../]

	[./H2O_IC]
		type = DGConcentrationIC
		variable = H2O
		initial_mole_frac = 0.0
 		initial_press = 101.35
 		initial_temp = 298.15
	[../]

 [] #END ICs

[Kernels]

 	[./accumN2]
 		type = RetardedTimeDerivative
 		variable = N2
 		index = 0
 	[../]

	[./diffN2]
		type = GColumnMassDispersion
		variable = N2
		index = 0
	[../]

	[./advN2]
		type = GColumnMassAdvection
		variable = N2
	[../]

 	[./accumO2]
 		type = RetardedTimeDerivative
 		variable = O2
 		index = 1
 	[../]

	[./diffO2]
		type = GColumnMassDispersion
		variable = O2
		index = 1
	[../]

	[./advO2]
		type = GColumnMassAdvection
		variable = O2
	[../]

 	[./accumH2O]
 		type = RetardedTimeDerivative
 		variable = H2O
 		index = 2
 	[../]

	[./diffH2O]
		type = GColumnMassDispersion
		variable = H2O
		index = 2
	[../]

	[./advH2O]
		type = GColumnMassAdvection
		variable = H2O
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


 [] #END Kernels

[DGKernels]

	[./dg_disp_N2]
		type = DGColumnMassDispersion
		variable = N2
		index = 0
	[../]

 	[./dg_adv_N2]
		type = DGColumnMassAdvection
		variable = N2
	[../]

	[./dg_disp_O2]
		type = DGColumnMassDispersion
		variable = O2
		index = 1
	[../]

 	[./dg_adv_O2]
		type = DGColumnMassAdvection
		variable = O2
	[../]

	[./dg_disp_H2O]
		type = DGColumnMassDispersion
		variable = H2O
		index = 2
	[../]

	[./dg_adv_H2O]
		type = DGColumnMassAdvection
		variable = H2O
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
		coupled_gases = 'N2 O2 H2O'
	[../]
 
	[./water_adsorption]
		type = MAGPIE_Adsorption
		variable = H2O_Adsorbed
		index = 2
	[../]

 [] #END AuxKernels

[BCs]

 	[./N2_Flux]
 		type = DGMassFluxLimitedBC
 		variable = N2
 		boundary = 'top bottom'
 		input_temperature = 298.15
 		input_pressure = 101.35
 		input_molefraction = 0.78863
 		index = 0
 	[../]

 	[./O2_Flux]
 		type = DGMassFluxLimitedBC
 		variable = O2
 		boundary = 'top bottom'
 		input_temperature = 298.15
 		input_pressure = 101.35
 		input_molefraction = 0.20974
 		index = 1
 	[../]

 	[./H2O_Flux]
 		type = DGMassFluxLimitedBC
 		variable = H2O
 		boundary = 'top bottom'
 		input_temperature = 298.15
 		input_pressure = 101.35
 		input_molefraction = 0.00163
 		index = 2
 	[../]

	[./Heat_Gas_Flux]
 		type = DGHeatFluxLimitedBC
 		variable = column_temp
 		boundary = 'top bottom'
 		input_temperature = 298.15
 	[../]

	[./Heat_Wall_Flux]
 		type = DGColumnWallHeatFluxLimitedBC
 		variable = column_temp
 		boundary = 'right'
 		wall_temp = wall_temp
 	[../]

 [] #END BCs

[Materials]

	[./BedMaterials]
		type = BedProperties
		block = 0
		length = 12.7
		din = 2.54
		dout = 3.54
		eb = 0.421
		Kz = 6.292E-05
		rhow = 1.0
		hw = 1.0
		Uw = 1.0
		Ua = 1.0
		temperature = column_temp
		coupled_gases = 'N2 O2 H2O'
	[../]

	[./FlowMaterials]
		type = FlowProperties
		block = 0
		molecular_wieght = '28.016 32 18'
		comp_heat_capacity = '1.04 0.919 1.97'
		comp_ref_viscosity = '0.0001781 0.0002018 0.0001043'
		comp_ref_temp = '300.55 292.25 298.16'
		comp_Sutherland_const = '111 127 784.72'
		flow_rate = 211680.0
		column_length = 12.7
		temperature = column_temp
 		total_pressure = total_pressure
		coupled_gases = 'N2 O2 H2O'
	[../]

	[./AdsorbentMaterials]
		type = AdsorbentProperties
		block = 0
		binder_fraction = 0.175
		binder_porosity = 0.27
		crystal_radius = 1.5
		pellet_diameter = 0.236
		macropore_radius = 3.5e-6
		rhos = 1.41
		hs = 1.045
		temperature = column_temp
		coupled_gases = 'N2 O2 H2O'
	[../]
	
	[./AdsorbateMaterials]
		type = AdsorbateProperties
		block = 0
		adsorption_type = 0
		temperature = column_temp
		total_pressure = total_pressure
		coupled_gases = 'N2 O2 H2O'
		number_sites = '0 0 4'
		maximum_capacity = '0 0 11.67'
		molar_volume = '0 0 13.91'
		enthalpy_site_1 = '0 0 -46597.5'
		enthalpy_site_2 = '0 0 -125024'
		enthalpy_site_3 = '0 0 -193619'
		enthalpy_site_4 = '0 0 -272228'
		enthalpy_site_5 = '0 0 0'
		enthalpy_site_6 = '0 0 0'
		
		entropy_site_1 = '0 0 -53.6994'
		entropy_site_2 = '0 0 -221.073'
		entropy_site_3 = '0 0 -356.728'
		entropy_site_4 = '0 0 -567.459'
		entropy_site_5 = '0 0 0'
		entropy_site_6 = '0 0 0'
	[../]

 [] #END Materials

[Postprocessors]

	[./N2_exit]
		type = SideAverageValue
		boundary = 'top'
		variable = N2
		execute_on = timestep_end
	[../]

	[./O2_exit]
		type = SideAverageValue
		boundary = 'top'
		variable = O2
		execute_on = timestep_end
	[../]

	[./H2O_exit]
		type = SideAverageValue
		boundary = 'top'
		variable = H2O
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
 
	[./H2O_avg_sorption]
		type = ElementAverageValue
		variable = H2O_Adsorbed
		execute_on = timestep_end
	[../]

 [] #END Postprocessors

[Executioner]

 	type = Transient
	scheme = bdf2

 	nl_rel_tol = 1e-6
 	picard_abs_tol = 1e-6
 	nl_abs_tol = 1e-6
 	nl_rel_step_tol = 1e-6
 	picard_rel_tol = 1e-6
 	nl_abs_step_tol = 1e-6
 	l_tol = 1e-4
 	l_max_its = 100

 	solve_type = PJFNK
	start_time = 0.0
	end_time = 60.0
	petsc_options_iname = '-pc_type -pc_hypre_type'
	petsc_options_value = 'hypre boomeramg'

	[./TimeStepper]
		type = SolutionTimeAdaptiveDT
		dt = 0.1
	[../]

 [] #END Executioner

[Outputs]

 output_initial = false
 exodus = true
 csv = true
 print_perf_log = true

 [] #END Outputs
