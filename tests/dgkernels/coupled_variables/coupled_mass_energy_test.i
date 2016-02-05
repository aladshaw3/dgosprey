[GlobalParams]

 [] #END GlobalParams

 [Problem]

 	coord_type = RZ

 [] #END Problem

[Mesh]

	type = GeneratedMesh
	dim = 2
	nx = 3
 	ny = 5
 	xmin = 0.0
	xmax = 2.54
 	ymin = 0.0
 	ymax = 12.7

 [] # END Mesh

[Variables]

	[./N2]
		order = FIRST
		family = L2_LAGRANGE
	[../]

	[./O2]
		order = FIRST
		family = L2_LAGRANGE
	[../]

	[./H2O]
		order = FIRST
		family = L2_LAGRANGE
	[../]

 	[./wall_temp]
 		order = FIRST
 		family = L2_LAGRANGE
 		initial_condition = 298.15
 	[../]

	[./column_temp]
 		order = FIRST
 		family = L2_LAGRANGE
 		initial_condition = 298.15
	[../]

 [] #END Variables

[AuxVariables]

	[./total_pressure]
		order = FIRST
		family = L2_LAGRANGE
		initial_condition = 101.35
	[../]

 	[./ambient_temp]
 		order = FIRST
 		family = L2_LAGRANGE
 		initial_condition = 298.15
 	[../]

	[./H2O_Adsorbed]
		order = CONSTANT
		family = MONOMIAL
		initial_condition = 0.0
	[../]

	[./N2_Adsorbed]
		order = CONSTANT
		family = MONOMIAL
		initial_condition = 0.0
	[../]

	[./O2_Adsorbed]
		order = CONSTANT
		family = MONOMIAL
		initial_condition = 0.0
	[../]

	[./H2O_Perturb]
		order = CONSTANT
		family = MONOMIAL
		initial_condition = 0.0
	[../]

	[./N2_Perturb]
		order = CONSTANT
		family = MONOMIAL
		initial_condition = 0.0
	[../]

	[./O2_Perturb]
		order = CONSTANT
		family = MONOMIAL
		initial_condition = 0.0
	[../]

 [] #END AuxVariables

[ICs]

	[./N2_IC]
		type = ConcentrationIC
		variable = N2
		initial_mole_frac = 0.79
		initial_press = 101.35
		initial_temp = 298.15
	[../]

	[./O2_IC]
		type = ConcentrationIC
		variable = O2
		initial_mole_frac = 0.21
 		initial_press = 101.35
 		initial_temp = 298.15
	[../]

	[./H2O_IC]
		type = ConcentrationIC
		variable = H2O
		initial_mole_frac = 0.0
 		initial_press = 101.35
 		initial_temp = 298.15
	[../]

 [] #END ICs

[Kernels]

	[./diffN2]
		type = GColumnMassDispersion
		variable = N2
		index = 0
	[../]

	[./advN2]
		type = GColumnMassAdvection
		variable = N2
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

	[./diffH2O]
		type = GColumnMassDispersion
		variable = H2O
		index = 2
	[../]

	[./advH2O]
		type = GColumnMassAdvection
		variable = H2O
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
		inner_diameter = 2.54
		outer_diameter = 3.54
		bulk_porosity = 0.421
		axial_conductivity = 6.292E-05
		wall_density = 1.0
		wall_heat_capacity = 1.0
		wall_heat_trans_coef = 1.0
		extern_heat_trans_coef = 1.0
		temperature = column_temp
		coupled_gases = 'N2 O2 H2O'
	[../]

	[./FlowMaterials]
		type = FlowProperties
		block = 0
		molecular_weight = '28.016 32 18'
		comp_heat_capacity = '1.04 0.919 1.97'
		comp_ref_viscosity = '0.0001781 0.0002018 0.0001043'
		comp_ref_temp = '300.55 292.25 298.16'
		comp_Sutherland_const = '111 127 784.72'
		flow_rate = 211680.0
		length = 12.7
		temperature = column_temp
 		total_pressure = total_pressure
		coupled_gases = 'N2 O2 H2O'
		coupled_adsorption = 'N2_Adsorbed O2_Adsorbed H2O_Adsorbed'
		coupled_perturbation = 'N2_Perturb O2_Perturb H2O_Perturb'
	[../]

	[./AdsorbentMaterials]
		type = AdsorbentProperties
		block = 0
		binder_fraction = 0.175
		binder_porosity = 0.27
		crystal_radius = 1.5
		pellet_diameter = 0.236
		macropore_radius = 3.5e-6
		pellet_density = 1.41
		pellet_heat_capacity = 1.045
		temperature = column_temp
		coupled_gases = 'N2 O2 H2O'
	[../]

 [] #END Materials

[Postprocessors]

	[./N2_exit]
		type = SideAverageValue
		boundary = 'top'
		variable = N2
	[../]

	[./O2_exit]
		type = SideAverageValue
		boundary = 'top'
		variable = O2
	[../]

	[./H2O_exit]
		type = SideAverageValue
		boundary = 'top'
		variable = H2O
	[../]

	[./temp_exit]
		type = SideAverageValue
		boundary = 'top'
		variable = column_temp
	[../]

	[./press_exit]
		type = SideAverageValue
		boundary = 'top'
		variable = total_pressure
	[../]

 	[./wall_temp]
 		type = SideAverageValue
 		boundary = 'right'
 		variable = wall_temp
 	[../]

 [] #END Postprocessors

[Executioner]

 	type = Steady
 	scheme = bdf2

 	nl_rel_tol = 1e-6
 	picard_abs_tol = 1e-6
 	nl_abs_tol = 1e-6
 	nl_rel_step_tol = 1e-6
 	picard_rel_tol = 1e-6
 	nl_abs_step_tol = 1e-6
 	l_tol = 0.01
 	l_max_its = 100

 	solve_type = PJFNK
	start_time = 0.0
	end_time = 60.0
	petsc_options_iname = '-pc_type -pc_hypre_type'
	petsc_options_value = 'hypre boomeramg'

	[./TimeStepper]
		type = ConstantDT
		dt = 0.1
	[../]

 [] #END Executioner

[Adaptivity]

 	[./Indicators]
 		[./error_mass_N2]
 			type = GradientJumpIndicator
 			variable = N2
 		[../]
 		[./error_mass_O2]
 			type = GradientJumpIndicator
 			variable = O2
 		[../]
 		[./error_mass_H2O]
 			type = GradientJumpIndicator
 			variable = H2O
 		[../]
 		[./error_energy]
 			type = GradientJumpIndicator
 			variable = column_temp
 		[../]
 	[../]


 	[./Markers]
		[./ef_mass_N2]
 			type = ErrorFractionMarker
 			indicator = error_mass_N2
 		[../]
 		[./ef_mass_O2]
 			type = ErrorFractionMarker
 			indicator = error_mass_O2
 		[../]
 		[./ef_mass_H2O]
 			type = ErrorFractionMarker
 			indicator = error_mass_H2O
 		[../]
 		[./ef_energy]
 			type = ErrorFractionMarker
 			indicator = error_energy
 		[../]
 	[../]

 [] #END Adaptivity

[Outputs]
  execute_on = 'timestep_end'

 exodus = true
 csv = true
 [] #END Outputs
