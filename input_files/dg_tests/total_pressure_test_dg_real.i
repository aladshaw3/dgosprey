[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 10
  xmax = 12.7
[]

[Variables] #NOTE: ALL VARIABLES NEED CUSTOM ICs and BCs to handle coupling appropriately

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

[]

[AuxVariables]
  [./ambient_temp]
    order = FIRST
    family = L2_LAGRANGE
    initial_condition = 298.15
  [../]
 	[./total_pressure]
 		order = FIRST
 		family = L2_LAGRANGE
 		initial_condition = 101.35
 	[../]
[]

[ICs]
 [./N2_IC]
		type = ConcentrationIC
    variable = N2
  	initial_mole_frac = 0.79
		initial_temp = 298.15
 		initial_press = 101.35
	[../]
	[./O2_IC]
		type = ConcentrationIC
  	variable = O2
  	initial_mole_frac = 0.21
 		initial_temp = 298.15
 		initial_press = 101.35
	[../]
  [./H2O_IC]
		type = ConcentrationIC
  	variable = H2O
  	initial_mole_frac = 0.0
 		initial_temp = 298.15
 		initial_press = 101.35
	[../]
[]

[Kernels]

  [./accumN2]
    type = RetardedTimeDerivative
    variable = N2
    index = 0
  [../]
  [./diffN2]
    type = Dispersion
    variable = N2
    index = 0
  [../]
  [./advN2]
    type = ConservativeAdvection
    variable = N2
    dir = 0
  [../]

  [./accumO2]
    type = RetardedTimeDerivative
    variable = O2
    index = 1
  [../]
  [./diffO2]
    type = Dispersion
    variable = O2
    index = 1
  [../]
  [./advO2]
    type = ConservativeAdvection
    variable = O2
    dir = 0
  [../]

  [./accumH2O]
    type = RetardedTimeDerivative
    variable = H2O
    index = 2
  [../]
  [./diffH2O]
    type = Dispersion
    variable = H2O
    index = 2
  [../]
  [./advH2O]
    type = ConservativeAdvection
    variable = H2O
    dir = 0
  [../]
  [./H2O_penalty]
  	type = MassPenalty
    variable = H2O
    penalty = 1.0
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

  [./bedAccum]
    type = BedHeatAccumulation
    variable = column_temp
  [../]
  [./wall_fluid_trans]
    type = FluidHeatTransfer
    variable = column_temp
    coupled = wall_temp
  [../]
  [./heat_conductivity]	#change
    type = HeatConductivity
    variable = column_temp
  [../]
  [./heat_convection]	#change
    type = HeatConvection
    variable = column_temp
    dir = 0
  [../]

[]

[DGKernels]
	[./dg_disp_N2]
		type = DGDispersion
		variable = N2
		epsilon = 1.0
		sigma = 1.0
 		index = 0
	[../]
 	[./dg_disp_O2]
 		type = DGDispersion
 		variable = O2
 		epsilon = 1.0
 		sigma = 1.0
 		index = 1
 	[../]
 	[./dg_disp_H2O]
 		type = DGDispersion
 		variable = H2O
 		epsilon = 1.0
 		sigma = 1.0
 		index = 2
 	[../]

  [./dg_adv_N2]
  	type = DGConservativeAdvection
    variable = N2
    dir = 0
  [../]
  [./dg_adv_O2]
  	type = DGConservativeAdvection
    variable = O2
    dir = 0
  [../]
  [./dg_adv_H2O]
  	type = DGConservativeAdvection
    variable = H2O
    dir = 0
  [../]
[]

[AuxKernels]
	[./column_pressure]
      type = TotalColumnPressure
      variable = total_pressure
      temperature = column_temp
      coupled_gases = 'N2 O2 H2O'
	[../]
[]

[BCs]

  [./left_N2]
    type = DGCoupledMassBC
    variable = N2
    boundary = 0
    inlet_temperature = 298.15
    inlet_pressure = 101.35
    input_mole_frac = 0.78863
    index = 0
    dir = 0
  [../]
  [./right_N2]
    type = DGMassOutFlowBC
    variable = N2
    boundary = 1
    value = 0
    dir = 0
  [../]
  [./left_O2]
    type = DGCoupledMassBC
    variable = O2
    boundary = 0
    inlet_temperature = 298.15
    inlet_pressure = 101.35
    input_mole_frac = 0.20974
    index = 1
    dir = 0
  [../]
  [./right_O2]
    type = DGMassOutFlowBC
    variable = O2
    boundary = 1
    value = 0
    dir = 0
  [../]
  [./left_H2O]
    type = DGCoupledMassBC
    variable = H2O
    boundary = 0
    inlet_temperature = 298.15
    inlet_pressure = 101.35
    input_mole_frac = 0.00163
    index = 2
    dir = 0
  [../]
  [./right_H2O]
    type = DGMassOutFlowBC
    variable = H2O
    boundary = 1
    value = 0
    dir = 0
  [../]

  [./left_temp]
    type = HeatCauchyBC
    variable = column_temp
    boundary = 0
    boundary_value = 298.15
  [../]
  [./right_temp]
    type = NeumannBC
    variable = column_temp
    boundary = 1
    value = 0
  [../]
[]

[Materials]
  [./BedMaterials]
    type = BedProperties
    block = 0
    length = 12.7 # same as xmax (always)
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
    molecular_wieght = '28.016 32 18' # note: all vectorized quantities must be in same order as variable coupling!!!
    comp_heat_capacity = '1.04 0.919 1.97'
    comp_ref_viscosity = '0.0001781 0.0002018 0.0001043'
    comp_ref_temp = '300.55 292.25 298.16'
    comp_Sutherland_const = '111 127 784.72'
    flow_rate = 211680.0
    temperature = column_temp
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
[]

[Postprocessors]
  [./N2_exit]
    #type = ElementalVariableValue
    #elementid = 0
    type = SideAverageValue
    boundary = 1
    variable = N2
  [../]
  [./O2_exit]
    #type = ElementalVariableValue
    #elementid = 0
    type = SideAverageValue
    boundary = 1
    variable = O2
  [../]
  [./H2O_exit]
    #type = ElementalVariableValue
    #elementid = 0
    type = SideAverageValue
    boundary = 1
    variable = H2O
  [../]
  [./temp_exit]
    #type = ElementalVariableValue
    #elementid = 0
    type = SideAverageValue
    boundary = 1
    variable = column_temp
  [../]
  [./press_exit]
    #type = ElementalVariableValue
    #elementid = 0
    type = SideAverageValue
    boundary = 1
    variable = total_pressure
  [../]
  [./wall_temp_exit]
    #type = ElementalVariableValue
    #elementid = 0
    type = SideAverageValue
    boundary = 1
    variable = wall_temp
  [../]
  [./max_wave_front]
 		type = WaveFront
 		execute_on = timestep #This line is necessary to prevent wave_front from becoming random and causing errors
 	[../]
	[./avg_elem_size]
		type = AverageElementSize
 		variable = H2O
	[../]
[]

[Executioner]
 	#NOTE: For better results, we must forcably relax ALL residual tolerances
  type = Transient
  scheme = bdf2

  nl_rel_tol = 1e-06
  picard_abs_tol = 1e-10
  nl_abs_tol = 1e-10
  nl_rel_step_tol = 1e-10
  picard_rel_tol = 1e-10
  nl_abs_step_tol = 1e-10

  solve_type = PJFNK
  start_time = 0.0
  end_time = 60.0
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  [./TimeStepper]
		#type = CFLConditionApproximation
		type = SolutionTimeAdaptiveDT  #Note: This works well now that we are mass conservative!!!
		dt = 0.01
		mesh_size = avg_elem_size
 		wave_front = max_wave_front
  [../]
[]

[Adaptivity]
  # NOTE: you may need one of these for each variable that has an advective term
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
[]

[Outputs]
  exodus = true
  csv = true
[]
