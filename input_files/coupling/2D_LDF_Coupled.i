[GlobalParams]

	vx = 0.0 #R-direction (positive leaves outward, negative pushes inward)
	vy = 2.0 #Z-direction (positive moves bottom to top, negative moves top to bottom)
	vz = 0.0 #Not used in RZ system

	Dxx = 0.01 #Radial Diffusion in cylinder
	Dxy = 0.0
	Dxz = 0.0

	Dyx = 0.0
	Dyy = 0.01 #Axial Diffusion in cylinder
	Dyz = 0.0

	Dzx = 0.0
	Dzy = 0.0
	Dzz = 0.0

[]

[Problem]

	coord_type = RZ

[] #END Problem



[Mesh]

	type = GeneratedMesh
	dim = 2
	nx = 10
	ny = 40
	xmin = 0.0
	xmax = 0.25 #R
	ymin = 0.0
	ymax = 1.0 #Z

[] # END Mesh

[Variables]

	[./u]
		order = CONSTANT
		family = MONOMIAL
		initial_condition = 1.0
	[../]

	[./q]
		order = CONSTANT
		family = MONOMIAL
		initial_condition = 11.0
	[../]

[]

[AuxVariables]

[]

[ICs]

[]

[Kernels]

	[./u_dot]
		type = CoefTimeDerivative
		variable = u
		Coefficient = 1.0
	[../]

	[./u_gadv]
		type = GAdvection
		variable = u
	[../]

	[./u_gdiff]
		type = GAnisotropicDiffusion
		variable = u
	[../]

	[./q_dot]
		type = CoefTimeDerivative
		variable = q
		Coefficient = 1.0
	[../]

	[./q_ldf]
		type = CoupledLDF
		variable = q
		coupled = q
		driving_var = u
		ldf_coef = 4.0
		gaining = true
		driving_coef = 10.0
	[../]
 
	[./u_ldf]
		type = CoupledLDF
		variable = u
		coupled = q
		driving_var = u
		ldf_coef = 8.0
		gaining = false
		driving_coef = 10.0
	[../]

[]

[DGKernels]

	[./u_dgadv]
		type = DGAdvection
		variable = u
	[../]

	[./u_dgdiff]
		type = DGAnisotropicDiffusion
		variable = u
	[../]

[]


[AuxKernels]

[]

[BCs]

	[./u_bc_bot]
		type = DGFluxLimitedBC
		variable = u
		boundary = 'bottom'
		u_input = 1.0
	[../]

	[./u_bc_top]
		type = DGFluxLimitedBC
		variable = u
		boundary = 'top'
	[../]

[]

[Materials]

[]

[Postprocessors]

	[./u_out]
		type = SideAverageValue
		boundary = 'top'
		variable = u
	[../]
 
	[./q_avg]
		type = ElementAverageValue
		variable = q
	[../]

[]

[Executioner]

	type = Transient
	scheme = implicit-euler
	#scheme = bdf2

	nl_rel_tol = 1e-6
	nl_abs_tol = 1e-6
	nl_rel_step_tol = 1e-10
	nl_abs_step_tol = 1e-10
	l_tol = 1e-6
	l_max_its = 100

	solve_type = PJFNK
	line_search = bt
	start_time = 0.0
	end_time = 60.0
	petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
	petsc_options_value = 'hypre boomeramg 100'

	[./TimeStepper]
		#type = ConstantDT
		type = SolutionTimeAdaptiveDT
		dt = 0.01
	[../]

[]

[Outputs]
	exodus = true
	csv = true
	print_linear_residuals = true
[]
