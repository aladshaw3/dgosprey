[GlobalParams]

    vy = 2.0

    Dxx = 1.0
    Dxy = 1.0
    Dyx = 1.0
    Dyy = 1.0

    u_input = 1.0

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
    xmax = 1.0
    ymin = 0.0
    ymax = 1.0

[] # END Mesh

[Variables]

    [./u]
        order = CONSTANT
        family = MONOMIAL
        intital_condition = 0
    [../]


[] #END Variables

[AuxVariables]


[] #END AuxVariables

[ICs]


[] #END ICs

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

[] #END Kernels

[DGKernels]

    [./u_dgadv]
        type = DGAdvection
        variable = u
    [../]

    [./u_dgdiff]
        type = DGAnisotropicDiffusion
        variable = u
    [../]

[] #END DGKernels

[AuxKernels]


[] #END AuxKernels

[BCs]

    [./u_Flux]
        type = DGFluxBC
        variable = u
        boundary = 'top bottom left right'

    [../]


[] #END BCs

[Materials]


[] #END Materials

[Postprocessors]

    [./u_exit]
        type = SideAverageValue
        boundary = 'top'
        variable = u
        execute_on = 'initial timestep_end'
    [../]

    [./u_enter]
        type = SideAverageValue
        boundary = 'bottom'
        variable = u
        execute_on = 'initial timestep_end'
    [../]

    [./u_avg]
        type = ElementAverageValue
        variable = u
        execute_on = 'initial timestep_end'
    [../]

[] #END Postprocessors

[Executioner]

    type = Transient
    scheme = implicit-euler

    # NOTE: The default tolerances are far to strict and cause the program to crawl
    nl_rel_tol = 1e-6
    nl_abs_tol = 1e-6
    nl_rel_step_tol = 1e-10
    nl_abs_step_tol = 1e-10
    l_tol = 1e-6
    l_max_its = 100
    nl_max_its = 10

    solve_type = newton
    line_search = none    # Options: default shell none basic l2 bt cp
    start_time = 0.0
    end_time = 1.0
    dtmax = 0.1
    petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
    petsc_options_value = 'hypre boomeramg 100'

    [./TimeStepper]
        #Need to write a custom TimeStepper to enforce a maximum allowable dt
        type = ConstantDT
        #type = SolutionTimeAdaptiveDT
        dt = 0.05
    [../]

[] #END Executioner

[Preconditioning]

[] #END Preconditioning

[Outputs]

    exodus = true
    csv = true
    print_linear_residuals = true

[] #END Outputs
