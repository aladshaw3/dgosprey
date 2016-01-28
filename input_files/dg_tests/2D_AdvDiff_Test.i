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
[]

[Mesh]
type = GeneratedMesh
dim = 2
nx = 10 #R-direction == x
ny = 40 #Z-direction == y
xmax = 0.25
ymax = 1.0
[]

[Variables]

[./u]
order = FIRST
family = L2_LAGRANGE
initial_condition = 0.0
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

[./u_bc_top_bot]
type = DGFluxBC
variable = u
boundary = 'top bottom'
u_input = 1.0
[../]

[./u_bc_left_right]
type = DGFluxBC
variable = u
boundary = 'right' #Note: RZ applies natural BC at left Boundary
u_input = 1.0
[../]

[]

[Materials]

[]

[Postprocessors]

[./u_top]
type = SideAverageValue
boundary = 'top'
variable = u
[../]

[./u_bottom]
type = SideAverageValue
boundary = 'bottom'
variable = u
[../]

[./u_right]
type = SideAverageValue
boundary = 'right'
variable = u
[../]

[./volume]
type = VolumePostprocessor
execute_on = 'timestep_begin initial'
[../]

[./u_total]
type = ElementIntegralVariablePostprocessor
variable = u
[../]

[]

[Executioner]

type = Transient
scheme = implicit-euler

nl_rel_tol = 1e-06
picard_abs_tol = 1e-6
nl_abs_tol = 1e-6
nl_rel_step_tol = 1e-6
picard_rel_tol = 1e-6
nl_abs_step_tol = 1e-6
l_tol = 0.01
l_max_its = 100

solve_type = PJFNK
start_time = 0.0
end_time = 1.0
petsc_options_iname = '-pc_type -pc_hypre_type'
petsc_options_value = 'hypre boomeramg'

[./TimeStepper]
type = ConstantDT
dt = 0.1
[../]

[]

[Adaptivity]

marker = ef

[./Indicators]
[./u_grad_error]
type = GradientJumpIndicator
variable = u
[../]

[../]

[./Markers]
[./ef]
type = ErrorFractionMarker
indicator = u_grad_error
[../]
[../]

[]

[Outputs]
exodus = true
csv = true
[]
